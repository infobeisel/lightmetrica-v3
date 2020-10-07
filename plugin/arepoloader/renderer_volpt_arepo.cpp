/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#include <lm/core.h>
#include <lm/renderer.h>
#include <lm/scene.h>
#include <lm/film.h>
#include <lm/scheduler.h>
#include <lm/path.h>
#include <lm/timer.h>
#include <lm/stats.h>
#include <lm/light.h>
#include "statstags.h"
#include "arepoloader.h"
#include "lm/accel.h"
#define VOLPT_IMAGE_SAMPLING 0

LM_NAMESPACE_BEGIN(LM_NAMESPACE)



void to_json(lm::Json& j, const LightToCameraRaySegmentCDF& p) {
    j = {
        {"weight" , {
            {"x",p.weight.x},
            {"y",p.weight.y},
            {"z",p.weight.z}
            }
        },
        {"cdfSoFar" , p.cdfSoFar},
        {"localcdf" , p.localcdf},
        {"t" , p.t},
        {"a" , p.a},
        {"b" , p.b}
    };
}


void from_json(const lm::Json& j, LightToCameraRaySegmentCDF& p) {
    //do nothing;
}

class Renderer_VolPT_Base : public Renderer {
protected:
    std::string strategy_;
    Float mis_power_;
    Scene* scene_;
    
    Volume_Arepo* volume_;
    Film* film_;
    int max_verts_;
    int knn_min_k_;
    Float rr_prob_;
    std::optional<unsigned int> seed_;
    Component::Ptr<scheduler::Scheduler> sched_;
    AccelKnn* vrlStorage_;
    AccelKnn* pointLightAccel_;
    long long spp_;
    int num_knn_queries_;
    bool sample_lights_;
    bool sample_vrls_;
    Float knn_min_percent_vrls_;
    Float knn_max_percent_vrls_;
    Float knn_min_percent_points_;
    Float knn_max_percent_points_;


public:
    LM_SERIALIZE_IMPL(ar) {
        ar(scene_, film_, max_verts_, rr_prob_, sched_,vrlStorage_);
    }

    virtual void foreach_underlying(const ComponentVisitor& visit) override {
        comp::visit(visit, scene_);
        comp::visit(visit, film_);
        comp::visit(visit, sched_);
        comp::visit(visit, vrlStorage_);
    }

public:
    virtual void construct(const Json& prop) override {
        
        knn_min_k_ = json::value<int>(prop, "knn_min_k",10);
        knn_min_percent_vrls_ = json::value<Float>(prop, "knn_min_percent_vrls",0.01);
        knn_max_percent_vrls_ = json::value<Float>(prop, "knn_max_percent_vrls",0.1);
        knn_min_percent_points_ = json::value<Float>(prop, "knn_min_percent_points",0.1);
        knn_max_percent_points_ = json::value<Float>(prop, "knn_max_percent_points",0.2);
        
        sample_lights_ = json::value<bool>(prop, "sample_lights",true);
        sample_vrls_= json::value<bool>(prop, "sample_vrls",true);
        num_knn_queries_ = json::value<int>(prop, "num_knn_queries",10);

        strategy_ = json::value<std::string>(prop, "strategy");
        mis_power_ = json::value<Float>(prop, "mis_power",2.0);
        scene_ = json::comp_ref<Scene>(prop, "scene");
        volume_= json::comp_ref<Volume_Arepo>(prop, "volume");
        film_ = json::comp_ref<Film>(prop, "output");
        max_verts_ = json::value<int>(prop, "max_verts");
        seed_ = json::value_or_none<unsigned int>(prop, "seed");
        rr_prob_ = json::value<Float>(prop, "rr_prob", .2_f);
        const auto sched_name = json::value<std::string>(prop, "scheduler");
        spp_ = json::value<long long>(prop, "spp");
        vrlStorage_ = json::comp_ref<AccelKnn>(prop, "vrl_accel");
        pointLightAccel_ = json::comp_ref<AccelKnn>(prop, "pointlight_accel");

        #if VOLPT_IMAGE_SAMPLING
        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spi::" + sched_name, make_loc("scheduler"), prop);
        #else
        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spp::" + sched_name, make_loc("scheduler"), prop);
        #endif

    }
};




    
class Renderer_VolPT_Arepo final : public Renderer_VolPT_Base {

public:
    virtual Json render() const override {
		scene_->require_renderable();

        stats::clearGlobal<lm::stats::SampleIdCacheHits,int,long long>( );
        stats::clearGlobal<lm::stats::SampleIdCacheMisses,int,long long>( );
        stats::clearGlobal<lm::stats::UsedCachedTetra,int,long long>( );
        stats::clearGlobal<lm::stats::UsedNeighborTetra,int,long long>( );
        stats::clearGlobal<lm::stats::ResampleAccel,int,long long>( );
        stats::clearGlobal<lm::stats::TotalTetraTests,int,long long>( );


        //film_->clear();
        const auto size = film_->size();
        timer::ScopedTimer st;



        const auto processed = sched_->run([&](long long pixel_index, long long sample_index, int threadid) {
            LM_KEEP_UNUSED(sample_index);
        
            //TODO implement :
            //    - over in vrls render: adding vrls to an accel knn structure (yes not optimal but ok)
            //    - here: sample vrl buffer to choose vrl segments that will shine wonderfully on this the distance 
            //    - think of pdf that direct light splatting gives ? 0 probabliy because its a point light...
            //    - implement pybinding helper for parsing the stars (takes hours in python)
            auto & tetraIToLightSegments =  lm::stats::getGlobalRefUnsafe<lm::stats::VRL,lm::stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>()[0]; //for the moment everything is stored in vector 0
            int num_vrls = tetraIToLightSegments.size();

            auto & tetraToPointLights = stats::getGlobalRefUnsafe<lm::stats::LightsInTetra,lm::stats::TetraIndex,std::vector<int>>();
            /*for(auto p : tetraToPointLights)  {
                for(auto v : p.second)  
                LM_INFO("tetra {}, light {}",p.first,v);
            }*/

            auto & equiContributions = 
            stats::getRef<stats::EquiContribution,int,std::vector<Vec3>>();
            //equiContributions.clear();
            auto & equipdfs = 
            stats::getRef<stats::EquiEquiPDF,int,std::vector<Float>>();
            //equipdfs.clear();

            auto & regularContributions = 
            stats::getRef<stats::RegularContribution,int,std::vector<Vec3>>();
            //regularContributions.clear();
            auto & regularpdfs = 
            stats::getRef<stats::RegularRegularPDF,int,std::vector<Float>>();
            //regularpdfs.clear();

            auto & regular_of_equi_pdfs = 
            stats::getRef<stats::RegularEquiPDF,int,std::vector<Float>>();
            //regular_of_equi_pdfs.clear();

            auto & equi_of_regular_pdfs = 
            stats::getRef<stats::EquiRegularPDF,int,std::vector<Float>>();
            //equi_of_regular_pdfs.clear();



            int num_pointlights = scene_->num_lights(); 



            // Per-thread random number generator
            thread_local Rng rng(seed_ ? *seed_ + threadid : math::rng_seed());
            
             //LM_INFO("current sample id : {}", std::to_string(stats::get<stats::CurrentSampleId,long long>(0)));
            //
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_0,0.0);
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_1,0.0);
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_1_0,0.0);
            stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_1_1,0.0);


            // ------------------------------------------------------------------------------------

            // Sample window
            Vec4 window(0,0,1,1);
            #if VOLPT_IMAGE_SAMPLING
            LM_UNUSED(pixel_index);
            #else
            {
                const int x = int(pixel_index % size.w);
                const int y = int(pixel_index / size.w);
                const auto dx = 1_f / size.w;
                const auto dy = 1_f / size.h;
                window = { dx * x, dy * y, dx, dy };
            }
            #endif


            

            

               



            std::vector<Vec3> EquiMeasurementContributions;
           // EquiMeasurementContributions.reserve(max_verts_);
            std::vector<Float> EquiPdfOfEquiSmpls;
            //EquiPdfOfEquiSmpls.reserve(max_verts_);
            std::vector<Float> EquiPdfOfRegularSmpls;
            //EquiPdfOfRegularSmpls.reserve(max_verts_);
            
            std::vector<Vec3> RegularMeasurementContributions;
            //RegularMeasurementContributions.reserve(max_verts_);
            std::vector<Float> RegularPdfOfEquiSmpls;
            //RegularPdfOfEquiSmpls.reserve(max_verts_);
            std::vector<Float> RegularPdfOfRegularSmpls;
            //RegularPdfOfRegularSmpls.reserve(max_verts_);
            
            
            //how many contributions are there, per path length?
            std::vector<int> contributionIndex_equi = {};
            std::vector<int> contributionIndex_regular = {};
            contributionIndex_equi.resize(max_verts_);
            contributionIndex_regular.resize(max_verts_);
            
            
            Vec2 raster_pos{};
            
            //prepare sample path routine as a lambda, then call it 2 times, once for each strategy
            std::function<Vec3(std::string)> samplePath = [&] (std::string usestrategy) {
                //LM_INFO("good day good Sir, i am thread {} and i will handle sample {} of pixel {}",threadid,sample_index,pixel_index);

                auto pathPDF = 1.0;
             
                Vec3 contributionRegularStrategy = Vec3(0.0);
                Vec3 contributionEquiStrategy = Vec3(0.0);


                
                for(int i = 0; i < max_verts_;i++) { 
                    contributionIndex_equi[i] = 0; // important : path of length 0. gets index, but is ignored later.
                    contributionIndex_regular[i] = 0;
                }

              
                Vec3 contribution = usestrategy == "delta" ? Vec3(0) : Vec3(1);

                // ------------------------------------------------------------------------------------

                // Sample initial vertex
                const auto sE = path::sample_position(rng, scene_, TransDir::EL);
                const auto sE_comp = path::sample_component(rng, scene_, sE->sp, {});
                auto sp = sE->sp;
                int comp = sE_comp.comp;
                auto throughput = sE->weight * sE_comp.weight;

                // ------------------------------------------------------------------------------------

                // Perform random walk
                Vec3 wi{};

                thread_local KNNResult vrl_knn_res;
                thread_local KNNResult point_knn_res;

                thread_local auto vrl_knns = std::priority_queue<Neighbour, std::vector<Neighbour>>();
                thread_local auto point_knns = std::priority_queue<Neighbour, std::vector<Neighbour>>();

                thread_local auto point_knns_wo_duplicates =
                std::unordered_map<int,Neighbour>();
                thread_local auto vrl_knns_wo_duplicates = 
                std::unordered_map<int,Neighbour>();

                int lastDistanceSampleTetraI = -1; //a guess where the next distance sample search will start

                for (int num_verts = 1; num_verts < max_verts_; num_verts++) {
                    vrl_knns = std::priority_queue<Neighbour, std::vector<Neighbour>>();
                    point_knns = std::priority_queue<Neighbour, std::vector<Neighbour>>();

                    
                    
                  


                    //store the sample id that this thread currently works on 
                    stats::set<stats::CachedSampleId,int,long long>(0,max_verts_ * spp_ * pixel_index + max_verts_ * sample_index + num_verts);
            

                    // Sample a NEE edge
                    #if VOLPT_IMAGE_SAMPLING
                    const auto samplable_by_nee = !path::is_specular_component(scene_, sp, comp);
                    #else
                    const auto samplable_by_nee = num_verts > 1 && !path::is_specular_component(scene_, sp, comp);
                    #endif
                    /*if (samplable_by_nee) [&] {
                        // Sample a light
                        auto sL = path::sample_direct(rng, scene_, sp, TransDir::LE);
                        int key = 0;
                        auto pdfLight = stats::get<stats::LastSampledPDF,int,Float>(key);
                        //LM_INFO("light pdf {} ", pdfLight);
                        if(pdfLight < std::numeric_limits<Float>::epsilon())
                            pdfLight = 1.0; //i dont know what the heck
                            
                        sL->weight *= pdfLight; //remove pdf from measurement
                        //
                        // Transmittance
                        const auto Tr = path::eval_transmittance(rng, scene_, sp, sL->sp);
                        //if (math::is_zero(Tr)) {
                        //    return;
                       // }

                        
                        if (!sL || math::is_zero(Tr)  ) {
                            if(usestrategy == "equiangular") {

                                EquiMeasurementContributions.push_back(contribution);

                                EquiPdfOfEquiSmpls.push_back(0);
                                RegularPdfOfEquiSmpls.push_back(0);
                            
                            }
                            if(usestrategy == "regular") {
                                RegularMeasurementContributions.push_back(contribution);

                                EquiPdfOfRegularSmpls.push_back(0);
                                RegularPdfOfRegularSmpls.push_back(0);
                                    
                            }
                            return;
                        }

                        
                        auto sp_from = sp;
                        auto sp_to = sL->sp;
                        auto area_measure_conversion = 1.0 / glm::abs(glm::dot(sp_from.geom.p, sp_to.geom.p));
                        // LM_INFO("area measure conv {} ", area_measure_conversion);
                        //sL->weight *= area_measure_conversion; //convert measurement to area measure
                        pdfLight *= area_measure_conversion; // convert pdf to area measure


                        //TODO not sure if  - sL->wo OR sL->wo
                        //auto solidangledirectionpdf = 1.0;//path::pdf_direct(scene_, sp_from, sp_to, -sL->wo,false);
                        
                        // Recompute raster position for the primary edge
                        Vec2 rp = raster_pos;
                        if (num_verts == 1) {
                            const auto rp_ = path::raster_position(scene_, -sL->wo);
                            if (!rp_) { return; }
                            rp = *rp_;
                        }

                        

                        // Evaluate BSDF
                        const auto wo = -sL->wo;
                        auto fs = path::eval_contrb_direction(scene_, sp, wi, wo, comp, TransDir::EL, true);

                        const auto& primitive = scene_->node_at(sp.primitive).primitive;
                        //todo directions correct ?
                        auto solidangledirectionpdf = primitive.medium->phase()->pdf_direction(sp.geom, wi, wo);
                        //eliminate pdf from phase again:
                        fs *= solidangledirectionpdf;
                        //LM_INFO("phase pdf {} ", solidangledirectionpdf);
                        //fs *= area_measure_conversion; //convert contribution to vertex area
                        solidangledirectionpdf *= area_measure_conversion; //convert pdf to vertex area


                        //if (math::is_zero(fs)) {
                        //    return;
                       // }

                        //int h = 0;
                        //lm::Float normFac = lm::stats::get<lm::stats::MaxTransmittance,int,lm::Float>(h);


                        // Evaluate and accumulate contribution
                        //leave light contribution as is 
                        //phase does not contain pdf anymore, check!
                        const auto C = throughput * Tr * fs * sL->weight * area_measure_conversion; 
                        contribution = usestrategy == "delta" ? contribution + C : contribution * C;

                        //LM_INFO("pdf light{}",pdfLight);

                        //convert to vertex area measure
                        //solidangledirectionpdf *= glm::dot(wi, wo) / glm::dot(sp.geom)

                        //TODO there must be some more pdf to make this correct
                        if(usestrategy == "equiangular") {

                            EquiMeasurementContributions.push_back(contribution);

                            EquiPdfOfEquiSmpls.push_back(pathPdfEquiStrategyEquiSamples * solidangledirectionpdf * pdfLight);
                            RegularPdfOfEquiSmpls.push_back(pathPdfRegularStrategyEquiSamples* solidangledirectionpdf * pdfLight);
                            
                        }
                        if(usestrategy == "regular") {
                            RegularMeasurementContributions.push_back(contribution);

                            EquiPdfOfRegularSmpls.push_back(pathPdfEquiStrategyRegularSamples* solidangledirectionpdf * pdfLight);
                            RegularPdfOfRegularSmpls.push_back(pathPdfRegularStrategyRegularSamples* solidangledirectionpdf * pdfLight);
                                
                        }
                        if(usestrategy == "delta")
                            contribution /= solidangledirectionpdf * pdfLight;
                        
                        //film_->splat(rp, C);
                    }();*/

                    // --------------------------------------------------------------------------------

                    // Sample direction
                    auto s = [&]() -> std::optional<path::DirectionSample> {
                        if (num_verts == 1) {
                            const auto [x, y, w, h] = window.data.data;
                            const auto ud = Vec2(x + w * rng.u(), y + h * rng.u());
                            return path::sample_direction({ ud, rng.next<Vec2>() }, scene_, sp, wi, comp, TransDir::EL);
                        }
                        else {
                            return path::sample_direction(rng, scene_, sp, wi, comp, TransDir::EL);
                        }
                    }();
                    
                    //if (!s) {
                    //    return contribution;
                    //}

                    auto directionpdf = 0.0;
                    {//remember direction pdf
                        if (s) {
                        
                            const auto& primitive = scene_->node_at(sp.primitive).primitive;
                            if(primitive.medium)
                                directionpdf = primitive.medium->phase()->pdf_direction(sp.geom, wi, s->wo);
                            else //camera 
                                directionpdf = 1.0;
                            
                            //eliminate pdf from s again
                            s->weight *= directionpdf;

                            
                        }

                        
                    }
                    // --------------------------------------------------------------------------------

                    // Compute and cache raster position
                    if (num_verts == 1) {
                        raster_pos = *path::raster_position(scene_, s->wo);
                    }

                    //if(num_verts > 1) { //not the camera vertex

                    // --------------------------------------------------------------------------------

                    //hash map to save nearest lights
                    //auto nearestLights = std::unordered_map<int,Float>();
                    //nearestLights.reserve(scene_->num_lights() * 0.1);
                    Float nearestD = std::numeric_limits<Float>::max();
                    int nearestI = 0;
                    std::function<void(int,Float)> nearestLightsVisitor = [&] (int nodeindex,Float distance) -> void {
                        //if(nearestLights[nodeindex] > distance)
                        //    nearestLights[nodeindex] = distance;
                        if(distance < nearestD)
                            nearestI = nodeindex;
                    };
                    stats::set<stats::LightKnnVisitor,int,std::function<void(int,Float)>*>(0,&nearestLightsVisitor);

                    //save nearest lights along current ray 
                    //std::function<void(Vec3,RaySegmentCDF const &)> raysegmentVisitor = [&] (lm::Vec3 boundarypos,lm::RaySegmentCDF const & tetrasegment) -> void {
                    //do nothing for the moment//    scene_->sample_light_selection_from_pos(0.0,boundarypos);
                    //};
                    
                    //stats::set<stats::BoundaryVisitor,int,std::function<void(Vec3,RaySegmentCDF const &)>>(0,raysegmentVisitor);


                    //ask nearest light to current vertex
                    scene_->sample_light_selection_from_pos(0.0,sp.geom.p);

                    lm::Scene::LightPrimitiveIndex lightprim; 

                    // Sample next scene interaction
                    
                    //importance sample distance following light source, strategy 0

                     //do we have vrls in this tetra ? if yes, sample them
                    
                    auto volumeTMin = 0.0;
                    auto volumeTMax = std::numeric_limits<Float>::max();
                    //no not ray dependent, need global! otherwise i am sacked
                    //auto isectret  = volume_->bound().isect_range({sp.geom.p,s->wo}, volumeTMin, volumeTMax);
                    volumeTMax = glm::max(10000.0, glm::distance(volume_->bound().min,volume_->bound().max));
                    
                    

                    //prepare queries
                    auto pdf_vrl_selection = 1.0;
                    auto pdf_light_selection = 1.0;

#ifdef USE_KNN_EMBREE

                    if(num_vrls > 0) {
                        auto min_percent = glm::min(
                            static_cast<Float>(num_vrls) * knn_min_percent_vrls_,static_cast<Float>(num_vrls));
                        auto max_percent = glm::min(
                            static_cast<Float>(num_vrls) * knn_max_percent_vrls_,static_cast<Float>(num_vrls));
                        min_percent = glm::max(min_percent,static_cast<Float>(knn_min_k_));
                        auto vrl_normFac = 1.0 - glm::exp(-(max_percent-min_percent));
                        auto logu = min_percent - glm::log(1 - rng.u() * vrl_normFac);
                        vrl_knn_res.k = static_cast<unsigned int>(logu);
                        vrl_knn_res.k++;

                        //vrl_knn_res.k = num_vrls; //TESST

                        vrl_knn_res.k = vrl_knn_res.k >= num_vrls ? num_vrls : vrl_knn_res.k;
                        vrl_knn_res.k = vrl_knn_res.k < 1 ? 1 : vrl_knn_res.k;
                        //auto pdf_vrl_selection = static_cast<Float>(vrl_knn_res.k) /  static_cast<Float>(num_vrls);
                        auto pdf_vrl_selection = glm::exp(min_percent - vrl_knn_res.k) / vrl_normFac;
  
                    }


                    if(num_pointlights > 0) {
                        auto min_percent = glm::min(static_cast<Float>(num_pointlights) * knn_min_percent_points_,
                            static_cast<Float>(num_pointlights));
                        min_percent = glm::max(min_percent,static_cast<Float>(knn_min_k_));

                        auto max_percent = glm::min(
                            static_cast<Float>(num_pointlights) * knn_max_percent_points_,static_cast<Float>(num_pointlights));

                        
                        //auto point_normFac = 1.0 - glm::exp(-(num_pointlights-min_percent));
                        auto point_normFac = 1.0 - glm::exp(-(max_percent-min_percent));
                        auto logu = min_percent - glm::log(1 - rng.u() * point_normFac);
                        point_knn_res.k = static_cast<unsigned int>(logu);
                        point_knn_res.k++;
                        point_knn_res.k = point_knn_res.k >= num_pointlights ? num_pointlights : point_knn_res.k;
                        point_knn_res.k = point_knn_res.k < 1 ? 1 : point_knn_res.k;

                        //TEST
                        //point_knn_res.k = num_pointlights;

                        auto pdf_light_selection = glm::exp(min_percent - point_knn_res.k) / point_normFac;
                        // static_cast<Float>(point_knn_res.k) /  static_cast<Float>(num_pointlights);

                        
                    }
#endif
                                       
                    //make sure it is nullptr before, so i can test later if we hit sth at all
                    stats::set<stats::LastBoundarySequence,int,std::vector<RaySegmentCDF>*>(0,nullptr);

                    stats::set<stats::TetraIdGuess,int,int>(0,lastDistanceSampleTetraI); //the next distance sampling will start where the last ended!
                
                    //importance sample distance following volume 
                    std::optional<path::DistanceSample> sd = path::sample_distance(rng, scene_, sp, s->wo);

                    

                    //now the knn res contain up to date information

                    int k = 0;
                    std::vector<RaySegmentCDF> * boundaries =  stats::get<stats::LastBoundarySequence,int,std::vector<RaySegmentCDF>*>(k);

                    lastDistanceSampleTetraI = stats::get<stats::RegularTrackingStrategyTetraIndex,int,int>(k);  

                    auto regularT = stats::get<stats::RegularTrackingStrategyDistanceSample,int,Float>(k);  
                    auto regularXi = stats::get<stats::RegularTrackingStrategyXi,int,Float>(k);//sample that was used for regular distance sample
                    auto totalT = stats::get<stats::RegularTrackingStrategyTotalT,int,Float>(k);//sample that was used for regular distance sample
                    auto totalEffT = stats::get<stats::RegularTrackingStrategyTotalEffT,int,Float>(k);
                    auto minT = stats::get<stats::RegularTrackingStrategyMinT,int,Float>(k);//sample that was used for regular distance sample
                    auto totalTau = stats::get<stats::RegularTrackingStrategyTotalTau,int,Float>(k);//sample that was used for regular distance sample
                    auto tauUntilRegularT = stats::get<stats::RegularTrackingStrategyTauUntilScatter,int,Float>(k);//sample that was used for regular distance sample
                    auto segmentCount = stats::get<stats::LastBoundarySequence,int,int>(k);  
                    auto lowDensityNormalizationFactor = stats::get<stats::RegularTrackingStrategyNormFac,int,Float>(k); 
                    auto regularMuT =  stats::get<stats::RegularTrackingStrategyMuT,int,Float>(k);  

                    auto _strat_smpl_key = stats::IJ::_1_1;
                    auto regularStratRegularSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);      

                    Vec3 currentContribution = Vec3(0); 
                    int contribCount = 0;
                    Float currentPdf = 1.0;
                    //bool nee = false;

                    //std::function<void(Vec3,RaySegmentCDF const &, int,Float, Float )> raysegmentVisitor = [&] (lm::Vec3 boundarypos,lm::RaySegmentCDF const & tetrasegment, int tetraI, Float currentTransmittanceDistanceSample, Float maxT) -> void {
                    auto a_d = s->wo;
                    auto a = sp.geom.p;// + a_d * travelT; 
                    //tetraIToLightSegments
                    //if(!nee && currentTransmittanceDistanceSample > travelT + tetrasegment.t * 0.01) {
                    //if( currentTransmittanceDistanceSample > travelT ) {
                        // nee = true;

                    if(boundaries != nullptr && totalTau > 0.0) { //the chance to have in-scattering                        
                        auto & cameraSegments = *boundaries;

                        std::vector<Float> queryTs;
                        std::vector<int> queryTetraInds;
                        {
                            //CHOOSE POINTS FOR KNN QUERIES ALONG RAY, BOTH FOR KNN EMBREE AND DELAUNAY BFS
                            int found_sample = 0;
                            std::vector<Float> zetas;

                            for(int i = 0; i < num_knn_queries_; i++) {
                                zetas.push_back(rng.u() * totalTau);
                                queryTs.push_back(0.0);
                                queryTetraInds.push_back(0);
                            }
                            //TODO PDFS ?!
                            
                            //auto zeta = rng.u() * totalTau; //a new sample within total 
                            auto zetaRegularPDF = 0.0;
                            auto zetaAccCdf = 0.0;
                            //warp Zeta according optical thickness
                            {
                                auto travelT = 0.0;
                                auto segmentThroughput = 1.0;
                                for(int segmentI = 0; segmentI < segmentCount && found_sample < num_knn_queries_; segmentI++) {
                                    auto & tetrasegment =  cameraSegments[segmentI];
                                    for(int i = 0; i < num_knn_queries_; i++) {
                                        if (zetaAccCdf  + tetrasegment.localcdf  > zetas[i]) {
                                            auto normcdf =  zetaAccCdf ;
                                            lm::Float t = sampleCachedICDF_andCDF( zetas[i],zetas[i] , tetrasegment.t ,
                                            normcdf ,  tetrasegment.a ,   tetrasegment.b );
                                            //normcdf = sampleCDF(t,tetrasegment.a,tetrasegment.b );
                                            //accumulate non normalized!
                                            //retAcc = accCdf + normcdf;
                                            //auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                            //auto particle_density = tetrasegment.b + tetrasegment.a * t;
                                            //auto mu_a = crosssection * particle_density;
                                            //auto phase_integrated = 1.0;//isotrope
                                            //auto mu_s = phase_integrated* particle_density;
                                            //auto mu_t = mu_a + mu_s;
                                            //auto scattering_albedo = mu_s / mu_t; 
                                            //contribution = mu_s * transmittance;
                                            //zetaRegularPDF = mu_t / tauUntilRegularT;  
                                            //zetaRegularPDF = mu_t / totalTau;  
                                            //zetaRegularPDF = 1.0;  
                                            queryTs[i] = travelT + t;   
                                            queryTetraInds[i] = tetrasegment.tetraI;
                                            //zetaAccCdf += normcdf;
                                            //break;
                                            found_sample++;
                                        }
                                    }
                                    travelT += tetrasegment.t;
                                    zetaAccCdf += tetrasegment.localcdf;
                                }
                            }

#ifdef USE_KNN_EMBREE
                            for(int i = 0; i < num_knn_queries_; i++) {
                                auto queryPos = a + a_d * queryTs[i]; // point based knn is potentially wrong, need ray based knn!!!

                                vrl_knn_res.knn = std::priority_queue<Neighbour, std::vector<Neighbour>>();
                                //vrl_knn_res.visited.clear();

                                point_knn_res.knn = std::priority_queue<Neighbour, std::vector<Neighbour>>();
                                //point_knn_res.visited.clear();

                                stats::set<stats::KNNLineComp,int,std::pair<Vec3,Vec3>>(0, 
                                    std::make_pair(sp.geom.p, sp.geom.p + s->wo * totalT ));
                                
                                if(num_vrls > 0 && sample_vrls_)
                                    vrlStorage_->queryKnn(queryPos.x,queryPos.y,queryPos.z,
                                        std::numeric_limits<Float>::max(), vrl_knn_res );
                                //LM_INFO("found {} vrls", vrl_knn_res.knn.size());
                                while (!vrl_knn_res.knn.empty())
                                {
                                    vrl_knns.push(vrl_knn_res.knn.top());
                                    vrl_knn_res.knn.pop();
                                    if(vrl_knns.size() > vrl_knn_res.k) {
                                        vrl_knns.pop();
                                    }
                                }
                                if(sample_lights_)
                                    pointLightAccel_->queryKnn(queryPos.x,queryPos.y,queryPos.z,
                                        std::numeric_limits<Float>::max(), point_knn_res );

                                while (!point_knn_res.knn.empty())
                                {
                                    point_knns.push(point_knn_res.knn.top());
                                    point_knn_res.knn.pop();
                                    if(point_knns.size() > point_knn_res.k) {
                                        point_knns.pop();
                                    }
                                }
                            }

                            //remove doubles
                            point_knns_wo_duplicates.clear();
                            while(!point_knns.empty()) {
                                auto & n = point_knns.top();
                                point_knns_wo_duplicates[n.nodeIndex] = n;
                                point_knns.pop();
                            }
                            vrl_knns_wo_duplicates.clear();
                            while(!vrl_knns.empty()) {
                                auto & n = vrl_knns.top();
                                vrl_knns_wo_duplicates[n.nodeIndex] = n;
                                vrl_knns.pop();
                            }


                            //now have w nearest neighbours at hand!
                            //pdf_light_selection = static_cast<Float>(point_knns_wo_duplicates.size()) /  static_cast<Float>(num_pointlights);
                            

#endif
                        }



                        //auto a = sp.geom.p + a_d * currentTransmittanceDistanceSample; //dont use camera ray but camera point
                        auto maxAllowedT = regularT;//tetrasegment.t;//glm::min(
                                //  currentTransmittanceDistanceSample - travelT,
                                    //  tetrasegment.t);

                        
                        
                               // LM_IlseNFO("accel k {}",res.knn.size());
                        //for(auto & vrls : tetraIToLightSegments){
                            //Float vrlCount = vrls.second.size();
                            //for( auto & vrl : vrls.second) {
                            //LM_INFO("have found {} vrls ", vrl_knns.size());
                            //LM_INFO("there are {} vrls ", vrl_knns_wo_duplicates.size() );
                            if(sample_vrls_)
#ifdef USE_KNN_EMBREE
                            for(auto & p : vrl_knns_wo_duplicates) {
                                auto & neighbor = p.second;//point_knns.top();
                            //while(false && num_vrls > 0 && !vrl_knns.empty() ) {

                                //auto & neighbor = vrl_knns.top();
                                //d_total += neighbor.d;
                                //sorted.push_back(neighbor);
                                //LM_INFO("visit primid {}",vrl_index);
                                auto & vrl = tetraIToLightSegments[neighbor.nodeIndex];
                                //vrl_knns.pop();
#else 
                            for( auto & vrl : tetraIToLightSegments) {

#endif


                                auto b_d = vrl.d;
                                auto b = vrl.p;
                                
                                auto n = glm::cross(a_d,b_d);
                                auto n2 = glm::cross(b_d,n);
                                auto a_t = glm::dot((b - a) , n2) / glm::dot(a_d,n2);

                                a_t = glm::max(minT,glm::min(a_t,totalT));

                                auto c1 = a + a_d *  a_t  ;                           
                                auto b_t = glm::dot((a - b) , n) / glm::dot(b_d,n);

                                b_t = glm::max(0.0,glm::min(b_t,vrl.t));

                                auto c2 = b + b_d * b_t;

                                auto h = glm::distance(c1,c2);
                                auto costheta = glm::dot(a_d,b_d);

                                Float sinTheta = glm::sqrt(1.0 - costheta*costheta);

                                Float vrl_scatter_pdf = 0.0;
                                Vec3 vrl_vlp;
                                Float vrl_vlp_t = 0.0;
                                //from vrl clustering paper mitsuba
                                if (sinTheta < std::numeric_limits<Float>::epsilon()) { // VRL and eye ray are (nearly) parallel
                                    // sample uniformly over the VRL
                                    vrl_vlp_t =  rng.u() * vrl.t;
                                    vrl_vlp = b + b_d * vrl_vlp_t;
                                    vrl_scatter_pdf =  1.0 / vrl.t;
                                } else {

                                    Float V0c = - b_t;
                                    Float V1c = vrl.t - b_t;

                                    auto A_end =  glm::asinh((V1c / h) * sinTheta);
                                    auto A_start =  glm::asinh((V0c / h) * sinTheta);

                                    Float newV = h * sinh(A_start + (rng.u() * (A_end - A_start)));
                                    newV = newV / sinTheta;

                                    vrl_scatter_pdf = 1.0f / sqrt(h*h + newV*newV*sinTheta*sinTheta);
                                    Float denom = (A_end - A_start) / sinTheta;
                                    vrl_scatter_pdf = vrl_scatter_pdf/denom;

                                    newV += b_t; //reparametrisation
                                    vrl_vlp_t = newV;
                                    vrl_vlp = b + newV * b_d;


                                }
                                auto vrl_area_measure_conv = 1.0/ (vrl.tSoFar + vrl_vlp_t) / (vrl.tSoFar + vrl_vlp_t);
                                //TODO whcih one ?!
                                //auto vrl_area_measure_conv = 1.0/ ( vrl_vlp_t) / ( vrl_vlp_t);
                                vrl_scatter_pdf *= vrl_area_measure_conv;

                                auto regularPDF_of_equiSample_vrl = 0.0;
                                /*{
                                    auto t = vrl_vlp_t;
                                    auto cdf = vrl.cdfSoFar + 0.5 * vrl.a * t * t + vrl.b * t;
                                    auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                    auto particle_density = vrl.b + vrl.a * t;
                                    auto mu_a = crosssection * particle_density;
                                    auto phase_integrated = 1.0;//isotrope
                                    auto mu_s = phase_integrated* particle_density;
                                    auto mu_t = mu_a + mu_s;
                                    //auto scattering_albedo = mu_s / mu_t; 
                                    regularPDF_of_equiSample_vrl = mu_t * glm::exp(-cdf) ;/// normFac; 
                                }*/

                                //have sampled virtual ray light, now sample point on camera ray
                                auto th = glm::dot((vrl_vlp - a), a_d); 
                                //th = glm::max(0.0, glm::min(maxAllowedT,th ));
                                auto shortest = a + a_d * th - vrl_vlp;
                                auto equih = glm::length(shortest);
                                //what if it the shortest d is not within the line segment? TODO
                                Float a_ = minT-th;//- th;//ray.o - lightPos;
                                Float b_ = totalT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                //Float b_ = maxT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                

                                // auto xi = rng_u_i_equi;

                               //auto xi = rng.u();



                                auto zetaRegularPDF = 1.0;
                                //auto zeta = rng.u();

                                /*auto zeta = rng.u() * lowDensityNormalizationFactor; //a new sample with same normFac
                                lm::Float logzeta = -gsl_log1p(-zeta);
                                auto zetaTransmittance = 1.0;
                                auto zetaT = 0.0;
                                //warp Zeta according transmittance
                                {
                                    auto travelT = 0.0;
                                    auto accCdf = 0.0;
                                    auto segmentThroughput = 1.0;
                                    for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                        auto & tetrasegment =  cameraSegments[segmentI];
                                        if (accCdf  + tetrasegment.localcdf  > logzeta) {
                                            auto normcdf =  accCdf ;
                                            lm::Float t = sampleCachedICDF_andCDF( logzeta,zeta , tetrasegment.t ,
                                            normcdf ,  tetrasegment.a ,   tetrasegment.b );
                                            normcdf = sampleCDF(t,tetrasegment.a,tetrasegment.b );
                                            zetaTransmittance *= glm::exp(-  normcdf );
                                            //accumulate non normalized!
                                            //retAcc = accCdf + normcdf;
                                            auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                            auto particle_density = tetrasegment.b + tetrasegment.a * t;
                                            auto mu_a = crosssection * particle_density;
                                            auto phase_integrated = 1.0;//isotrope
                                            auto mu_s = phase_integrated* particle_density;
                                            auto mu_t = mu_a + mu_s;
                                            //auto scattering_albedo = mu_s / mu_t; 
                                            //contribution = mu_s * transmittance;
                                            zetaRegularPDF = mu_t * zetaTransmittance / lowDensityNormalizationFactor;  
                                            zetaT = travelT + t;   
                                            break;
                                        }
                                        travelT += tetrasegment.t;
                                        accCdf += tetrasegment.localcdf;
                                        zetaTransmittance *= glm::exp(-  accCdf );
                                    }
                                }


                                auto zeta = rng.u() * totalTau; //a new sample within total 
                                auto zetaTransmittance = 1.0;
                                auto zetaT = totalT * zeta / totalTau; //will get replaced in loop
                                auto zetaAccCdf = 0.0;
                                //warp Zeta according optical thickness
                                {
                                    auto travelT = 0.0;
                                    auto segmentThroughput = 1.0;
                                    for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                        auto & tetrasegment =  cameraSegments[segmentI];
                                        if (zetaAccCdf  + tetrasegment.localcdf  > zeta) {
                                            auto normcdf =  zetaAccCdf ;
                                            lm::Float t = sampleCachedICDF_andCDF( zeta,zeta , tetrasegment.t ,
                                            normcdf ,  tetrasegment.a ,   tetrasegment.b );
                                            normcdf = sampleCDF(t,tetrasegment.a,tetrasegment.b );
                                            zetaTransmittance *= glm::exp(-  normcdf );
                                            //accumulate non normalized!
                                            //retAcc = accCdf + normcdf;
                                            auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                            auto particle_density = tetrasegment.b + tetrasegment.a * t;
                                            auto mu_a = crosssection * particle_density;
                                            auto phase_integrated = 1.0;//isotrope
                                            auto mu_s = phase_integrated* particle_density;
                                            auto mu_t = mu_a + mu_s;
                                            //auto scattering_albedo = mu_s / mu_t; 
                                            //contribution = mu_s * transmittance;
                                            //zetaRegularPDF = mu_t / tauUntilRegularT;  
                                            zetaRegularPDF = mu_t / totalTau;  
                                            //zetaRegularPDF = 1.0;  
                                            zetaT = travelT + t;   
                                            zetaAccCdf += normcdf;
                                            break;
                                        }
                                        travelT += tetrasegment.t;
                                        zetaAccCdf += tetrasegment.localcdf;
                                        zetaTransmittance *= glm::exp(-  zetaAccCdf );
                                    }
                                }*/

                                auto zeta = rng.u();// * totalEffT; //a new sample within total 
                                auto zetaTransmittance = 1.0;
                                auto zetaT = totalT * zeta / totalTau; //will get replaced in loop
                                auto zetaAccCdf = 0.0;
                                //warp Zeta according nonzero particle density
                                {
                                    auto travelT = 0.0;
                                    auto effTravelT = 0.0;
                                    auto segmentThroughput = 1.0;
                                    for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                        auto & tetrasegment =  cameraSegments[segmentI];
                                        auto tToAdd = tetrasegment.localcdf > std::numeric_limits<Float>::epsilon() 
                                        ? tetrasegment.t : 0.0;

                                        if (effTravelT + tToAdd >  zeta) {
                                            auto lastBit = zeta - effTravelT;
                                            auto normcdf = sampleCDF(lastBit,tetrasegment.a,tetrasegment.b );
                                            zetaTransmittance *= glm::exp(-  normcdf );
                                            auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                            auto particle_density = tetrasegment.b + tetrasegment.a * lastBit;
                                            auto mu_a = crosssection * particle_density;
                                            auto phase_integrated = 1.0;//isotrope
                                            auto mu_s = phase_integrated* particle_density;
                                            auto mu_t = mu_a + mu_s;
                                            //zetaRegularPDF = (totalT-minT)/ totalEffT;  //debatable TODO
                                            zetaT = travelT + lastBit;   //debatable TODO 

                                            //a_ = travelT -th; 
                                            //b_ = travelT + tetrasegment.t - th;
                                            break;

                                        }
                                        effTravelT += tToAdd;
                        
                                        travelT += tetrasegment.t;

                                        zetaAccCdf += tetrasegment.localcdf;
                                        zetaTransmittance *= glm::exp(-  zetaAccCdf );
                                    }
                                }

                                //warp zeta
                                //zeta =  (zetaT)/totalT ; //zetaT / regularT; //now is between 0 and 1, transmittance-distributed
                                
                                auto theta_a = glm::atan(a_,equih);
                                auto theta_b = glm::atan(b_,equih);

                                //auto equiT = th + h * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                auto equiT = equih * glm::tan((1.0 - zeta) * theta_a + zeta * theta_b);
                                auto t = th +  equiT ;
                                
                                //auto cam_area_measure_conv = 1.0 / (travelT + t) / (travelT + t) ;
                                auto cam_area_measure_conv = 1.0 / t / t;
                                auto cam_scatter_pdf = equih / (equih*equih+glm::pow(equiT,2.0) *(theta_b-theta_a));
                                cam_scatter_pdf *= cam_area_measure_conv;
                                cam_scatter_pdf *= zetaRegularPDF; //conditional pdf
                                
                                auto camera_sample_point = a + a_d * t;


                                auto connectiondir =  glm::normalize(camera_sample_point-vrl_vlp);
                                auto connection_area_measure = 1.0 / glm::abs(glm::dot(camera_sample_point-vrl_vlp,camera_sample_point-vrl_vlp));

                                auto sp1 = SceneInteraction::make_medium_interaction(
                                            scene_->medium_node(),
                                            PointGeometry::make_degenerated(camera_sample_point)
                                        );

                                auto sp2 = SceneInteraction::make_medium_interaction(
                                            scene_->medium_node(),
                                            PointGeometry::make_degenerated(vrl_vlp)
                                        );

                                
                                //virtual point light emittance
                                auto VRL_C = vrl.weight;

                                {
                                    auto sigma_t = (vrl.b + vrl.a * vrl_vlp_t);
                                    auto tau = vrl.cdfSoFar + 0.5 * vrl.a *  vrl_vlp_t * vrl_vlp_t + vrl.b * vrl_vlp_t;
                                    VRL_C *= Vec3(
                                    A_R_A_V_S * sigma_t * glm::exp(-A_R_A_V_T*tau),
                                    A_G_A_V_S * sigma_t * glm::exp(-A_G_A_V_T*tau),
                                    A_B_A_V_S * sigma_t * glm::exp(-A_B_A_V_T*tau)
                                    );
                                }
                                

                                //camera throughput until sample point: 
                                //transmittance!
                                auto travelT = 0.0;
                                auto accCdf = 0.0;
                                auto segmentThroughput = Vec3(1.0);
                                int tetraIndexLandedIn = -1;
                                for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                    auto & tetrasegment =  cameraSegments[segmentI];
                                    if(t < travelT + tetrasegment.t) { //found sample point
                                        auto lastBit = t - travelT;
                                        auto sigma_t = (tetrasegment.b + tetrasegment.a * lastBit);
                                        auto tau = (accCdf 
                                                + 0.5 * tetrasegment.a * lastBit * lastBit + tetrasegment.b * lastBit);
                                        segmentThroughput = Vec3(
                                        A_R_A_V_S * sigma_t * glm::exp(-A_R_A_V_T*tau),
                                        A_G_A_V_S * sigma_t * glm::exp(-A_G_A_V_T*tau),
                                        A_B_A_V_S * sigma_t * glm::exp(-A_B_A_V_T*tau)
                                        );

                                        tetraIndexLandedIn = tetrasegment.tetraI;
                                        break;
                                    }
                                    travelT += tetrasegment.t;
                                    accCdf += tetrasegment.localcdf;
                                }
                                
                                auto throughputCam = throughput * segmentThroughput *
                                cam_area_measure_conv;

                                

                                stats::set<stats::TetraIdGuess,int,int>(0,tetraIndexLandedIn);
                                //now evaluate transmittance
                                //const auto Tr = path::eval_transmittance(rng, scene_, sp1,sp2);
                                //path::eval_transmittance(rng, scene_, sp1,sp2);

                                //channelwise, receive result
                                auto accCDF = 0.0;
                                int k = 0;
                                accCDF = stats::get<stats::OpticalThickness,int,Float>(k);

                                throughputCam.r *= A_R_A_V_S * glm::exp(-accCDF * A_R_A_V_T);
                                throughputCam.g *= A_G_A_V_S * glm::exp(-accCDF * A_G_A_V_T);
                                throughputCam.b *= A_B_A_V_S * glm::exp(-accCDF * A_B_A_V_T);


                                // Evaluate BSDF
                                const auto wo = connectiondir;
                                
                                //TODO evaluate what is wo and what is wi, is it EL or LE ?! 
                                auto fs1 = path::eval_contrb_direction(scene_, sp1, a_d, wo, comp, TransDir::EL, true);
                                //fs1 *= tetrasegment.b + tetrasegment.a * t; //phase times \mu_s
                                auto fs2 = path::eval_contrb_direction(scene_, sp2, b_d, -wo, comp, TransDir::LE, true);
                                //fs2 *= vrl.b + vrl.a * vrl_vlp_t;//phase times \mu_s
                                const auto& primitive = scene_->node_at(scene_->medium_node()).primitive;
                                //todo directions correct ?
                                auto solidangledirectionpdf1 = primitive.medium->phase()->pdf_direction(sp1.geom, a_d, wo);
                                auto solidangledirectionpdf2 = primitive.medium->phase()->pdf_direction(sp2.geom, b_d, wo);
                                fs1 *= solidangledirectionpdf1;
                                fs2 *= solidangledirectionpdf2;
                                solidangledirectionpdf1 *= cam_area_measure_conv; //convert pdf to vertex area
                                solidangledirectionpdf2 *= vrl_area_measure_conv; //convert pdf to vertex area
                                
                                
                                auto segment_contribution = throughputCam  * connection_area_measure * fs1 * fs2 * VRL_C;
                                auto segment_pdf = 
                                pathPDF 
                                * vrl_scatter_pdf
                                * cam_scatter_pdf
                                * pdf_vrl_selection
                                // / vrlCount
                                //* solidangledirectionpdf1
                                //* solidangledirectionpdf2
                                // connection probability stated to be one * connection_area_measure
                                ;
                                
                                //auto otherPdf = 
                                //    pathPdfRegularStrategyEquiSamples
                                //*   regularPDF_of_equiSample_cam
                                //*   regularPDF_of_equiSample_vrl;


                                


                                if(! std::isnan(segment_pdf) && !std::isinf(segment_pdf) && segment_pdf > std::numeric_limits<Float>::epsilon()) {
                                    
                                    if(contributionIndex_equi[num_verts] >= equiContributions[num_verts].size()) {
                                        equiContributions[num_verts].
                                        push_back(segment_contribution );
                                        equipdfs[num_verts].
                                        push_back( segment_pdf  );
                                    } else {
                                        equiContributions[num_verts][contributionIndex_equi[num_verts]] = segment_contribution;
                                        equipdfs[num_verts][contributionIndex_equi[num_verts]] = segment_pdf ;
                                    }
                                    


                                    currentContribution += segment_contribution;
                                    currentPdf *= segment_pdf;  

                                    contributionIndex_equi[num_verts]++;
                                    //contribCount++;

                                    //EquiMeasurementContributions.push_back(segment_contribution);
                                    //EquiPdfOfEquiSmpls.push_back(segment_pdf);
                                    //RegularPdfOfEquiSmpls.push_back(otherPdf);

                                //EquiPdfOfRegularSmpls.push_back(1);
                                //RegularPdfOfRegularSmpls.push_back(segment_pdf);
                                }
                            

                            

                            }

                            //while(true && num_pointlights > 0 && !point_knns.empty()) {
                            if(sample_lights_)
#ifdef USE_KNN_EMBREE
                            for(auto & p : point_knns_wo_duplicates) {
                                auto & neighbor = p.second;//point_knns.top();
                                //d_total += neighbor.d;
                                //sorted.push_back(neighbor);
                                int point_scene_nodeIndex = neighbor.nodeIndex; 
                                
                                //point_knns.pop();

                                //LM_INFO("visit node index {}",point_scene_nodeIndex);
                                auto & pointNode = scene_->node_at(point_scene_nodeIndex);
                                //LM_INFO("{}, {}, {}",scene_->num_lights(),scene_->num_nodes(),point_scene_nodeIndex);
                                auto point_light_index = scene_->light_index_at(point_scene_nodeIndex);
                                auto lightprim = scene_->light_primitive_index_at(point_light_index); //need to ask for "light primitive", given "light index"
                                assert(pointNode.primitive.light != nullptr);

                                Mat4 global_transform = lightprim.global_transform.M;
                                auto & light = *pointNode.primitive.light;
                                auto possample = rng.next<Light::PositionSampleU>();
                                auto lightPos = light.sample_position(possample, Transform(global_transform)).value().geom.p;    

#else 
                            
                            stats::clear<stats::DuplicateWatchdog,int,int>();
                            int lightsPerQuery = glm::max(1,static_cast<int>(point_knn_res.k / num_knn_queries_));
                            for(int i = 0; i < num_knn_queries_; i++) { 
                                auto queryPos = a + a_d * queryTs[i]; 
                                int visitedLightCount = 0;
                                int acceptedBFSLayer = 0;

                                //LM_INFO("knn query {}, expect k {}", i, lightsPerQuery);
                                //perform bfs search

                                
                                stats::set<stats::TetraIdGuess,int,int>(0,queryTetraInds[i]);

                                volume_->visitBFS(queryPos,[&] (int tetraI, int bfsLayer) -> bool {
                                    bool continueBFS = true;
                                    auto & lightNodeIndicesInThisTetra = tetraToPointLights[tetraI];
                                    //LM_INFO("lights associated {} with visit tetra {}, bfs {}",lightNodeIndicesInThisTetra.size(),tetraI, bfsLayer);
                                    
                                    for(auto nodeI : lightNodeIndicesInThisTetra) {
                                        if(!stats::has<stats::DuplicateWatchdog,int,int>(nodeI)){ //if this light hasnt been handled yet, perform lighting!
                                            //LM_INFO("have not seen{} yet", nodeI);
                                            stats::set<stats::DuplicateWatchdog,int,int>(nodeI,1); 
                                            visitedLightCount++;
                                            //continue bfs, or stop, if
                                            continueBFS = true;
                                                    //visitedLightCount > lightsPerQuery && //we have seen enough lights AND
                                                     //       acceptedBFSLayer != bfsLayer //we did just come to the next layer. always handle a layer entirely.
                                                    //? false : true; 
                                            acceptedBFSLayer = bfsLayer;
                                            auto & pointNode = scene_->node_at(nodeI);
                                            //LM_INFO("{}, {}, {}",scene_->num_lights(),scene_->num_nodes(),point_scene_nodeIndex);
                                            auto point_light_index = scene_->light_index_at(nodeI);
                                            auto lightprim = scene_->light_primitive_index_at(point_light_index); //need to ask for "light primitive", given "light index"
                                            assert(pointNode.primitive.light != nullptr);

                                            Mat4 global_transform = lightprim.global_transform.M;
                                            auto & light = *pointNode.primitive.light;
                                            auto possample = rng.next<Light::PositionSampleU>();
                                            auto lightPos = light.sample_position(possample, Transform(global_transform)).value().geom.p;  
                                            //now, handle the light



                                
                            //scene_->traverse_primitive_nodes([&](const SceneNode& node, Mat4 global_transform) {
                            //    if (node.type != SceneNodeType::Primitive) {
                            //        return;
                            //    }
                            //    if (!node.primitive.light) {
                            //        return;
                            //    }
                            //    auto & light = *node.primitive.light;
                            //    auto lightPos = light.sample_position({{0.0,0.0},0.0}, Transform(global_transform)).value().geom.p;  
#endif



                                
                               

                                auto b = lightPos;

                                auto th = glm::dot((b - a), a_d); 
                                auto shortest = a + a_d * th - b;
                                auto equih = glm::length(shortest);

                                
                               
                                //what if it the shortest d is not within the line segment? TODO
                                //Float a_ = 4.0 - th;//- th;//ray.o - lightPos;
                                Float a_ =  minT - th;//- th;//ray.o - lightPos;
                                
                                //Float b_ = maxAllowedT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                //Float b_ = 10.0 - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                Float b_ =  totalT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                
                                //Float b_ = totalT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                //Float b_ = maxT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                

                                // auto xi = rng_u_i_equi;

                               //auto xi = rng.u();

                               

                                auto cam_scatter_pdf_regular = regularStratRegularSmplPDF;
                                auto segmentThroughput_regular = Vec3(
                                            A_R_A_V_S * regularMuT * glm::exp(-A_R_A_V_T*tauUntilRegularT),
                                            A_G_A_V_S * regularMuT * glm::exp(-A_G_A_V_T*tauUntilRegularT),
                                            A_B_A_V_S * regularMuT * glm::exp(-A_B_A_V_T*tauUntilRegularT)
                                            );

                                auto t_regular = regularT;
                                int tetraILandedIn_regular = lastDistanceSampleTetraI;
                                //auto zeta = rng.u(); //a new sample WITHIN the current distance sample
                                /*auto zeta_regular = rng.u() * lowDensityNormalizationFactor; //a new sample with same normFac
                                lm::Float logzeta_regular = -gsl_log1p(-zeta_regular);
                                auto zetaTransmittance_regular = 1.0;
                                //warp Zeta according transmittance
                                {
                                    auto travelT = 0.0;
                                    auto accCdf = 0.0;
                                    for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                        auto & tetrasegment =  cameraSegments[segmentI];
                                        if (accCdf  + tetrasegment.localcdf  > logzeta_regular) {
                                            auto normcdf =  accCdf ;
                                            lm::Float t = sampleCachedICDF_andCDF( logzeta_regular,zeta_regular , tetrasegment.t ,
                                            normcdf ,  tetrasegment.a ,   tetrasegment.b );
                                            normcdf = sampleCDF(t,tetrasegment.a,tetrasegment.b );
                                            zetaTransmittance_regular *= glm::exp(-  normcdf );
                                            //accumulate non normalized!
                                            //retAcc = accCdf + normcdf;
                                            //auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                            //auto particle_density = tetrasegment.b + tetrasegment.a * t;
                                            //auto mu_a = crosssection * particle_density;
                                            //auto phase_integrated = 1.0;//isotrope
                                            //auto mu_s = phase_integrated* particle_density;
                                            //auto mu_t = mu_a + mu_s;
                                            //take the easy one for the moment
                                            auto mu_t = tetrasegment.b + tetrasegment.a * t;

                                            //auto scattering_albedo = mu_s / mu_t; 
                                            //contribution = mu_s * transmittance;
                                            cam_scatter_pdf_regular = mu_t * zetaTransmittance_regular / lowDensityNormalizationFactor;  
                                            t_regular = travelT + t;   

                                            //also save per channel transmittance
                                            segmentThroughput_regular = Vec3(
                                            A_R_A_V_S * mu_t * glm::exp(-A_R_A_V_T*(normcdf+accCdf)),
                                            A_G_A_V_S * mu_t * glm::exp(-A_G_A_V_T*(normcdf+accCdf)),
                                            A_B_A_V_S * mu_t * glm::exp(-A_B_A_V_T*(normcdf+accCdf))
                                            );

                                            tetraILandedIn_regular = tetrasegment.tetraI;
                                            break;
                                        }
                                        travelT += tetrasegment.t;
                                        accCdf += tetrasegment.localcdf;
                                        zetaTransmittance_regular *= glm::exp(-  accCdf );
                                    }
                                }*/


                                

                                //auto zeta = rng.u() * tauUntilRegularT; //a new sample WITHIN the current distance sample
                                /*auto zeta = rng.u() * totalTau; //a new sample within total 
                                auto zetaTransmittance = 1.0;
                                auto zetaT = 0.0;
                                auto zetaAccCdf = 0.0;
                                {
                                    auto travelT = 0.0;
                                    auto segmentThroughput = 1.0;
                                    for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                        auto & tetrasegment =  cameraSegments[segmentI];
                                        if (zetaAccCdf  + tetrasegment.localcdf  > zeta) {
                                            auto normcdf =  zetaAccCdf ;
                                            lm::Float t = sampleCachedICDF_andCDF( zeta,zeta , tetrasegment.t ,
                                            normcdf ,  tetrasegment.a ,   tetrasegment.b );
                                            normcdf = sampleCDF(t,tetrasegment.a,tetrasegment.b );
                                            zetaTransmittance *= glm::exp(-  normcdf );
                                            //accumulate non normalized!
                                            //retAcc = accCdf + normcdf;
                                            auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                            auto particle_density = tetrasegment.b + tetrasegment.a * t;
                                            auto mu_a = crosssection * particle_density;
                                            auto phase_integrated = 1.0;//isotrope
                                            auto mu_s = phase_integrated* particle_density;
                                            auto mu_t = mu_a + mu_s;
                                            //auto scattering_albedo = mu_s / mu_t; 
                                            //contribution = mu_s * transmittance;
                                            //zetaRegularPDF = mu_t / tauUntilRegularT;  
                                            zetaRegularPDF = mu_t / totalTau;  
                                            //zetaRegularPDF = 1.0;  
                                            zetaT = travelT + t;   
                                            zetaAccCdf += normcdf;
                                            break;
                                        }
                                        travelT += tetrasegment.t;
                                        zetaAccCdf += tetrasegment.localcdf;
                                        zetaTransmittance *= glm::exp(-  zetaAccCdf );
                                    }
                                }*/


                                auto theta_a = glm::atan(a_,equih);
                                auto theta_b = glm::atan(b_,equih);
                                auto zeta_equi = rng.u();
                                //auto equiT = th + h * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                auto equiT = equih * glm::tan((1.0 - zeta_equi) * theta_a + zeta_equi * theta_b);
                                auto t_equi = th +  equiT ;
                                //auto cam_area_measure_conv = 1.0 / (travelT + t) / (travelT + t) ;
                                auto cam_area_measure_conv_equi = 1.0 / t_equi / t_equi;
                                auto cam_scatter_pdf_equi = equih / (equih*equih+glm::pow(equiT,2.0) *(theta_b-theta_a));
                                cam_scatter_pdf_equi *= cam_area_measure_conv_equi;
                                auto camera_sample_point_equi = a + a_d * t_equi;
                                auto connectiondir_equi =  glm::normalize(camera_sample_point_equi-b);
                                auto connection_area_measure_equi = 1.0 / glm::abs(glm::dot(camera_sample_point_equi-b,camera_sample_point_equi-b));
                                
                                /*auto zeta = rng.u();// * totalEffT; //a new sample within total 
                                auto zetaTransmittance = 1.0;
                                auto zetaT = totalT * zeta / totalEffT; //will get replaced in loop
                                auto zetaAccCdf = 0.0;
                                //warp Zeta according nonzero particle density
                                {
                                    auto travelT = 0.0;
                                    auto effTravelT = 0.0;
                                    auto segmentThroughput = 1.0;
                                    for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                        auto & tetrasegment =  cameraSegments[segmentI];
                                        auto tToAdd = tetrasegment.localcdf > std::numeric_limits<Float>::epsilon() 
                                        ? tetrasegment.t : 0.0;

                                        if (effTravelT + tToAdd >  zeta) {
                                            auto lastBit = zeta - effTravelT;
                                            auto normcdf = sampleCDF(lastBit,tetrasegment.a,tetrasegment.b );
                                            zetaTransmittance *= glm::exp(-  normcdf );
                                            auto crosssection = 1.0;//327.0/1000.0; //barn ...? but has to be in transmittance as well ugh
                                            auto particle_density = tetrasegment.b + tetrasegment.a * lastBit;
                                            auto mu_a = crosssection * particle_density;
                                            auto phase_integrated = 1.0;//isotrope
                                            auto mu_s = phase_integrated* particle_density;
                                            auto mu_t = mu_a + mu_s;
                                            zetaRegularPDF = 1.0;// (totalT-minT) / totalEffT;  //debatable TODO
                                            zetaT = travelT + lastBit;   //debatable TODO 
                                            //a_ = travelT -th; 
                                            //b_ = travelT + tetrasegment.t - th; 
                                            break;

                                        }
                                        effTravelT += tToAdd;
                        
                                        travelT += tetrasegment.t;

                                        zetaAccCdf += tetrasegment.localcdf;
                                        zetaTransmittance *= glm::exp(-  zetaAccCdf );
                                    }
                                }*/

                                //warp zeta
                                //zeta =  (zetaT)/totalT ; //zetaT / regularT; //now is between 0 and 1, transmittance-distributed
                                //only equiangular sampling
                                //zeta = rng.u();
                                //zetaRegularPDF = 1.0;

                                //zetaRegularPDF *= totalT / (totalT - minT); //evaluate only fraction of path segment by equiangular sampling!
                                
                                /*auto theta_a = glm::atan(a_,equih);
                                auto theta_b = glm::atan(b_,equih);
                                auto zeta_equi = rng.u();
                                //auto equiT = th + h * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                auto equiT = equih * glm::tan((1.0 - zeta_equi) * theta_a + zeta_equi * theta_b);
                                auto t_equi = th +  equiT ;
                                //auto cam_area_measure_conv = 1.0 / (travelT + t) / (travelT + t) ;
                                auto cam_area_measure_conv_equi = 1.0 / t_equi / t_equi;
                                auto cam_scatter_pdf_equi = equih / (equih*equih+glm::pow(equiT,2.0) *(theta_b-theta_a));
                                cam_scatter_pdf_equi *= cam_area_measure_conv_equi;
                                auto camera_sample_point_equi = a + a_d * t_equi;
                                auto connectiondir_equi =  glm::normalize(camera_sample_point_equi-b);
                                auto connection_area_measure_equi = 1.0 / glm::abs(glm::dot(camera_sample_point_equi-b,camera_sample_point_equi-b));*/



                                auto sp1_equi = SceneInteraction::make_medium_interaction(
                                            scene_->medium_node(),
                                            PointGeometry::make_degenerated(camera_sample_point_equi)
                                        );

                                auto sp2 = SceneInteraction::make_medium_interaction(
                                            scene_->medium_node(),
                                            PointGeometry::make_degenerated(b)
                                        );

                                //REGULAR PDF
                                //auto cam_area_measure_conv = 1.0 / (travelT + t) / (travelT + t) ;
                                auto cam_area_measure_conv_regular = 1.0 / t_regular / t_regular;
                                cam_scatter_pdf_regular *= cam_area_measure_conv_regular;
                                auto camera_sample_point_regular = a + a_d * t_regular;
                                auto connectiondir_regular =  glm::normalize(camera_sample_point_regular-b);
                                auto connection_area_measure_regular = 1.0 / glm::abs(glm::dot(camera_sample_point_regular-b,camera_sample_point_regular-b));
                                auto sp1_regular = SceneInteraction::make_medium_interaction(
                                            scene_->medium_node(),
                                            PointGeometry::make_degenerated(camera_sample_point_regular)
                                        );

                                //EQUI PDF OF REGULAR
                                auto cam_scatter_pdf_equi_of_regular = equih / (equih*equih+glm::pow(t_regular-th,2.0) *(theta_b-theta_a));
                                cam_scatter_pdf_equi_of_regular *= cam_area_measure_conv_regular;


                                
                                //virtual point light emittance


                                //auto pdf = light->pdf_direct()




                                //camera throughput until sample points: 
                                //transmittance!
                                auto travelT = 0.0;
                                auto accCdf = 0.0;
                                auto segmentThroughput_equi = Vec3(1.0);
                                int tetraILandedIn_equi = -1;
                                auto cam_scatter_pdf_regular_of_equi = 1.0;
                                bool foundEqui = false;
                                for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                    auto & tetrasegment =  cameraSegments[segmentI];
                                    if(!foundEqui && t_equi < travelT + tetrasegment.t) { //found sample point
                                        auto lastBit = t_equi - travelT;
                                        auto sigma_t = (tetrasegment.b + tetrasegment.a * lastBit);
                                        auto tau = (accCdf 
                                                + 0.5 * tetrasegment.a * lastBit * lastBit + tetrasegment.b * lastBit);
                                        segmentThroughput_equi = Vec3(
                                        A_R_A_V_S * sigma_t * glm::exp(-A_R_A_V_T*tau),
                                        A_G_A_V_S * sigma_t * glm::exp(-A_G_A_V_T*tau),
                                        A_B_A_V_S * sigma_t * glm::exp(-A_B_A_V_T*tau)
                                        );

                                        //also, save regular pdf , but not with per channel
                                        cam_scatter_pdf_regular_of_equi = sigma_t * glm::exp(-tau) / lowDensityNormalizationFactor;  
                                        cam_scatter_pdf_regular_of_equi *= cam_area_measure_conv_equi;   

                                        tetraILandedIn_equi = tetrasegment.tetraI;
                                        foundEqui = true;
                                        break;
                                    }
                                  

                                    travelT += tetrasegment.t;
                                    accCdf += tetrasegment.localcdf;
                                }
                                
                                auto throughputCam_equi = throughput * segmentThroughput_equi *
                                cam_area_measure_conv_equi;


                                stats::set<stats::TetraIdGuess,int,int>(0,tetraILandedIn_equi);
                                //now evaluate transmittance
                                //const auto Tr = path::eval_transmittance(rng, scene_, sp1,sp2);
                                path::eval_transmittance(rng, scene_, sp1_equi,sp2);

                                //channelwise, receive result
                                auto accCDF = 0.0;
                                int k = 0;
                                accCDF = stats::get<stats::OpticalThickness,int,Float>(k);

                                

                                // Evaluate BSDF
                                const auto wo_equi = glm::normalize(b - camera_sample_point_equi);

                                //TODO überprüfen: alle PDFs beim light sampling
                                
                                //TODO evaluate what is wo and what is wi, is it EL or LE ?! 
                                auto fs1_equi = path::eval_contrb_direction(scene_, sp1_equi, a_d, wo_equi, comp, TransDir::EL, true);
                                //fs1 *= tetrasegment.b + tetrasegment.a * t; //phase times \mu_s
                                //point light has no phase auto fs2 = path::eval_contrb_direction(scene_, sp2, b_d, -wo, comp, TransDir::LE, true);
                                //fs2 *= vrl.b + vrl.a * vrl_vlp_t;//phase times \mu_s
                                const auto& primitive = scene_->node_at(scene_->medium_node()).primitive;
                                //todo directions correct ?
                                auto solidangledirectionpdf1_equi = primitive.medium->phase()->pdf_direction(sp1_equi.geom, a_d, wo_equi);
                                //auto solidangledirectionpdf2 = primitive.medium->phase()->pdf_direction(sp2.geom, b_d, wo);
                                fs1_equi *= solidangledirectionpdf1_equi;
                                //fs2 *= solidangledirectionpdf2;
                                //solidangledirectionpdf1_equi *= cam_area_measure_conv_equi; //convert pdf to vertex area
                                //solidangledirectionpdf2 *= vrl_area_measure_conv; //convert pdf to vertex area
                                
                                auto C_equi =  light.eval(sp2.geom,wo_equi,false);

                                
                                C_equi.r *= A_R_A_V_S * glm::exp(-accCDF * A_R_A_V_T);
                                C_equi.g *= A_G_A_V_S * glm::exp(-accCDF * A_G_A_V_T);
                                C_equi.b *= A_B_A_V_S * glm::exp(-accCDF * A_B_A_V_T);


                                //auto segment_contribution = throughputCam * Tr * connection_area_measure * fs1 * fs2 * VRL_C;
                                auto segment_contribution_equi = throughputCam_equi  * connection_area_measure_equi * fs1_equi * C_equi;
                                auto segment_pdf_equi = 
                                pathPDF 
                                //* vrl_scatter_pdf
                                * cam_scatter_pdf_equi
                                * pdf_light_selection
                                // / vrlCount
                                //* solidangledirectionpdf1
                                //* solidangledirectionpdf2
                                // connection probability stated to be one * connection_area_measure
                                ;
                                auto segment_pdf_equi_of_regular = 
                                pathPDF * cam_scatter_pdf_equi_of_regular * pdf_light_selection;
                                
                                //auto otherPdf = 
                                //    pathPdfRegularStrategyEquiSamples
                                //*   regularPDF_of_equiSample_cam
                                //*   regularPDF_of_equiSample_vrl;




                                auto throughputCam_regular = throughput * segmentThroughput_regular *
                                cam_area_measure_conv_regular;


                                stats::set<stats::TetraIdGuess,int,int>(0,tetraILandedIn_regular);
                                //now evaluate transmittance
                                //const auto Tr = path::eval_transmittance(rng, scene_, sp1,sp2);
                                path::eval_transmittance(rng, scene_, sp1_regular,sp2);

                                //channelwise, receive result
                                accCDF = 0.0;
                                k = 0;
                                accCDF = stats::get<stats::OpticalThickness,int,Float>(k);

                                
                                // Evaluate BSDF
                                const auto wo_regular = glm::normalize(b - camera_sample_point_regular);

                                //TODO überprüfen: alle PDFs beim light sampling
                                
                                //TODO evaluate what is wo and what is wi, is it EL or LE ?! 
                                auto fs1_reqular = path::eval_contrb_direction(scene_, sp1_regular, a_d, wo_regular, comp, TransDir::EL, true);
                                //fs1 *= tetrasegment.b + tetrasegment.a * t; //phase times \mu_s
                                //point light has no phase auto fs2 = path::eval_contrb_direction(scene_, sp2, b_d, -wo, comp, TransDir::LE, true);
                                //fs2 *= vrl.b + vrl.a * vrl_vlp_t;//phase times \mu_s
                                //const auto& primitive = scene_->node_at(scene_->medium_node()).primitive;
                                //todo directions correct ?
                                auto solidangledirectionpdf1_regular = primitive.medium->phase()->pdf_direction(sp1_regular.geom, a_d, wo_regular);
                                //auto solidangledirectionpdf2 = primitive.medium->phase()->pdf_direction(sp2.geom, b_d, wo);
                                fs1_reqular *= solidangledirectionpdf1_regular;
                                //fs2 *= solidangledirectionpdf2;
                                solidangledirectionpdf1_regular *= cam_area_measure_conv_regular; //convert pdf to vertex area
                                //solidangledirectionpdf2 *= vrl_area_measure_conv; //convert pdf to vertex area
                                
                                auto C_regular =  light.eval(sp2.geom,wo_regular,false);

                                C_regular.r *= A_R_A_V_S * glm::exp(-accCDF * A_R_A_V_T);
                                C_regular.g *= A_G_A_V_S * glm::exp(-accCDF * A_G_A_V_T);
                                C_regular.b *= A_B_A_V_S * glm::exp(-accCDF * A_B_A_V_T);

                                //auto segment_contribution = throughputCam * Tr * connection_area_measure * fs1 * fs2 * VRL_C;
                                auto segment_contribution_regular = throughputCam_regular  * connection_area_measure_regular * fs1_reqular * C_regular;
                                auto segment_pdf_regular = 
                                pathPDF //pathPdfEquiStrategyEquiSamples  //yes it is the main path part
                                //* vrl_scatter_pdf
                                * cam_scatter_pdf_regular
                                * pdf_light_selection
                                // / vrlCount
                                //* solidangledirectionpdf1
                                //* solidangledirectionpdf2
                                // connection probability stated to be one * connection_area_measure
                                ;
                                auto segment_pdf_regular_of_equi = 
                                    pathPDF * cam_scatter_pdf_regular_of_equi * pdf_light_selection;
                                //auto otherPdf = 
                                //    pathPdfRegularStrategyEquiSamples
                                //*   regularPDF_of_equiSample_cam
                                //*   regularPDF_of_equiSample_vrl;


                                //invalid, fill with appropriate values
                                if( std::isnan(segment_pdf_equi) || std::isinf(segment_pdf_equi) || segment_pdf_equi < std::numeric_limits<Float>::epsilon()) {
                                    segment_contribution_equi = Vec3(0.0);
                                    segment_pdf_equi = 0.0;
                                    segment_pdf_equi_of_regular = 0.0;
                                }
                                    if(contributionIndex_equi[num_verts] >= equiContributions[num_verts].size()) {
                                        equiContributions[num_verts].
                                        push_back(segment_contribution_equi );
                                        equipdfs[num_verts].
                                        push_back( segment_pdf_equi  );
                                        equi_of_regular_pdfs[num_verts].
                                        push_back(segment_pdf_equi_of_regular);
                                        
                                    } else {
                                        equiContributions[num_verts][contributionIndex_equi[num_verts]] = segment_contribution_equi;
                                        equipdfs[num_verts][contributionIndex_equi[num_verts]] = segment_pdf_equi ;
                                        equi_of_regular_pdfs[num_verts][contributionIndex_equi[num_verts]] = segment_pdf_equi_of_regular;
                                    }
                                    currentContribution += segment_contribution_equi;
                                    currentPdf *= segment_pdf_equi;  
                                    contributionIndex_equi[num_verts]++;
                                




                                if(std::isnan(segment_pdf_regular) || std::isinf(segment_pdf_regular) || segment_pdf_regular < std::numeric_limits<Float>::epsilon()) {
                                    segment_contribution_regular = Vec3(0.0);
                                    segment_pdf_regular = 0.0;
                                    segment_pdf_regular_of_equi = 0.0;
                                }
                                if(contributionIndex_regular[num_verts] >= regularContributions[num_verts].size()) {
                                    regularContributions[num_verts].
                                    push_back(segment_contribution_regular );
                                    regularpdfs[num_verts].
                                    push_back( segment_pdf_regular  );
                                    regular_of_equi_pdfs[num_verts].
                                    push_back(segment_pdf_regular_of_equi);
                                } else {
                                    regularContributions[num_verts][contributionIndex_regular[num_verts]] = segment_contribution_regular;
                                    regularpdfs[num_verts][contributionIndex_regular[num_verts]] = segment_pdf_regular ;
                                    regular_of_equi_pdfs[num_verts][contributionIndex_regular[num_verts]] = segment_pdf_regular_of_equi;
                                }
                                currentContribution += segment_contribution_regular;
                                currentPdf *= segment_pdf_regular;  
                                contributionIndex_regular[num_verts]++;
                                


                                


                            

#ifdef USE_KNN_EMBREE
                            
                            }
#else
                            }//if light 
                            //for each light in tetra
                            }
                            return continueBFS;
                            });//for bfs tetras
                            
                            

                            }
#endif

                       // }

                    }
                    //TODO DEBATE!
                    //throughput = throughput / regularT / regularT;
                    //pathPDF = pathPDF / regularT / regularT;
                    
                    // --------------------------------------------------------------------------------

                    // Update throughput
                    throughput *= s->weight * sd->weight; //todo convert to vertex area measure!
                    
                    // --------------------------------------------------------------------------------

                    // Accumulate contribution from emissive interaction
                    /* doesnt happen here
                    if (!samplable_by_nee && scene_->is_light(sd->sp)) {
                        const auto spL = sd->sp.as_type(SceneInteraction::LightEndpoint);
                        const auto woL = -s->wo;
                        const auto Le = path::eval_contrb_direction(scene_, spL, {}, woL, comp, TransDir::LE, true);
                        const auto C = throughput * Le;
                        film_->splat(raster_pos, C);
                    }*/
                    
                    // --------------------------------------------------------------------------------

                    // Termination on a hit with environment
                    //if (sd->sp.geom.infinite) {
                    //    break;
                   // }

                    // Russian roulette
                    /*if (num_verts > 5) {
                        const auto q = glm::max(rr_prob_, 1_f - glm::compMax(throughput));
                        if (rng.u() < q) {
                            break;
                        }
                        throughput /= 1_f - q;
                    }*/

                    // --------------------------------------------------------------------------------

                    // Sample component
                    const auto s_comp = path::sample_component(rng, scene_, sd->sp, -s->wo);
                    throughput *= s_comp.weight;

                    

                    // --------------------------------------------------------------------------------

                    // Update information
                    wi = -s->wo;
                    sp = sd->sp;
                    comp = s_comp.comp;

                    //} //endif num_vertex > 1 ()
                }

                return contribution;
                
            };


            //constructed 2 paths, one with equiangular, one with regular strategy,
            //combine contributions with MIS

            if(strategy_ == "equiangular") {
                samplePath("equiangular");
                auto contributionEquiStrategy = Vec3(0.0);

                /*auto & cis = stats::getRef<stats::EquiContribution,int,std::vector<Vec3>>();
                for(auto p : cis) {
                    for(auto v : p.second)
                        LM_INFO("contrib path length {}: {},{},{}",p.first,v.x,v.y,v.z );
                }*/


                for(int i = 0; i < max_verts_; i++) { //for paths of all lengths
                    auto & cs = stats::getRef<stats::EquiContribution,int,std::vector<Vec3>>(i);
                    auto & pdfs = stats::getRef<stats::EquiEquiPDF,int,std::vector<Float>>(i);
                    //Float segments = stats::get<stats::EquiContribution,int,int>(i);
                    auto pathsWithLengtI = contributionIndex_equi[i];
                    //cs.size() can be anything from previous samples on this thread...


                    Vec3 ci = Vec3(0);
                    Float js = 0.0;
                    //for(int j = 0; j < cs.size(); j++) { //for pahts with length i
                    for(int j = 0; j < pathsWithLengtI; j++) { //for pahts with length i
                        //LM_INFO("path length {}, sample {} , pdf {} ,  contr {},{},{}", pathsWithLengtI,j,pdfs[j],cs[j].x,cs[j].y,cs[j].z );
                        if(pdfs[j] > std::numeric_limits<Float>::epsilon()) { //a valid sample
                            ci += cs[j] / pdfs[j];
                            js += 1.0;
                        }
                    }
                    if(js != 0.0) {
                        ci /= static_cast<Float>(js); //divided by how many paths
                       //  LM_INFO("divide by: {}",js);
                    }
                    //if(segments != 0.0)
                    //     ci /= segments; //divided by how many segments in camera path

                    contributionEquiStrategy += ci;
                }
                
                //contributionEquiStrategy /= static_cast<Float>(EquiMeasurementContributions.size()); //divided by how many paths
                film_->splat(raster_pos, contributionEquiStrategy / static_cast<Float>(spp_));
                
            }
            else if(strategy_ == "regular") {
                samplePath("regular");
                auto contributionRegularStrategy = Vec3(0.0);

                for(int i = 0; i < max_verts_; i++) { //for paths of all lengths
                    auto & cs = stats::getRef<stats::RegularContribution,int,std::vector<Vec3>>(i);
                    auto & pdfs = stats::getRef<stats::RegularRegularPDF,int,std::vector<Float>>(i);
                    //Float segments = stats::get<stats::EquiContribution,int,int>(i);
                    auto pathsWithLengtI = contributionIndex_regular[i];
                    //cs.size() can be anything from previous samples on this thread...

                    Vec3 ci = Vec3(0);
                    Float js = 0.0;
                    //for(int j = 0; j < cs.size(); j++) { //for pahts with length i
                    for(int j = 0; j < pathsWithLengtI; j++) { //for pahts with length i
                        if(pdfs[j] > std::numeric_limits<Float>::epsilon()) { //a valid sample
                            ci += cs[j] / pdfs[j];
                            js += 1.0;
                        }
                    }
                    if(js != 0.0) {
                        ci /= static_cast<Float>(js); //divided by how many paths
                       //  LM_INFO("divide by: {}",js);
                    }
                    //if(segments != 0.0)
                    //     ci /= segments; //divided by how many segments in camera path

                    contributionRegularStrategy += ci;
                }
                
                //contributionEquiStrategy /= static_cast<Float>(EquiMeasurementContributions.size()); //divided by how many paths
                film_->splat(raster_pos, contributionRegularStrategy / static_cast<Float>(spp_));
            }
            else if(strategy_ == "mis") {
                samplePath("mis");
                //samplePath("regular");

                

                auto totalContribution = Vec3(0.0);
                

                for(int path_length = 0; path_length < max_verts_; path_length++) { //for paths of all lengths

                    auto & reg_cs = stats::getRef<stats::RegularContribution,int,std::vector<Vec3>>(path_length);
                    auto & reg_pdfs = stats::getRef<stats::RegularRegularPDF,int,std::vector<Float>>(path_length);

                    auto & equ_cs = stats::getRef<stats::EquiContribution,int,std::vector<Vec3>>(path_length);
                    auto & equ_pdfs = stats::getRef<stats::EquiEquiPDF,int,std::vector<Float>>(path_length);

                    
                    auto & reg_equi_pdfs = stats::getRef<stats::RegularEquiPDF,int,std::vector<Float>>(path_length);
                    auto & equ_regu_pdfs = stats::getRef<stats::EquiRegularPDF,int,std::vector<Float>>(path_length);

                    std::vector<Float> regular_weights = {};
                    std::vector<Float> equi_weights = {};

                    Vec3 weightedContributionRegular = Vec3(0.0);
                    Vec3 weightedContributionEqui =  Vec3(0.0);

                    auto pathsWithLengtI = contributionIndex_regular[path_length];
                    if(contributionIndex_regular[path_length] != contributionIndex_equi[path_length])
                        LM_ERROR("must be same amount of samples");
                    //regular strategy
                    auto weightSum = 0.0;
                    auto validSamplesRegular = 0;
                    for(int path_i = 0; path_i < pathsWithLengtI; path_i++) { //each sample i
                        if(reg_pdfs[path_i] > std::numeric_limits<Float>::epsilon()) { //a valid sample
                            auto weight = glm::pow(reg_pdfs[path_i],mis_power_) / (glm::pow(reg_pdfs[path_i],mis_power_) + glm::pow(equ_regu_pdfs[path_i],mis_power_));
                            weightedContributionRegular += weight * reg_cs[path_i] / reg_pdfs[path_i];// 
                            validSamplesRegular += 1;
                            regular_weights.push_back(weight);                        
                            weightSum += weight;
                        }
                    }
                    //LM_INFO("{} paths of length {}, weight sum regular {}, valid paths {} ",pathsWithLengtI,path_length, weightSum,validSamplesRegular);

                    //equi strategy
                    //weightSum = 0.0;
                    auto validSamplesEqui = 0;
                    for(int path_i = 0; path_i < pathsWithLengtI; path_i++) { //each sample i
                        if(equ_pdfs[path_i] > std::numeric_limits<Float>::epsilon()) { //a valid sample
                            auto weight = glm::pow(equ_pdfs[path_i],mis_power_) / (glm::pow(equ_pdfs[path_i],mis_power_) + glm::pow(reg_equi_pdfs[path_i],mis_power_));
                            weightedContributionEqui += weight * equ_cs[path_i] / equ_pdfs[path_i];// * ;
                            validSamplesEqui += 1;
                            equi_weights.push_back(weight);                        
                            weightSum += weight;
                        }
                    }
                    if(weightSum < 1.0 - 0.000001 && validSamplesEqui != 0) {
                        //LM_ERROR("weights are not one: {}", weightSum);
                        weightedContributionEqui = Vec3(0.0);
                        weightedContributionRegular =  Vec3(0.0);
                        

                    }
                    //LM_INFO("{} paths of length {}, weight sum equi {}, valid paths {} ",pathsWithLengtI,path_length, weightSum,validSamplesEqui);

                    if(validSamplesEqui != validSamplesRegular)
                        LM_ERROR("do not have same amound of VALID samples: {} equi, {} regular", validSamplesEqui,validSamplesRegular );

                    if(validSamplesEqui != 0)
                        totalContribution += weightedContributionEqui / static_cast<Float>(validSamplesEqui);// + weightedContributionRegular / static_cast<Float>(validSamplesRegular);
                    if(validSamplesRegular != 0)
                        totalContribution += weightedContributionRegular / static_cast<Float>(validSamplesRegular);// + weightedContributionRegular / static_cast<Float>(validSamplesRegular);

                }

                film_->splat(raster_pos, totalContribution/ static_cast<Float>(spp_) );/// static_cast<Float>(spp_));

                
            } else { //onesample_mis or delta strategy 
                auto c = samplePath("delta");
                film_->splat(raster_pos, c/ static_cast<Float>(spp_));
            }


            /*for(int i = 0; i < max_verts_; i++) { //for paths of all lengths
                auto & cs = stats::getRef<stats::EquiContribution,int,std::vector<Vec3>>(i);
                auto & pdfs = stats::getRef<stats::EquiEquiPDF,int,std::vector<Float>>(i);
                //Float segments = stats::get<stats::EquiContribution,int,int>(i);
                cs.clear();
                pdfs.clear();
                cs.shrink_to_fit();
                pdfs.shrink_to_fit();
                
            }
            int k = 0;
            std::vector<RaySegmentCDF> * boundaries =  stats::get<stats::LastBoundarySequence,int,std::vector<RaySegmentCDF>*>(k);
            if(boundaries != nullptr) {
                boundaries->clear();
                boundaries->shrink_to_fit();
            }

            lm::stats::set<int,int,std::vector<lm::RaySegmentCDF>>(0,{}); */


        
            
        },  
        [&](auto pxlindx,auto smplindx,auto threadid) {
            stats::clear<lm::stats::SampleIdCacheHits,int,long long>( );
            stats::clear<lm::stats::SampleIdCacheMisses,int,long long>( );
            stats::clear<lm::stats::UsedCachedTetra,int,long long>( );
            stats::clear<lm::stats::UsedNeighborTetra,int,long long>( );
            stats::clear<lm::stats::ResampleAccel,int,long long>( );
            stats::clear<lm::stats::TotalTetraTests,int,long long>( );

        } , 
        [&](auto pxlindx,auto smplindx,auto threadid) {
            stats::mergeToGlobal<lm::stats::SampleIdCacheHits,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::SampleIdCacheMisses,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::UsedCachedTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::UsedNeighborTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::ResampleAccel,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<lm::stats::TotalTetraTests,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );

        }
        );

        auto smplhits = stats::getGlobal<lm::stats::SampleIdCacheHits,int,long long>(0 );
        auto smplmisses = stats::getGlobal<lm::stats::SampleIdCacheMisses,int,long long>(0 );
        auto tetrahits =  stats::getGlobal<lm::stats::UsedCachedTetra,int,long long>( 0);
        auto tetraneighborhits = stats::getGlobal<lm::stats::UsedNeighborTetra,int,long long>(0 );
        auto accelsmpls = stats::getGlobal<lm::stats::ResampleAccel,int,long long>( 0);
        auto totaltetratests = stats::getGlobal<lm::stats::TotalTetraTests,int,long long>(0 );

        LM_INFO("sample hits: {}, misses : {}, tetra hits {}, tetra neighbor hits {}, accel smpls {} . total tetra probes {}", 
         smplhits,
         smplmisses,
         tetrahits,
         tetraneighborhits,
         accelsmpls,
         totaltetratests
         );

        // Rescale film
        #if VOLPT_IMAGE_SAMPLING
        film_->rescale(Float(size.w* size.h) / processed);
        #else
       // film_->rescale(1_f / processed);
        #endif

        return { {"processed", processed}, {"elapsed", st.now()} };
    }
};

LM_COMP_REG_IMPL(Renderer_VolPT_Arepo, "renderer::volpt_arepo");

LM_NAMESPACE_END(LM_NAMESPACE)

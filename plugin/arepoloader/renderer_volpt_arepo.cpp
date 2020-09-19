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
    Volume* volume_;
    Film* film_;
    int max_verts_;
    Float rr_prob_;
    std::optional<unsigned int> seed_;
    Component::Ptr<scheduler::Scheduler> sched_;
    long long spp_;


public:
    LM_SERIALIZE_IMPL(ar) {
        ar(scene_, film_, max_verts_, rr_prob_, sched_);
    }

    virtual void foreach_underlying(const ComponentVisitor& visit) override {
        comp::visit(visit, scene_);
        comp::visit(visit, film_);
        comp::visit(visit, sched_);
    }

public:
    virtual void construct(const Json& prop) override {
        strategy_ = json::value<std::string>(prop, "strategy");
        mis_power_ = json::value<Float>(prop, "mis_power",2.0);
        scene_ = json::comp_ref<Scene>(prop, "scene");
        volume_= json::comp_ref<Volume>(prop, "volume");
        film_ = json::comp_ref<Film>(prop, "output");
        max_verts_ = json::value<int>(prop, "max_verts");
        seed_ = json::value_or_none<unsigned int>(prop, "seed");
        rr_prob_ = json::value<Float>(prop, "rr_prob", .2_f);
        const auto sched_name = json::value<std::string>(prop, "scheduler");
        spp_ = json::value<long long>(prop, "spp");

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


        film_->clear();
        const auto size = film_->size();
        timer::ScopedTimer st;



        const auto processed = sched_->run([&](long long pixel_index, long long sample_index, int threadid) {
            LM_KEEP_UNUSED(sample_index);

            //TODO implement :
            //    - over in vrls render: adding vrls to a accel knn structure (yes not optimal but ok)
            //    - here: sample vrl buffer to choose vrl segments that will shine wonderfully on this the distance 
            //    - think of pdf that direct light splatting gives ? 0 probabliy because its a point light...

            //stats::getGlobalRefUnsafe<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>(0); 

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


            

            auto pathPdfEquiStrategyEquiSamples = 1.0;
            auto pathPdfEquiStrategyRegularSamples = 1.0;
            auto pathPdfRegularStrategyRegularSamples = 1.0;
            auto pathPdfRegularStrategyEquiSamples = 1.0;


            auto pathPdfDeltaStrategy = 1.0;

            Vec3 contributionRegularStrategy = Vec3(0.0);
            Vec3 contributionEquiStrategy = Vec3(0.0);


            //generate uniform distributed samples here to make the same samples accessible 
            //for both paths

            /*for(int i = 0; i < 2 * max_verts_;i++) { //we have two strategies therefore twice the amount
                stats::set<stats::DistanceSampleRandomValues,int,Float>(i,rng.u());
            }*/

            std::vector<Vec3> EquiMeasurementContributions;
            EquiMeasurementContributions.reserve(max_verts_);
            std::vector<Float> EquiPdfOfEquiSmpls;
            EquiPdfOfEquiSmpls.reserve(max_verts_);
            std::vector<Float> EquiPdfOfRegularSmpls;
            EquiPdfOfRegularSmpls.reserve(max_verts_);
            
            std::vector<Vec3> RegularMeasurementContributions;
            RegularMeasurementContributions.reserve(max_verts_);
            std::vector<Float> RegularPdfOfEquiSmpls;
            RegularPdfOfEquiSmpls.reserve(max_verts_);
            std::vector<Float> RegularPdfOfRegularSmpls;
            RegularPdfOfRegularSmpls.reserve(max_verts_);
            
            
            
            Vec2 raster_pos{};
            //prepare sample path routine as a lambda, then call it 2 times, once for each strategy
            std::function<Vec3(std::string)> samplePath = [&] (std::string usestrategy) {
                int vertexIndexStatsKey = 0;
                stats::set<stats::EquiDistanceSampleRandomValueVertexIndex,int,int>(vertexIndexStatsKey,0);
                stats::set<stats::RegularDistanceSampleRandomValueVertexIndex,int,int>(vertexIndexStatsKey,max_verts_);

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


                for (int num_verts = 1; num_verts < max_verts_; num_verts++) {
                    
                    
                    //update current random sample index
                    stats::set<stats::EquiDistanceSampleRandomValueVertexIndex,int,int>(vertexIndexStatsKey,num_verts);
                    stats::set<stats::RegularDistanceSampleRandomValueVertexIndex,int,int>(vertexIndexStatsKey,num_verts+ max_verts_);

                    //receive current random sample
                                    

                    //auto rng_u_i_equi = stats::get<stats::DistanceSampleRandomValues,int,Float>(num_verts);

                    //auto rng_u_i_regular = stats::get<stats::DistanceSampleRandomValues,int,Float>(num_verts + max_verts_);
                    


                    //store the sample id that this thread currently works on 
                    stats::set<stats::CachedSampleId,int,long long>(0,max_verts_ * spp_ * pixel_index + max_verts_ * sample_index + num_verts);
            

                    // Sample a NEE edge
                    #if VOLPT_IMAGE_SAMPLING
                    const auto samplable_by_nee = !path::is_specular_component(scene_, sp, comp);
                    #else
                    const auto samplable_by_nee = num_verts > 1 && !path::is_specular_component(scene_, sp, comp);
                    #endif
                    if (samplable_by_nee) [&] {
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
                    }();

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
                    std::function<void(Vec3,RaySegmentCDF const &)> raysegmentVisitor = [&] (lm::Vec3 boundarypos,lm::RaySegmentCDF const & tetrasegment) -> void {
                    //do nothing for the moment//    scene_->sample_light_selection_from_pos(0.0,boundarypos);
                    };
                    
                    stats::set<stats::BoundaryVisitor,int,std::function<void(Vec3,RaySegmentCDF const &)>>(0,raysegmentVisitor);

                    //ask nearest light to current vertex
                    scene_->sample_light_selection_from_pos(0.0,sp.geom.p);

                    lm::Scene::LightPrimitiveIndex lightprim; 

                    // Sample next scene interaction
                    
                    //importance sample distance following light source, strategy 0


                    
                    auto lightDistanceSample = path::DistanceSample();
                    lightDistanceSample.weight = Vec3(0);
                    auto equiangularT = 0.0;
                    if(usestrategy != "delta") {


                        lightprim = scene_->light_primitive_index_at(nearestI);


                        auto mediumindex = scene_->medium_node();
                        auto & medium = scene_->node_at(mediumindex);
                        //medium.primitive.medium->
                        Ray ray = {
                        sp.geom.p,
                        s->wo};
                        
                        //sample next scene interaction following light sources situation
                       // scene_->traverse_primitive_nodes(
                        //    [&] (const SceneNode& node, Mat4 global_transform) {
                        const auto & node = scene_->node_at(lightprim.index);
                        Mat4 global_transform = lightprim.global_transform.M;
                                if(node.primitive.light != nullptr && !node.primitive.light->is_infinite()) {
                                    auto & light = *node.primitive.light;
                                    auto possample = rng.next<Light::PositionSampleU>();
                                    auto lightPos = light.sample_position(possample, Transform(global_transform)).value().geom.p;
                                    //LM_INFO("light pos {},{},{}",lightPos.x,lightPos.y,lightPos.z);
                                    auto th = glm::dot((lightPos - ray.o), ray.d); 
                                    auto shortest = ray.o + ray.d * th - lightPos;
                                    auto h = glm::length(shortest);
                                    Float a = -th;//- th;//ray.o - lightPos;
                                    Float b = 1000.9 ;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                    auto theta_a = glm::atan(a,h);
                                    auto theta_b = glm::atan(b,h);

                                   // auto xi = rng_u_i_equi;

                                    auto xi = rng.u();
                                    //auto equiT = th + h * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                    auto equiT = h * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                    auto pdf = h / (h*h+glm::pow(equiT,2.0) *(theta_b-theta_a));
                                    auto t = th +  equiT ;
                                    equiangularT = t;
                                    
                                    //store sample of strategy 0
                                    int key = 0;
                                    stats::set<stats::EquiangularStrategyDistanceSample,int,Float>(key,t);
                                    //store pdf for sample of strategy 0
                                    stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_0,pdf);

                                    auto T = medium.primitive.medium->eval_transmittance(
                                    rng,  ray,0.0, t);
                                    auto mu_s = volume_->eval_scalar(ray.o + ray.d * t);//,ray.d);
                                    //assume mu_s with isotropic phase, i.e.,particle_dens *  1

                                    auto albedo = Vec3( mu_s* T );



                                    lightDistanceSample = path::DistanceSample {
                                        SceneInteraction::make_medium_interaction(
                                            scene_->medium_node(),
                                            PointGeometry::make_degenerated(ray.o + ray.d * t)
                                        ),
                                        albedo
                                    };
                                    
                                }
                       //});
                    }

                    

                    
                    //importance sample distance following volume 
                    std::optional<path::DistanceSample> sd = path::sample_distance(rng, scene_, sp, s->wo);
                    
                    if(!sd && usestrategy == "delta")
                        return contribution ;

                    
                    //2 distance samples 0 and 1, 2 strategies 0 and 1,
                    // we miss one pdf value: equiangular pdf 0 of regular distance sample 1
                    int key = 0;
                    Float regularT = 0.0;
                    
                    if(usestrategy != "delta") {
                        regularT = stats::get<stats::RegularTrackingStrategyDistanceSample,int,Float>(key);

                        //auto mediumindex = scene_->medium_node();
                        //auto & medium = scene_->node_at(mediumindex);
                        //medium.primitive.medium->
                        Ray ray = {
                        sp.geom.p,
                        s->wo};
                                    
                        //sample next scene interaction following light sources situation
                        //scene_->traverse_primitive_nodes(
                        //    [&] (const SceneNode& node, Mat4 global_transform) {
                        const auto & node = scene_->node_at(lightprim.index);
                        Mat4 global_transform = lightprim.global_transform.M;

                                if(node.primitive.light != nullptr && !node.primitive.light->is_infinite()) {
                                    auto & light = *node.primitive.light;
                                    auto possample = rng.next<Light::PositionSampleU>();
                                    auto lightPos = light.sample_position(possample, Transform(global_transform)).value().geom.p;
                                    //LM_INFO("light pos {},{},{}",lightPos.x,lightPos.y,lightPos.z);
                                    auto th = glm::dot((lightPos - ray.o), ray.d); 
                                    auto shortest = ray.o + ray.d * th - lightPos;
                                    auto h = glm::length(shortest);
                                    Float a = -th ;//- th;//ray.o - lightPos;
                                    Float b = 1000.9 ;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                    auto theta_a = glm::atan(a,h);
                                    auto theta_b = glm::atan(b,h);
                                    //auto xi = rng.u();
                                    //auto equiT = th + h * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                    //auto equiT = h * glm::tan((1.0 - xi) * theta_a + xi * theta_b);
                                    
                                    //auto xi = rng.u();
                                    auto equiT = regularT - th;//todo wrong?
                                    auto pdf = h / (h*h+glm::pow(equiT,2.0) *(theta_b-theta_a));
                                    //auto pdf = D / ((theta_b - theta_a) * (D*D+equiT*equiT));
                                    //store pdf for sample of strategy 0
                                    stats::set<stats::DistanceSamplesPDFs,stats::IJ,Float>(stats::IJ::_0_1,pdf);                                
                                }
                        //});
                    }
                    auto _strat_smpl_key = stats::IJ::_0_0;
                    auto equiStratEquiSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                
                    _strat_smpl_key = stats::IJ::_0_1;
                    auto equiStratRegularSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                                            
                    _strat_smpl_key = stats::IJ::_1_0;
                    auto regularStratEquiSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                                             
                    _strat_smpl_key = stats::IJ::_1_1;
                    auto regularStratRegularSmplPDF = stats::get<stats::DistanceSamplesPDFs,stats::IJ,Float>(_strat_smpl_key);                                                          
                    auto w_Equi = glm::pow(1.0 * equiStratEquiSmplPDF,1.0) / (glm::pow(1.0 * equiStratEquiSmplPDF,1.0) + glm::pow(1.0 * regularStratEquiSmplPDF,1.0));
                    auto w_Regular = glm::pow(1.0 * regularStratRegularSmplPDF,1.0) / (glm::pow(1.0 * equiStratRegularSmplPDF,1.0) + glm::pow(1.0 * regularStratRegularSmplPDF,1.0));
                    
                
                    if(usestrategy == "equiangular") {

                        contribution = contribution / equiangularT / equiangularT; //convert measurement

                        //convert directionpdf to vertex area measure
                        pathPdfEquiStrategyEquiSamples *= equiStratEquiSmplPDF * directionpdf
                       / equiangularT / equiangularT;
                        pathPdfRegularStrategyEquiSamples *= regularStratEquiSmplPDF * directionpdf  
                           / regularT / regularT;
                    }
                    if(usestrategy == "regular") {
                        contribution = contribution / regularT / regularT; //convert measurement
                        pathPdfRegularStrategyRegularSamples *= regularStratRegularSmplPDF* directionpdf
                              / regularT / regularT;
                        pathPdfEquiStrategyRegularSamples *= equiStratRegularSmplPDF * directionpdf  
                             / equiangularT / equiangularT;

                    }
                    //first choice: leave sd as is

                    //second choice: sample equiangular
                    //sd = lightDistanceSample;


                    

                    if(usestrategy == "equiangular" ) {
                        sd = lightDistanceSample;
                        sd->weight = sd->weight;
                    } else if(usestrategy == "regular") {
                        sd->weight = sd->weight;
                        //already done
                    } else if(usestrategy == "onesample_mis") {
                        auto chooseStrategy =  rng.u();
                        if(chooseStrategy < 0.5) {
                            sd->weight = w_Regular * sd->weight / regularStratRegularSmplPDF / 0.5;
                        } else {
                            sd = lightDistanceSample;
                            sd->weight = w_Equi * sd->weight / equiStratEquiSmplPDF / 0.5;
                        }
                    } else if(usestrategy == "delta") {
                        //everythin is as it should be
                    }
                    //LM_INFO("weights sum {}" ,Vec3(sd->weight  + lightDistanceSample.weight )[0]);
                    //lightDistanceSample.weight * lightDistanceSample.sp.geom.p 
                    //+ sd.weight * sd->sp.geom.p 
                    //linear combine points
                    //sd->sp.geom.p = lightDistanceSample.sp.geom.p * lightDistanceSample.weight * 0.5 +
                    //                sd->sp.geom.p * sd->weight * 0.5;
                    //draw another uniform value and take one
                    //sd->sp.geom.p = rng.u() > 0.5 ? sd->sp.geom.p : lightDistanceSample.sp.geom.p;

                    //always combine all weights together (right now with 0.5 weight)
                    //sd->weight = 0.5 * lightDistanceSample.weight + 0.5 * sd->weight;

                    

                    //sd->weight = sd->weight / regularStratRegularSmplPDF;
                    //sd->weight = sd->weight / equiStratEquiSmplPDF;
                    
                    
                    // --------------------------------------------------------------------------------

                    // Update throughput
                    throughput *= s->weight * sd->weight;

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
                }

                return contribution;
                
            };


            //constructed 2 paths, one with equiangular, one with regular strategy,
            //combine contributions with MIS

            if(strategy_ == "equiangular") {
                samplePath("equiangular");
                contributionEquiStrategy = Vec3(0.0);
                for(int i = 0; i < EquiMeasurementContributions.size(); i++) {
                    
                    auto boolvec3 =  glm::isinf(EquiMeasurementContributions[i]);
                    if(boolvec3.x || boolvec3.y || boolvec3.z || EquiPdfOfEquiSmpls[i] < std::numeric_limits<Float>::epsilon()) {
                        //LM_INFO("contr: is inf or pdf is zero, pdf {}",  RegularPdfOfRegularSmpls[i]);
                    } else 
                        contributionEquiStrategy += EquiMeasurementContributions[i] / EquiPdfOfEquiSmpls[i];
                }
                contributionEquiStrategy /= static_cast<Float>(EquiMeasurementContributions.size()); //divided by how many paths
                film_->splat(raster_pos, contributionEquiStrategy / static_cast<Float>(spp_));
            }
            else if(strategy_ == "regular") {
                samplePath("regular");
                contributionRegularStrategy = Vec3(0.0);
                for(int i = 0; i < RegularMeasurementContributions.size(); i++) {
                    auto boolvec3 =  glm::isinf(RegularMeasurementContributions[i]);
                    if(boolvec3.x || boolvec3.y || boolvec3.z || RegularPdfOfRegularSmpls[i] < std::numeric_limits<Float>::epsilon()) {
                        //LM_INFO("contr: is inf or pdf is zero, pdf {}",  RegularPdfOfRegularSmpls[i]);
                    } else 
                        contributionRegularStrategy += RegularMeasurementContributions[i] / RegularPdfOfRegularSmpls[i];
                }
                contributionRegularStrategy /= static_cast<Float>(RegularMeasurementContributions.size()); //divided by how many paths
                film_->splat(raster_pos, contributionRegularStrategy / static_cast<Float>(spp_));
            }
            else if(strategy_ == "mis") {
                samplePath("equiangular");
                samplePath("regular");
                auto contribution = Vec3(0.0);
                if(RegularMeasurementContributions.size() != EquiMeasurementContributions.size())
                    LM_INFO("not the same path lengths");

                for(int i = 0; i < RegularMeasurementContributions.size(); i++) {

                    

                    Float weightRegularStrategy = 0.0;
                    auto boolvec3 =  glm::isinf(RegularMeasurementContributions[i]) || glm::isnan(RegularMeasurementContributions[i]);
                    if(boolvec3.x || boolvec3.y || boolvec3.z || RegularPdfOfRegularSmpls[i] < std::numeric_limits<Float>::epsilon()) {
                        //LM_INFO("contr: is inf or pdf is zero, pdf {}",  RegularPdfOfRegularSmpls[i]);
                    } else {
                        weightRegularStrategy = glm::pow(1.0 * RegularPdfOfRegularSmpls[i],mis_power_) /
                            (glm::pow(1.0 * RegularPdfOfRegularSmpls[i],mis_power_) + glm::pow(1.0 * EquiPdfOfRegularSmpls[i],mis_power_));
                    }

                    Float weightEquiStrategy = 0.0;
                    boolvec3 =  glm::isinf(EquiMeasurementContributions[i]) || glm::isnan(EquiMeasurementContributions[i]);
                    if(boolvec3.x || boolvec3.y || boolvec3.z || EquiPdfOfEquiSmpls[i] < std::numeric_limits<Float>::epsilon()) {
                        //LM_INFO("contr: is inf or pdf is zero, pdf {}",  RegularPdfOfRegularSmpls[i]);
                    } else {
                        weightEquiStrategy = glm::pow(1.0 * EquiPdfOfEquiSmpls[i],mis_power_) /
                            (glm::pow(1.0 * EquiPdfOfEquiSmpls[i],mis_power_) + glm::pow(1.0 * RegularPdfOfEquiSmpls[i],mis_power_));

                    }
                    weightRegularStrategy /= weightRegularStrategy + weightEquiStrategy;
                    weightEquiStrategy /= weightRegularStrategy + weightEquiStrategy;
                    contribution += weightRegularStrategy * RegularMeasurementContributions[i] / RegularPdfOfRegularSmpls[i] ;
                    contribution += weightEquiStrategy * EquiMeasurementContributions[i] / EquiPdfOfEquiSmpls[i] ;
                    //contribution += 0.5 * RegularMeasurementContributions[i] / RegularPdfOfRegularSmpls[j] ;
                    //contribution += 0.5 * EquiMeasurementContributions[i] / EquiPdfOfEquiSmpls[j] ;
                    //LM_INFO("weights {} , {}",weightEquiStrategy ,weightRegularStrategy);
                }
                film_->splat(raster_pos,contribution);  
            } else { //onesample_mis or delta strategy 
                auto c = samplePath("delta");
                film_->splat(raster_pos, c/ static_cast<Float>(spp_));
            }
        
            
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

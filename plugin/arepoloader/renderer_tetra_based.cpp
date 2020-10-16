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



//from 200 to 40000 in steps of 200.
static const Vec3 blackBodyCIEs[199] = 
{{1.26427583e-50,2.05026700e-51,4.35619247e-60}
,{1.54858470e-29,5.29443495e-30,1.32918710e-34}
,{6.93284707e-22,2.76749913e-22,8.96796945e-26}
,{6.94430507e-18,3.12903671e-18,4.97602273e-21}
,{1.97606340e-15,9.95025707e-16,5.33968163e-18}
,{9.00083832e-14,4.99691969e-14,6.64752656e-16}
,{1.41331327e-12,8.53160132e-13,2.24566394e-14}
,{1.13178168e-11,7.33639893e-12,3.26356073e-13}
,{5.76373920e-11,3.96882168e-11,2.66962191e-12}
,{2.13386836e-10,1.54668474e-10,1.45164246e-11}
,{6.25770747e-10,4.73816178e-10,5.84610553e-11}
,{1.53970745e-09,1.21010071e-09,1.87605967e-10}
,{3.30842038e-09,2.68456531e-09,5.04956431e-10}
,{6.38894076e-09,5.32853304e-09,1.18285023e-09}
,{1.13249238e-08,9.67171734e-09,2.47813619e-09}
,{1.87214615e-08,1.63198439e-08,4.74005344e-09}
,{2.92168021e-08,2.59266091e-08,8.40966379e-09}
,{4.34549303e-08,3.91644720e-08,1.40114320e-08}
,{6.20616020e-08,5.66976596e-08,2.21387006e-08}
,{8.56251649e-08,7.91594547e-08,3.34352976e-08}
,{1.14682516e-07,1.07134697e-07,4.85754892e-08}
,{1.49709918e-07,1.41147575e-07,6.82441383e-08}
,{1.91118031e-07,1.81654259e-07,9.31184476e-08}
,{2.39250379e-07,2.29039626e-07,1.23852167e-07}
,{2.94384449e-07,2.83617251e-07,1.61062715e-07}
,{3.56734704e-07,3.45631836e-07,2.05321319e-07}
,{4.26456878e-07,4.15263353e-07,2.57146047e-07}
,{5.03653060e-07,4.92632285e-07,3.16997458e-07}
,{5.88377152e-07,5.77805471e-07,3.85276491e-07}
,{6.80640436e-07,6.70802183e-07,4.62324216e-07}
,{7.80417013e-07,7.71600159e-07,5.48423061e-07}
,{8.87648979e-07,8.80141391e-07,6.43799158e-07}
,{1.00225126e-06,9.96337560e-07,7.48625498e-07}
,{1.12411603e-06,1.12007502e-06,8.63025628e-07}
,{1.25311670e-06,1.25121933e-06,9.87077671e-07}
,{1.38911151e-06,1.38961927e-06,1.12081850e-06}
,{1.53194668e-06,1.53511043e-06,1.26424791e-06}
,{1.68145917e-06,1.68751831e-06,1.41733273e-06}
,{1.83747909e-06,1.84666105e-06,1.58001071e-06}
,{1.99983179e-06,2.01235171e-06,1.75219428e-06}
,{2.16833962e-06,2.18440029e-06,1.93377398e-06}
,{2.34282344e-06,2.36261541e-06,2.12462166e-06}
,{2.52310393e-06,2.54680566e-06,2.32459342e-06}
,{2.70900262e-06,2.73678083e-06,2.53353222e-06}
,{2.90034277e-06,2.93235284e-06,2.75127034e-06}
,{3.09695010e-06,3.13333649e-06,2.97763146e-06}
,{3.29865342e-06,3.33955015e-06,3.21243256e-06}
,{3.50528501e-06,3.55081618e-06,3.45548566e-06}
,{3.71668106e-06,3.76696135e-06,3.70659923e-06}
,{3.93268191e-06,3.98781711e-06,3.96557951e-06}
,{4.15313226e-06,4.21321978e-06,4.23223167e-06}
,{4.37788133e-06,4.44301071e-06,4.50636073e-06}
,{4.60678290e-06,4.67703631e-06,4.78777240e-06}
,{4.83969541e-06,4.91514809e-06,5.07627384e-06}
,{5.07648193e-06,5.15720268e-06,5.37167422e-06}
,{5.31701014e-06,5.40306175e-06,5.67378522e-06}
,{5.56115230e-06,5.65259196e-06,5.98242149e-06}
,{5.80878515e-06,5.90566483e-06,6.29740099e-06}
,{6.05978986e-06,6.16215671e-06,6.61854521e-06}
,{6.31405190e-06,6.42194857e-06,6.94567947e-06}
,{6.57146097e-06,6.68492594e-06,7.27863307e-06}
,{6.83191084e-06,6.95097874e-06,7.61723936e-06}
,{7.09529927e-06,7.22000116e-06,7.96133592e-06}
,{7.36152785e-06,7.49189150e-06,8.31076455e-06}
,{7.63050192e-06,7.76655203e-06,8.66537132e-06}
,{7.90213039e-06,8.04388885e-06,9.02500658e-06}
,{8.17632565e-06,8.32381176e-06,9.38952494e-06}
,{8.45300345e-06,8.60623410e-06,9.75878523e-06}
,{8.73208274e-06,8.89107262e-06,1.01326505e-05}
,{9.01348558e-06,9.17824733e-06,1.05109878e-05}
,{9.29713700e-06,9.46768140e-06,1.08936683e-05}
,{9.58296491e-06,9.75930101e-06,1.12805673e-05}
,{9.87089997e-06,1.00530352e-05,1.16715637e-05}
,{1.01608754e-05,1.03488159e-05,1.20665403e-05}
,{1.04528272e-05,1.06465775e-05,1.24653838e-05}
,{1.07466935e-05,1.09462570e-05,1.28679843e-05}
,{1.10424148e-05,1.12477940e-05,1.32742355e-05}
,{1.13399342e-05,1.15511302e-05,1.36840344e-05}
,{1.16391964e-05,1.18562097e-05,1.40972814e-05}
,{1.19401486e-05,1.21629786e-05,1.45138802e-05}
,{1.22427398e-05,1.24713852e-05,1.49337375e-05}
,{1.25469208e-05,1.27813796e-05,1.53567631e-05}
,{1.28526443e-05,1.30929138e-05,1.57828696e-05}
,{1.31598649e-05,1.34059419e-05,1.62119727e-05}
,{1.34685385e-05,1.37204192e-05,1.66439904e-05}
,{1.37786230e-05,1.40363030e-05,1.70788438e-05}
,{1.40900775e-05,1.43535521e-05,1.75164563e-05}
,{1.44028627e-05,1.46721268e-05,1.79567538e-05}
,{1.47169408e-05,1.49919889e-05,1.83996648e-05}
,{1.50322752e-05,1.53131014e-05,1.88451197e-05}
,{1.53488306e-05,1.56354289e-05,1.92930515e-05}
,{1.56665730e-05,1.59589370e-05,1.97433952e-05}
,{1.59854696e-05,1.62835928e-05,2.01960878e-05}
,{1.63054888e-05,1.66093642e-05,2.06510685e-05}
,{1.66265999e-05,1.69362206e-05,2.11082783e-05}
,{1.69487734e-05,1.72641323e-05,2.15676601e-05}
,{1.72719807e-05,1.75930706e-05,2.20291585e-05}
,{1.75961944e-05,1.79230078e-05,2.24927199e-05}
,{1.79213878e-05,1.82539172e-05,2.29582926e-05}
,{1.82475351e-05,1.85857729e-05,2.34258261e-05}
,{1.85746114e-05,1.89185501e-05,2.38952717e-05}
,{1.89025928e-05,1.92522247e-05,2.43665823e-05}
,{1.92314560e-05,1.95867732e-05,2.48397120e-05}
,{1.95611783e-05,1.99221733e-05,2.53146165e-05}
,{1.98917382e-05,2.02584030e-05,2.57912527e-05}
,{2.02231145e-05,2.05954414e-05,2.62695790e-05}
,{2.05552868e-05,2.09332680e-05,2.67495549e-05}
,{2.08882354e-05,2.12718631e-05,2.72311410e-05}
,{2.12219412e-05,2.16112076e-05,2.77142995e-05}
,{2.15563856e-05,2.19512831e-05,2.81989933e-05}
,{2.18915508e-05,2.22920715e-05,2.86851865e-05}
,{2.22274193e-05,2.26335556e-05,2.91728444e-05}
,{2.25639743e-05,2.29757186e-05,2.96619333e-05}
,{2.29011996e-05,2.33185441e-05,3.01524202e-05}
,{2.32390792e-05,2.36620165e-05,3.06442735e-05}
,{2.35775978e-05,2.40061204e-05,3.11374621e-05}
,{2.39167406e-05,2.43508410e-05,3.16319560e-05}
,{2.42564932e-05,2.46961640e-05,3.21277262e-05}
,{2.45968414e-05,2.50420753e-05,3.26247441e-05}
,{2.49377718e-05,2.53885614e-05,3.31229823e-05}
,{2.52792711e-05,2.57356093e-05,3.36224140e-05}
,{2.56213265e-05,2.60832062e-05,3.41230132e-05}
,{2.59639257e-05,2.64313397e-05,3.46247546e-05}
,{2.63070564e-05,2.67799978e-05,3.51276135e-05}
,{2.66507071e-05,2.71291689e-05,3.56315661e-05}
,{2.69948663e-05,2.74788417e-05,3.61365891e-05}
,{2.73395229e-05,2.78290050e-05,3.66426598e-05}
,{2.76846662e-05,2.81796484e-05,3.71497563e-05}
,{2.80302857e-05,2.85307613e-05,3.76578570e-05}
,{2.83763713e-05,2.88823337e-05,3.81669412e-05}
,{2.87229131e-05,2.92343558e-05,3.86769886e-05}
,{2.90699014e-05,2.95868180e-05,3.91879793e-05}
,{2.94173269e-05,2.99397110e-05,3.96998941e-05}
,{2.97651805e-05,3.02930259e-05,4.02127144e-05}
,{3.01134534e-05,3.06467539e-05,4.07264218e-05}
,{3.04621370e-05,3.10008864e-05,4.12409985e-05}
,{3.08112227e-05,3.13554150e-05,4.17564273e-05}
,{3.11607026e-05,3.17103319e-05,4.22726913e-05}
,{3.15105687e-05,3.20656289e-05,4.27897741e-05}
,{3.18608131e-05,3.24212986e-05,4.33076595e-05}
,{3.22114285e-05,3.27773333e-05,4.38263321e-05}
,{3.25624073e-05,3.31337260e-05,4.43457765e-05}
,{3.29137426e-05,3.34904694e-05,4.48659779e-05}
,{3.32654273e-05,3.38475566e-05,4.53869218e-05}
,{3.36174546e-05,3.42049811e-05,4.59085942e-05}
,{3.39698179e-05,3.45627361e-05,4.64309812e-05}
,{3.43225107e-05,3.49208154e-05,4.69540694e-05}
,{3.46755268e-05,3.52792126e-05,4.74778457e-05}
,{3.50288601e-05,3.56379218e-05,4.80022973e-05}
,{3.53825044e-05,3.59969370e-05,4.85274117e-05}
,{3.57364541e-05,3.63562525e-05,4.90531768e-05}
,{3.60907034e-05,3.67158626e-05,4.95795806e-05}
,{3.64452468e-05,3.70757618e-05,5.01066115e-05}
,{3.68000788e-05,3.74359448e-05,5.06342582e-05}
,{3.71551941e-05,3.77964063e-05,5.11625096e-05}
,{3.75105876e-05,3.81571413e-05,5.16913548e-05}
,{3.78662543e-05,3.85181448e-05,5.22207835e-05}
,{3.82221891e-05,3.88794119e-05,5.27507851e-05}
,{3.85783874e-05,3.92409378e-05,5.32813496e-05}
,{3.89348445e-05,3.96027180e-05,5.38124672e-05}
,{3.92915556e-05,3.99647479e-05,5.43441283e-05}
,{3.96485165e-05,4.03270231e-05,5.48763234e-05}
,{4.00057226e-05,4.06895393e-05,5.54090434e-05}
,{4.03631697e-05,4.10522922e-05,5.59422791e-05}
,{4.07208537e-05,4.14152778e-05,5.64760220e-05}
,{4.10787704e-05,4.17784919e-05,5.70102633e-05}
,{4.14369159e-05,4.21419308e-05,5.75449947e-05}
,{4.17952863e-05,4.25055905e-05,5.80802078e-05}
,{4.21538778e-05,4.28694673e-05,5.86158948e-05}
,{4.25126865e-05,4.32335574e-05,5.91520477e-05}
,{4.28717090e-05,4.35978574e-05,5.96886588e-05}
,{4.32309415e-05,4.39623637e-05,6.02257207e-05}
,{4.35903807e-05,4.43270728e-05,6.07632258e-05}
,{4.39500231e-05,4.46919814e-05,6.13011671e-05}
,{4.43098654e-05,4.50570862e-05,6.18395374e-05}
,{4.46699043e-05,4.54223840e-05,6.23783298e-05}
,{4.50301366e-05,4.57878716e-05,6.29175376e-05}
,{4.53905592e-05,4.61535459e-05,6.34571542e-05}
,{4.57511690e-05,4.65194040e-05,6.39971731e-05}
,{4.61119631e-05,4.68854429e-05,6.45375879e-05}
,{4.64729384e-05,4.72516597e-05,6.50783924e-05}
,{4.68340922e-05,4.76180515e-05,6.56195804e-05}
,{4.71954216e-05,4.79846157e-05,6.61611462e-05}
,{4.75569239e-05,4.83513494e-05,6.67030837e-05}
,{4.79185963e-05,4.87182500e-05,6.72453873e-05}
,{4.82804362e-05,4.90853150e-05,6.77880513e-05}
,{4.86424410e-05,4.94525418e-05,6.83310703e-05}
,{4.90046082e-05,4.98199278e-05,6.88744389e-05}
,{4.93669354e-05,5.01874706e-05,6.94181517e-05}
,{4.97294200e-05,5.05551680e-05,6.99622037e-05}
,{5.00920597e-05,5.09230174e-05,7.05065897e-05}
,{5.04548521e-05,5.12910165e-05,7.10513047e-05}
,{5.08177950e-05,5.16591633e-05,7.15963439e-05}
,{5.11808861e-05,5.20274553e-05,7.21417025e-05}
,{5.15441232e-05,5.23958905e-05,7.26873758e-05}
,{5.19075041e-05,5.27644668e-05,7.32333591e-05}
,{5.22710268e-05,5.31331820e-05,7.37796481e-05}
,{5.26346891e-05,5.35020341e-05,7.43262382e-05}
,{5.29984890e-05,5.38710211e-05,7.48731251e-05}};


//#define ALL_LIGHTS


    
class RendererTetraBased final : public Renderer {
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
    Float impact_threshold_;

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

        impact_threshold_ = json::value<Float>(prop, "impact_threshold",1.0);
        LM_INFO("impact_threshold: {}", impact_threshold_);


        #if VOLPT_IMAGE_SAMPLING
        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spi::" + sched_name, make_loc("scheduler"), prop);
        #else
        sched_ = comp::create<scheduler::Scheduler>(
            "scheduler::spp::" + sched_name, make_loc("scheduler"), prop);
        #endif

    }
    virtual Json render() const override {
		scene_->require_renderable();

        stats::clearGlobal<stats::SampleIdCacheHits,int,long long>( );
        stats::clearGlobal<stats::SampleIdCacheMisses,int,long long>( );
        stats::clearGlobal<stats::UsedCachedTetra,int,long long>( );
        stats::clearGlobal<stats::UsedNeighborTetra,int,long long>( );
        stats::clearGlobal<stats::ResampleAccel,int,long long>( );
        stats::clearGlobal<stats::TotalTetraTests,int,long long>( );
        stats::clearGlobal<stats::CachedSampleId,int,long long>( );


        //film_->clear();
        const auto size = film_->size();
        timer::ScopedTimer st;

        int vrls = 0;
        auto & tetraIToLightSegments =  stats::getGlobalRefUnsafe<stats::VRL,stats::TetraIndex,std::vector<LightToCameraRaySegmentCDF>>();
        for(auto p : tetraIToLightSegments)  {
            for(auto v : p.second)  
                vrls++;
            //
        }
        LM_INFO("vrls {}",vrls);

        int lightspertetra = 0;
        //auto & tetraToPointLights = stats::getGlobalRefUnsafe<stats::LightsInTetra,stats::TetraIndex,std::vector<StarSource>>();
        auto & tetraToPointLights = stats::getGlobalRefUnsafe<stats::LightsInTetra,stats::TetraIndex,std::vector<int>>();
        for(auto p : tetraToPointLights)  {
            for(auto v : p.second) 
                lightspertetra++; 
            //LM_INFO("tetra {}, light {}",p.first,v.index);
        }
        LM_INFO("lights per tetra {}",lightspertetra);

        auto & lightSet = stats::getGlobalRefUnsafe<stats::LightSet,int,std::vector<StarSource>>()[0];
        int numlights = lightSet.size();
       
        LM_INFO("total num lights {}",numlights);

        


        LM_INFO("num threads: {}",parallel::num_threads());
        auto sqrtthreads = glm::sqrt(static_cast<Float>(parallel::num_threads()));
        int blockedge = static_cast<int>(sqrtthreads);
        int blocksize = blockedge*blockedge;
        
        /*int localX = pixel_index % blockedge;
        int localY = pixel_index / blockedge;
        int globalX = pixel_index / blocksize;
        int blocksX = size.w / blockedge;
        int blocksY = size.h / blockedge;*/
        
        const auto processed = sched_->run([&](long long pixel_index, long long sample_index, int threadid) {
            LM_KEEP_UNUSED(sample_index);
        
            
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

            auto & emissiveContributions = 
            stats::getRef<stats::EmissiveContribution,int,std::vector<Vec3>>();
            //equiContributions.clear();
            



            int num_pointlights = scene_->num_lights(); 


            thread_local KNNResult point_knn_res;

            // Per-thread random number generator
            thread_local Rng rng(seed_ ? *seed_ + threadid : math::rng_seed());
            

              //to sRGB
            auto m0 = Vec3(3.2404542,-0.9692660,0.0556434);     
            auto m1 = Vec3( -1.5371385, 1.8760108,-0.2040259);  
            auto m2 = Vec3(  -0.4985314, 0.0415560, 1.0572252);
            auto XYZ_TO_sRGB = Mat3(m0,m1,m2);
            
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
            std::vector<int> contributionIndex_emissive = {};
            contributionIndex_equi.resize(max_verts_);
            contributionIndex_regular.resize(max_verts_);
            contributionIndex_emissive.resize(max_verts_);
            
            
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
                    contributionIndex_emissive[i] = 0;
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

       
                    //ask nearest light to current vertex
                    scene_->sample_light_selection_from_pos(0.0,sp.geom.p);

                    Scene::LightPrimitiveIndex lightprim; 

                    auto volumeTMin = 0.0;
                    auto volumeTMax = std::numeric_limits<Float>::max();
                    //no not ray dependent, need global! otherwise i am sacked
                    //auto isectret  = volume_->bound().isect_range({sp.geom.p,s->wo}, volumeTMin, volumeTMax);
                    volumeTMax = glm::max(10000.0, glm::distance(volume_->bound().min,volume_->bound().max));
                    
                    

                    //prepare queries
                    auto pdf_vrl_selection = 1.0;
                    auto pdf_light_selection = 1.0;

                                       
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


                    if(boundaries != nullptr && totalTau > 0.0) { //the chance to have in-scattering                        
                        auto & cameraSegments = *boundaries;


                        //CHOOSE POINTS FOR KNN QUERIES ALONG RAY
                        std::vector<Float> queryTs;
                        std::vector<Float> queryPDFs;
                        std::vector<int> queryTetraInds;
                        {
                            int found_sample = 0;
                            std::vector<Float> zetas;

                            for(int i = 0; i < num_knn_queries_; i++) {
                                //zetas.push_back(rng.u() * totalTau);
                                //zetas.push_back(rng.u());
                                zetas.push_back(-gsl_log1p(-rng.u() * lowDensityNormalizationFactor));

                                queryTs.push_back(0.0);
                                queryTetraInds.push_back(0);
                                queryPDFs.push_back(1);
                            }
                            //TODO PDFS ?!
                            
                            //auto zeta = rng.u() * totalTau; //a new sample within total 
                            //auto zetaRegularPDF = 0.0;
                            auto zetaAccCdf = 0.0;
                            auto zetaTransmittance = 1.0;

                            {
                                auto travelT = 0.0;
                                auto segmentThroughput = 1.0;
                                for(int segmentI = 0; segmentI < segmentCount && found_sample < num_knn_queries_; segmentI++) {
                                    auto & tetrasegment =  cameraSegments[segmentI];
                                    for(int i = 0; i < num_knn_queries_; i++) {
                                        if (zetaAccCdf  + tetrasegment.localcdf  > zetas[i]) {
                                        //if (travelT - minT  > zetas[i] * (totalT-minT)) {
                                            auto normcdf =  zetaAccCdf ;
                                            lm::Float t = sampleCachedICDF_andCDF( zetas[i],zetas[i] , tetrasegment.t ,
                                            normcdf ,  tetrasegment.a ,   tetrasegment.b );
                                            zetaTransmittance *= glm::exp(-  normcdf );
                                            queryTs[i] = travelT + t;   
                                            //queryTs[i] = minT + zetas[i] * (totalT-minT); 
                                            //queryPDFs[i] /=  totalTau;
                                            queryPDFs[i] = 1.0;
                                            /*debatable: queryPDFs[i] = tetrasegment.b + tetrasegment.t * tetrasegment.a; 
                                            queryPDFs[i] *= zetaTransmittance; 
                                            queryPDFs[i] /=  lowDensityNormalizationFactor;*/
                                            queryTetraInds[i] = tetrasegment.tetraI;
                                            zetaAccCdf += normcdf;
                                            //break;
                                            found_sample++;
                                        }
                                    }
                                    travelT += tetrasegment.t;
                                    zetaAccCdf += tetrasegment.localcdf;
                                    zetaTransmittance *= glm::exp(-zetaAccCdf );

                                }
                            }
                        }

                        //auto a = sp.geom.p + a_d * currentTransmittanceDistanceSample; //dont use camera ray but camera point
                        auto maxAllowedT = regularT;//tetrasegment.t;//glm::min(
                            
                        if(sample_vrls_) {
                            for(int queryI = 0; queryI < num_knn_queries_; queryI++) {
                                int tetraI = queryTetraInds[queryI];

                                //int acceptedBFSLayer = 10;
                                
                                stats::set<stats::TetraIdGuess,int,int>(0,tetraI);

                                bool continueBFS = true;
                                if (tetraIToLightSegments.find(tetraI) == tetraIToLightSegments.end())
                                    continue;

                                auto & vrlsInThisTetra = tetraIToLightSegments[tetraI];


                                for(auto & vrl : vrlsInThisTetra) {
                                        
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
                                    //auto vrl_area_measure_conv = 1.0/ (vrl.tSoFar + vrl_vlp_t) / (vrl.tSoFar + vrl_vlp_t);
                                    //TODO whcih one ?!
                                    auto vrl_area_measure_conv = 1.0/ ( vrl_vlp_t) / ( vrl_vlp_t);
                                    vrl_scatter_pdf *= vrl_area_measure_conv;

                           

                                    //have sampled virtual ray light, now sample point on camera ray
                                    auto th = glm::dot((vrl_vlp - a), a_d); 
                                    //th = glm::max(0.0, glm::min(maxAllowedT,th ));
                                    auto shortest = a + a_d * th - vrl_vlp;
                                    auto equih = glm::length(shortest);
                                    //what if it the shortest d is not within the line segment? TODO
                                    Float a_ = minT-th;//- th;//ray.o - lightPos;
                                    Float b_ = totalT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                    //Float b_ = maxT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                    



                                    auto cam_scatter_pdf_regular = regularStratRegularSmplPDF;
                                    auto segmentThroughput_regular = Vec3(
                                                A_R_A_V_S * regularMuT * glm::exp(-A_R_A_V_T*tauUntilRegularT),
                                                A_G_A_V_S * regularMuT * glm::exp(-A_G_A_V_T*tauUntilRegularT),
                                                A_B_A_V_S * regularMuT * glm::exp(-A_B_A_V_T*tauUntilRegularT)
                                                );

                                    auto t_regular = regularT;
                                    int tetraILandedIn_regular = lastDistanceSampleTetraI;
                                    

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
                                    auto connectiondir_equi =  glm::normalize(camera_sample_point_equi-vrl_vlp);
                                    auto connection_area_measure_equi = 1.0 / glm::abs(glm::dot(camera_sample_point_equi-vrl_vlp,camera_sample_point_equi-vrl_vlp));
                                    


                                    auto sp1_equi = SceneInteraction::make_medium_interaction(
                                                scene_->medium_node(),
                                                PointGeometry::make_degenerated(camera_sample_point_equi)
                                            );

                                    auto sp2 = SceneInteraction::make_medium_interaction(
                                                scene_->medium_node(),
                                                PointGeometry::make_degenerated(vrl_vlp)
                                            );

                                    //REGULAR PDF
                                    //auto cam_area_measure_conv = 1.0 / (travelT + t) / (travelT + t) ;
                                    auto cam_area_measure_conv_regular = 1.0 / t_regular / t_regular;
                                    cam_scatter_pdf_regular *= cam_area_measure_conv_regular;
                                    auto camera_sample_point_regular = a + a_d * t_regular;
                                    auto connectiondir_regular =  glm::normalize(camera_sample_point_regular-vrl_vlp);
                                    auto connection_area_measure_regular = 1.0 / glm::abs(glm::dot(camera_sample_point_regular-vrl_vlp,camera_sample_point_regular-vrl_vlp));
                                    auto sp1_regular = SceneInteraction::make_medium_interaction(
                                                scene_->medium_node(),
                                                PointGeometry::make_degenerated(camera_sample_point_regular)
                                            );

                                    //EQUI PDF OF REGULAR
                                    auto cam_scatter_pdf_equi_of_regular = equih / (equih*equih+glm::pow(t_regular-th,2.0) *(theta_b-theta_a));
                                    cam_scatter_pdf_equi_of_regular *= cam_area_measure_conv_regular;


                                            
                                    //camera throughput until sample points: 
                                    //transmittance!
                                    auto travelT = 0.0;
                                    auto accCdf = 0.0;
                                    auto segmentThroughput_equi = Vec3(1.0);
                                    int tetraILandedIn_equi = -1;
                                    auto cam_scatter_pdf_regular_of_equi = 1.0;
                                    bool foundEqui = false;

                                    Vec3 accumulatedEmission = Vec3(0);

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
                                            //break;
                                        }

                                        auto fromTemp = tetrasegment.b_kelv;
                                        auto toTemp = tetrasegment.b_kelv + tetrasegment.t * tetrasegment.a_kelv;
                                        Vec3 fromrgb = Vec3(0);
                                        Vec3 torgb = Vec3(0);
                                            //LM_INFO("from temp {}",fromTemp);
                                            //LM_INFO("to temp {}",toTemp);
                                        if(fromTemp > 200.0) {
                                            auto lookupI = static_cast<int>((200.0 * (fromTemp - 200.0) / (40000.0 - 200.0)));
                                            fromrgb = tetrasegment.b * XYZ_TO_sRGB * blackBodyCIEs[glm::min(199,glm::max(0,lookupI))];
                                        // LM_INFO("from temp {}",fromTemp);
                                        }
                                        if(toTemp > 200.0) {
                                            auto lookupI = static_cast<int>((200.0 * (toTemp - 200.0) / (40000.0 - 200.0)));
                                            torgb =  (tetrasegment.b + tetrasegment.a * tetrasegment.t) * XYZ_TO_sRGB * blackBodyCIEs[glm::min(199,glm::max(0,lookupI))];
                                            //LM_INFO("temp {},{},{}",torgb.x,torgb.y,torgb.z);
                                        }
                                        accumulatedEmission += 10000.0 * (0.5*torgb * tetrasegment.t *  tetrasegment.t + fromrgb * tetrasegment.t);
                                        
                                    

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
                                    const auto wo_equi = glm::normalize(vrl_vlp - camera_sample_point_equi);

                                    //TODO überprüfen: alle PDFs beim light sampling
                                    
                                    //TODO evaluate what is wo and what is wi, is it EL or LE ?! 
                                    auto fs1_equi = path::eval_contrb_direction(scene_, sp1_equi, a_d, wo_equi, comp, TransDir::EL, true);

                                    auto fs2_equi = path::eval_contrb_direction(scene_, sp2, b_d, -wo_equi, comp, TransDir::LE, true);

                                    //fs1 *= tetrasegment.b + tetrasegment.a * t; //phase times \mu_s
                                    //point light has no phase auto fs2 = path::eval_contrb_direction(scene_, sp2, b_d, -wo, comp, TransDir::LE, true);
                                    //fs2 *= vrl.b + vrl.a * vrl_vlp_t;//phase times \mu_s
                                    const auto& primitive = scene_->node_at(scene_->medium_node()).primitive;
                                    //todo directions correct ?
                                    auto solidangledirectionpdf1_equi = primitive.medium->phase()->pdf_direction(sp1_equi.geom, a_d, wo_equi);
                                    auto solidangledirectionpdf2_equi = primitive.medium->phase()->pdf_direction(sp2.geom, b_d, wo_equi);
                                    //auto solidangledirectionpdf2 = primitive.medium->phase()->pdf_direction(sp2.geom, b_d, wo);
                                    fs1_equi *= solidangledirectionpdf1_equi;
                                    fs2_equi *= solidangledirectionpdf2_equi;
                                    //fs2 *= solidangledirectionpdf2;
                                    //solidangledirectionpdf1_equi *= cam_area_measure_conv_equi; //convert pdf to vertex area
                                    //solidangledirectionpdf2 *= vrl_area_measure_conv; //convert pdf to vertex area
                                    
                                    Vec3 VRL_C_equi =  vrl.weight;
                                    //transmittance up to scatter event on vrl AND scatter coefficient 
                                    {
                                        auto sigma_t = (vrl.b + vrl.a * vrl_vlp_t);
                                        auto tau = vrl.cdfSoFar + 0.5 * vrl.a *  vrl_vlp_t * vrl_vlp_t + vrl.b * vrl_vlp_t;
                                        VRL_C_equi.r *= A_R_A_V_S *sigma_t* glm::exp(-A_R_A_V_T*tau);
                                        VRL_C_equi.g *= A_G_A_V_S *sigma_t* glm::exp(-A_G_A_V_T*tau);
                                        VRL_C_equi.b *= A_B_A_V_S *sigma_t* glm::exp(-A_B_A_V_T*tau);
                                    }

                                    
                                    
                                    //connection transmittance
                                    VRL_C_equi.r *= A_R_A_V_S * glm::exp(-accCDF * A_R_A_V_T);
                                    VRL_C_equi.g *= A_G_A_V_S * glm::exp(-accCDF * A_G_A_V_T);
                                    VRL_C_equi.b *= A_B_A_V_S * glm::exp(-accCDF * A_B_A_V_T);


                                    //auto segment_contribution = throughputCam * Tr * connection_area_measure * fs1 * fs2 * VRL_C;
                                    auto segment_contribution_equi = throughputCam_equi  * connection_area_measure_equi * fs1_equi * fs2_equi * vrl_area_measure_conv * VRL_C_equi;
                                    auto segment_pdf_equi = 
                                    queryPDFs[queryI] *
                                    pathPDF 
                                    * vrl_scatter_pdf
                                    * cam_scatter_pdf_equi
                                    * pdf_vrl_selection
                                    // / vrlCount
                                    //* solidangledirectionpdf1
                                    //* solidangledirectionpdf2
                                    // connection probability stated to be one * connection_area_measure
                                    ;
                                    auto segment_pdf_equi_of_regular = 
                                    queryPDFs[queryI] *
                                    pathPDF * vrl_scatter_pdf * cam_scatter_pdf_equi_of_regular * pdf_vrl_selection;
                                    
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
                                    const auto wo_regular = glm::normalize(vrl_vlp - camera_sample_point_regular);

                                    //TODO überprüfen: alle PDFs beim light sampling
                                    
                                    //TODO evaluate what is wo and what is wi, is it EL or LE ?! 
                                    auto fs1_reqular = path::eval_contrb_direction(scene_, sp1_regular, a_d, wo_regular, comp, TransDir::EL, true);
                                    auto fs2_reqular = path::eval_contrb_direction(scene_, sp2, b_d, -wo_regular, comp, TransDir::EL, true);
                                    //fs1 *= tetrasegment.b + tetrasegment.a * t; //phase times \mu_s
                                    //point light has no phase auto fs2 = path::eval_contrb_direction(scene_, sp2, b_d, -wo, comp, TransDir::LE, true);
                                    //fs2 *= vrl.b + vrl.a * vrl_vlp_t;//phase times \mu_s
                                    //const auto& primitive = scene_->node_at(scene_->medium_node()).primitive;
                                    //todo directions correct ?
                                    auto solidangledirectionpdf1_regular = primitive.medium->phase()->pdf_direction(sp1_regular.geom, a_d, wo_regular);
                                    auto solidangledirectionpdf2_regular = primitive.medium->phase()->pdf_direction(sp2.geom, b_d, wo_regular);
                                    //auto solidangledirectionpdf2 = primitive.medium->phase()->pdf_direction(sp2.geom, b_d, wo);
                                    fs1_reqular *= solidangledirectionpdf1_regular;
                                    fs2_reqular *= solidangledirectionpdf2_regular;
                                    //fs2 *= solidangledirectionpdf2;
                                    //solidangledirectionpdf1_regular *= cam_area_measure_conv_regular; //convert pdf to vertex area
                                    //solidangledirectionpdf2 *= vrl_area_measure_conv; //convert pdf to vertex area
                                    
                                    Vec3 VRL_C_regular =  vrl.weight;
                                    //transmittance up to scatter event on vrl AND scatter coefficient 
                                    {
                                        auto sigma_t = (vrl.b + vrl.a * vrl_vlp_t);
                                        auto tau = vrl.cdfSoFar + 0.5 * vrl.a *  vrl_vlp_t * vrl_vlp_t + vrl.b * vrl_vlp_t;
                                        VRL_C_regular.r *= A_R_A_V_S *sigma_t* glm::exp(-A_R_A_V_T*tau);
                                        VRL_C_regular.g *= A_G_A_V_S *sigma_t* glm::exp(-A_G_A_V_T*tau);
                                        VRL_C_regular.b *= A_B_A_V_S *sigma_t* glm::exp(-A_B_A_V_T*tau);
                                    }

                                    
                                    
                                    //connection transmittance
                                    

                                    VRL_C_regular.r *= A_R_A_V_S * glm::exp(-accCDF * A_R_A_V_T);
                                    VRL_C_regular.g *= A_G_A_V_S * glm::exp(-accCDF * A_G_A_V_T);
                                    VRL_C_regular.b *= A_B_A_V_S * glm::exp(-accCDF * A_B_A_V_T);

                                    //auto segment_contribution = throughputCam * Tr * connection_area_measure * fs1 * fs2 * VRL_C;
                                    auto segment_contribution_regular = throughputCam_regular  * connection_area_measure_regular * fs1_reqular * fs2_reqular * vrl_area_measure_conv * VRL_C_regular;
                                    auto segment_pdf_regular = 
                                    queryPDFs[queryI] *
                                    pathPDF //pathPdfEquiStrategyEquiSamples  //yes it is the main path part
                                    * vrl_scatter_pdf
                                    * cam_scatter_pdf_regular
                                    * pdf_light_selection
                                    // / vrlCount
                                    //* solidangledirectionpdf1
                                    //* solidangledirectionpdf2
                                    // connection probability stated to be one * connection_area_measure
                                    ;
                                    auto segment_pdf_regular_of_equi = 
                                        queryPDFs[queryI] *
                                        pathPDF * vrl_scatter_pdf * cam_scatter_pdf_regular_of_equi * pdf_light_selection;
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
                                    

                                    if(contributionIndex_emissive[num_verts] >= emissiveContributions[num_verts].size()) {
                                        emissiveContributions[num_verts].
                                        push_back(accumulatedEmission);
                                        
                                    } else {
                                        emissiveContributions[num_verts][contributionIndex_emissive[num_verts]] = accumulatedEmission;
                                    }
                                    
                                    contributionIndex_emissive[num_verts]++;
                                    
                                }
                            }
                        }

                        if(sample_lights_) {
                        
                            stats::clear<stats::DuplicateWatchdog,int,int>();

                            //for(int segmentI = 0; segmentI < segmentCount; segmentI++) {
                                //auto & tetrasegment =  cameraSegments[segmentI];
                                //int tetraI = tetrasegment.tetraI;
#ifdef ALL_LIGHTS
                            for(auto & starSource : lightSet) {
                                int queryI = 0;
#else 

#ifdef USE_KNN_EMBREE
                            
                            for(int queryI = 0; queryI < num_knn_queries_; queryI++) {
                                auto queryPos = a + a_d * queryTs[queryI];
                                //LM_INFO("going to query knn ");
                                point_knn_res.k = 0;

                                pointLightAccel_->queryKnn(queryPos.x,queryPos.y,queryPos.z,
                                            0.001, point_knn_res );
                                //LM_INFO("queried {} objects ",point_knn_res.numResults);

                                for(int knnI = 0; knnI < point_knn_res.numResults; knnI ++) {
                                    auto starI = point_knn_res.results[knnI];
                                    auto & starSource = lightSet[starI];
#else
                            for(int queryI = 0; queryI < num_knn_queries_; queryI++) {
                                int tetraI = queryTetraInds[queryI];

                                //int acceptedBFSLayer = 10;
                                
                                stats::set<stats::TetraIdGuess,int,int>(0,tetraI);

                                bool continueBFS = true;
                                if (tetraToPointLights.find(tetraI) == tetraToPointLights.end())
                                    continue;

                                auto & lightNodeIndicesInThisTetra = tetraToPointLights[tetraI];

                                for(auto & starSourceI : lightNodeIndicesInThisTetra) {
#endif                              
                                    auto & starSource = lightSet[starSourceI];
                                    if(!stats::has<stats::DuplicateWatchdog,int,int>(starSource.index)){ //if this light hasnt been handled yet, perform lighting!
                                        stats::set<stats::DuplicateWatchdog,int,int>(starSource.index,1); 
#endif
                                        Vec3 lightPos = starSource.position;

                                        auto b = lightPos;

                                        auto th = glm::dot((b - a), a_d); 
                                        auto shortest = a + a_d * th - b;
                                        
                                        auto equih = glm::length(shortest);
                                        Float starintens = glm::max(starSource.intensity[0],glm::max(starSource.intensity[1],starSource.intensity[2]));

                                        if(starintens/equih < impact_threshold_) //BIAS
                                            continue;

                                        Float a_ =  minT - th;//- th;//ray.o - lightPos;

                                        Float b_ =  totalT - th;//-th;//ray.o + ray.d * 999999.0 - lightPos;
                                        

                                        auto cam_scatter_pdf_regular = regularStratRegularSmplPDF;
                                        auto segmentThroughput_regular = Vec3(
                                                    A_R_A_V_S * regularMuT * glm::exp(-A_R_A_V_T*tauUntilRegularT),
                                                    A_G_A_V_S * regularMuT * glm::exp(-A_G_A_V_T*tauUntilRegularT),
                                                    A_B_A_V_S * regularMuT * glm::exp(-A_B_A_V_T*tauUntilRegularT)
                                                    );

                                        auto t_regular = regularT;
                                        int tetraILandedIn_regular = lastDistanceSampleTetraI;
                                        

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


                                                
                                        //camera throughput until sample points: 
                                        //transmittance!
                                        auto travelT = 0.0;
                                        auto accCdf = 0.0;
                                        auto segmentThroughput_equi = Vec3(1.0);
                                        int tetraILandedIn_equi = -1;
                                        auto cam_scatter_pdf_regular_of_equi = 1.0;
                                        bool foundEqui = false;

                                        Vec3 accumulatedEmission = Vec3(0);

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
                                                //break;
                                            }

                                            auto fromTemp = tetrasegment.b_kelv;
                                            auto toTemp = tetrasegment.b_kelv + tetrasegment.t * tetrasegment.a_kelv;
                                            Vec3 fromrgb = Vec3(0);
                                            Vec3 torgb = Vec3(0);
                                                //LM_INFO("from temp {}",fromTemp);
                                                //LM_INFO("to temp {}",toTemp);
                                            if(fromTemp > 200.0) {
                                                auto lookupI = static_cast<int>((200.0 * (fromTemp - 200.0) / (40000.0 - 200.0)));
                                                fromrgb = tetrasegment.b * XYZ_TO_sRGB * blackBodyCIEs[glm::min(199,glm::max(0,lookupI))];
                                            // LM_INFO("from temp {}",fromTemp);
                                            }
                                            if(toTemp > 200.0) {
                                                auto lookupI = static_cast<int>((200.0 * (toTemp - 200.0) / (40000.0 - 200.0)));
                                                torgb =  (tetrasegment.b + tetrasegment.a * tetrasegment.t) * XYZ_TO_sRGB * blackBodyCIEs[glm::min(199,glm::max(0,lookupI))];
                                                //LM_INFO("temp {},{},{}",torgb.x,torgb.y,torgb.z);
                                            }
                                            accumulatedEmission += 10000.0 * (0.5*torgb * tetrasegment.t *  tetrasegment.t + fromrgb * tetrasegment.t);
                                            
                                        

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
                                        
                                        Vec3 C_equi =  starSource.intensity;
                                        
                                        
                                        C_equi.r *= A_R_A_V_S * glm::exp(-accCDF * A_R_A_V_T);
                                        C_equi.g *= A_G_A_V_S * glm::exp(-accCDF * A_G_A_V_T);
                                        C_equi.b *= A_B_A_V_S * glm::exp(-accCDF * A_B_A_V_T);


                                        //auto segment_contribution = throughputCam * Tr * connection_area_measure * fs1 * fs2 * VRL_C;
                                        auto segment_contribution_equi = throughputCam_equi  * connection_area_measure_equi * fs1_equi * C_equi;
                                        auto segment_pdf_equi = 
                                        queryPDFs[queryI] *
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
                                        queryPDFs[queryI] *
                                        pathPDF * cam_scatter_pdf_equi_of_regular * pdf_light_selection;
                                        
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
                                        
                                        Vec3 C_regular =  starSource.intensity;
                                        

                                        C_regular.r *= A_R_A_V_S * glm::exp(-accCDF * A_R_A_V_T);
                                        C_regular.g *= A_G_A_V_S * glm::exp(-accCDF * A_G_A_V_T);
                                        C_regular.b *= A_B_A_V_S * glm::exp(-accCDF * A_B_A_V_T);

                                        //auto segment_contribution = throughputCam * Tr * connection_area_measure * fs1 * fs2 * VRL_C;
                                        auto segment_contribution_regular = throughputCam_regular  * connection_area_measure_regular * fs1_reqular * C_regular;
                                        auto segment_pdf_regular = 
                                        queryPDFs[queryI] *
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
                                            queryPDFs[queryI] *
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
                                        

                                        if(contributionIndex_emissive[num_verts] >= emissiveContributions[num_verts].size()) {
                                            emissiveContributions[num_verts].
                                            push_back(accumulatedEmission);
                                            
                                        } else {
                                            emissiveContributions[num_verts][contributionIndex_emissive[num_verts]] = accumulatedEmission;
                                        }
                                        
                                        contributionIndex_emissive[num_verts]++;
#ifdef ALL_LIGHTS
#else
                                    }//if not yet seen

                                }//for every light that has impact
#endif
                            }

                        }
                            

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

                    //if(validSamplesEqui != validSamplesRegular)
                    //    LM_ERROR("do not have same amound of VALID samples: {} equi, {} regular", validSamplesEqui,validSamplesRegular );

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


            {//emissive
                auto emissiveContr = Vec3(0);
                for(int i = 0; i < max_verts_; i++) { //for paths of all lengths
                    auto & cs = stats::getRef<stats::EmissiveContribution,int,std::vector<Vec3>>(i);
                    auto pathsWithLengtI = contributionIndex_emissive[i];
                    Vec3 ci = Vec3(0);
                    Float js = 0.0;
                    //for(int j = 0; j < cs.size(); j++) { //for pahts with length i
                    for(int j = 0; j < pathsWithLengtI; j++) { //for pahts with length i
                        //LM_INFO("path length {}, sample {} , pdf ,  contr {},{},{}", pathsWithLengtI,j,cs[j].x,cs[j].y,cs[j].z );
                        ci += cs[j] ;
                        js += 1.0;
                    }
                    if(js != 0.0) {
                        ci /= static_cast<Float>(js); //divided by how many paths
                    }
                    emissiveContr += ci;
                }
                //contributionEquiStrategy /= static_cast<Float>(EquiMeasurementContributions.size()); //divided by how many paths
                film_->splat(raster_pos, emissiveContr / static_cast<Float>(spp_));
            }

            
        
            
        },  
        [&](auto pxlindx,auto smplindx,auto threadid) {
            stats::clear<stats::SampleIdCacheHits,int,long long>( );
            stats::clear<stats::SampleIdCacheMisses,int,long long>( );
            stats::clear<stats::UsedCachedTetra,int,long long>( );
            stats::clear<stats::UsedNeighborTetra,int,long long>( );
            stats::clear<stats::ResampleAccel,int,long long>( );
            stats::clear<stats::TotalTetraTests,int,long long>( );
            
            stats::clear<stats::DuplicateWatchdog,int,int>();
            stats::clear<stats::CachedSampleId,int,long long>();


        } , 
        [&](auto pxlindx,auto smplindx,auto threadid) {
            stats::mergeToGlobal<stats::SampleIdCacheHits,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::SampleIdCacheMisses,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::UsedCachedTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::UsedNeighborTetra,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::ResampleAccel,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::mergeToGlobal<stats::TotalTetraTests,int,long long>( 
                [](long long & v0,long long & v1 ) { return v0 + v1;} );
            stats::clear<stats::DuplicateWatchdog,int,int>();

        }
        );

        auto smplhits = stats::getGlobal<stats::SampleIdCacheHits,int,long long>(0 );
        auto smplmisses = stats::getGlobal<stats::SampleIdCacheMisses,int,long long>(0 );
        auto tetrahits =  stats::getGlobal<stats::UsedCachedTetra,int,long long>( 0);
        auto tetraneighborhits = stats::getGlobal<stats::UsedNeighborTetra,int,long long>(0 );
        auto accelsmpls = stats::getGlobal<stats::ResampleAccel,int,long long>( 0);
        auto totaltetratests = stats::getGlobal<stats::TotalTetraTests,int,long long>(0 );

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

LM_COMP_REG_IMPL(RendererTetraBased, "renderer::tetra_based");

LM_NAMESPACE_END(LM_NAMESPACE)

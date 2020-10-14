#include <lm/core.h>
#include <lm/lm.h>
#include <lm/volume.h>
#include <lm/stats.h>


//#include "../../../ArepoVTK/arepo/include/mesh/voronoi/voronoi.h"
#include "voronoi_3db.h" //shit didnt want to include any arepo code here

//#define USE_KNN_EMBREE

//scaling factors for optical thickness
#define A_B_A_V_T 1.3034926470588235//1.324
#define A_G_A_V_T 1.0346507352941177//1.324
#define A_R_A_V_T 0.8905061025223759//0.748


//scalling factors for scattering coefficient
#define A_B_A_V_S 1.3034926470588235
#define A_G_A_V_S 1.0346507352941177
#define A_R_A_V_S 0.8905061025223759


#define INSIDE_TOLERANCE 1.0 * std::numeric_limits<lm::Float>::epsilon()

inline lm::Float sampleCachedICDF_andCDF(lm::Float logxi,lm::Float xi, lm::Float tmax, lm::Float & out_cdf,lm::Float a, lm::Float b){
    

    //use tau*_t1 (t) which is the integral from t1 to t minus the integral from 0 to t1
    //auto y = logxi + a *  0.5 * tmin*tmin   + b * tmin  - out_cdf;
    auto y = logxi - out_cdf;
    lm::Float freeT;

    if (abs(a) < INSIDE_TOLERANCE
    && abs(b) < INSIDE_TOLERANCE) {
        freeT = std::numeric_limits<lm::Float>::max();
    } else if(abs(a) < INSIDE_TOLERANCE) {
        //freeT = 2.0 * y / (glm::sqrt( b*b +  2.0 * y * a  ) + b);
        freeT = y / b;
    }
    else if(abs(b) < INSIDE_TOLERANCE) {
        //freeT = (glm::sqrt( b*b +  2.0 * y * a  ) - b) / a;
        freeT = glm::sqrt(2.0 * y * a) / a;
    } else if (abs(a) < abs(b)) {
        freeT = 2.0 * y / (glm::sqrt( b*b +  2.0 * y * a  ) + b);
    } else if (abs(b) < abs(a)) {
        freeT = (glm::sqrt( b*b +  2.0 * y * a  ) - b) / a;
    }
  
    freeT = isnan(freeT) ? std::numeric_limits<lm::Float>::max() : freeT;
    //for evaluating tau, limit free path to tmax
    lm::Float t = glm::min(freeT , tmax);
    lm::Float acc_cdf = t  * b  +  a  * 0.5 * t * t;
    
    out_cdf += glm::max(0.0,
        acc_cdf) ; //cdf within tmin and  min of (t , tmax) 
    return freeT;//returns sth between tmin - tmin (so 0) and tmax - tmin 
}
inline lm::Float sampleCDF(  lm::Float toT,lm::Float a, lm::Float b) {
    return ( b * toT  +  a * 0.5 *toT * toT   );
}

namespace ArepoLoaderInternals {
    struct IArepoMeshMock {
        virtual tetra * getDT() = 0;
        virtual point * getDP() = 0;
        virtual std::vector<lm::Float> & getdensities() = 0;
        virtual int getNdt() = 0;
        virtual int getNdp() = 0;

        virtual BBox WorldBound() = 0;
        virtual lm::Float getDensity(int index) = 0;
        virtual lm::Float getTemperature(int index) = 0;

        virtual lm::Float max_density() = 0;
    };





    static const lm::Float g_band_wavelengths_nm[89] = 
    {
        3630.0,3655.0,3680.0,3705.0,3730.0,3755.0,3780.0,3805.0,3830.0,3855.0,3880.0,3905.0
        ,3930.0,3955.0,3980.0,4005.0,4030.0,4055.0,4080.0,4105.0,4130.0,4155.0,4180.0,4205.0
        ,4230.0,4255.0,4280.0,4305.0,4330.0,4355.0,4380.0,4405.0,4430.0,4455.0,4480.0,4505.0
        ,4530.0,4555.0,4580.0,4605.0,4630.0,4655.0,4680.0,4705.0,4730.0,4755.0,4780.0,4805.0
        ,4830.0,4855.0,4880.0,4905.0,4930.0,4955.0,4980.0,5005.0,5030.0,5055.0,5080.0,5105.0
        ,5130.0,5155.0,5180.0,5205.0,5230.0,5255.0,5280.0,5305.0,5330.0,5355.0,5380.0,5405.0
        ,5430.0,5455.0,5480.0,5505.0,5530.0,5555.0,5580.0,5605.0,5630.0,5655.0,5680.0,5705.0
        ,5730.0,5755.0,5780.0,5805.0,5830.0
    };

    static const lm::Float g_band_response[89] = 
    {
        0.000e+00,3.000e-04,8.000e-04,1.300e-03,1.900e-03,2.400e-03,3.400e-03
        ,5.500e-03,1.030e-02,1.940e-02,3.260e-02,4.920e-02,6.860e-02,9.000e-02
        ,1.123e-01,1.342e-01,1.545e-01,1.722e-01,1.873e-01,2.003e-01,2.116e-01
        ,2.214e-01,2.301e-01,2.378e-01,2.448e-01,2.513e-01,2.574e-01,2.633e-01
        ,2.691e-01,2.747e-01,2.801e-01,2.852e-01,2.899e-01,2.940e-01,2.979e-01
        ,3.016e-01,3.055e-01,3.097e-01,3.141e-01,3.184e-01,3.224e-01,3.257e-01
        ,3.284e-01,3.307e-01,3.327e-01,3.346e-01,3.364e-01,3.383e-01,3.403e-01
        ,3.425e-01,3.448e-01,3.472e-01,3.495e-01,3.519e-01,3.541e-01,3.562e-01
        ,3.581e-01,3.597e-01,3.609e-01,3.613e-01,3.609e-01,3.595e-01,3.581e-01
        ,3.558e-01,3.452e-01,3.194e-01,2.807e-01,2.339e-01,1.839e-01,1.352e-01
        ,9.110e-02,5.480e-02,2.950e-02,1.660e-02,1.120e-02,7.700e-03,5.000e-03
        ,3.200e-03,2.100e-03,1.500e-03,1.200e-03,1.000e-03,9.000e-04,8.000e-04
        ,6.000e-04,5.000e-04,3.000e-04,1.000e-04,0.000e+00
    };

    static const lm::Float g_band_weights_sum = 16.740800000000004;

    static const lm::Vec3 cie_color_match[81] = {
        {0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
        {0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
        {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
        {0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
        {0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
        {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
        {0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
        {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
        {0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
        {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
        {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
        {0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
        {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
        {0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
        {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
        {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
        {0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
        {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
        {0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
        {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
        {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
        {0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
        {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
        {0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
        {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
        {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
        {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
    };

    static const lm::Float sun_solid_angle = 6.793970367229902e-05;


    static const lm::Float g_to_temp_params[10] = {
         9.17721455e+02,  2.00722937e+05, -7.69752165e+06,  1.57132181e+08,
        -1.76560740e+09,  1.17234417e+10, -4.71587498e+10,  1.12931504e+11,
        -1.48156059e+11,  8.20461651e+10
    };

    static void mag_to_flux(lm::Float & u, lm::Float & g, lm::Float & r, lm::Float & i, lm::Float & z) {

        u = glm::pow(10.0,-u/2.5) * (3631.0*1e-26);
        g = glm::pow(10.0,-g/2.5) * (3631.0*1e-26);
        r = glm::pow(10.0,-r/2.5) * (3631.0*1e-26);
        i = glm::pow(10.0,-i/2.5) * (3631.0*1e-26);
        z = glm::pow(10.0,-z/2.5) * (3631.0*1e-26);

    }


    static void flux_to_mag(lm::Float & u, lm::Float & g, lm::Float & r, lm::Float & i, lm::Float & z) {

        u = -2.5*log10(u/(3631.0*1e-26));
        g = -2.5*log10(g/(3631.0*1e-26));
        r = -2.5*log10(r/(3631.0*1e-26));
        i = -2.5*log10(i/(3631.0*1e-26));
        z = -2.5*log10(z/(3631.0*1e-26));

    }

    static void flux_to_mag_g(lm::Float & g) {

        g = -2.5*log10(g/(3631.0*1e-26));
       
    }



    static lm::Float poly10th(lm::Float x) {
        auto ret = 0.0;
        lm::Float x_ = 1.0;
        for(int i = 0; i < 10; i++) {
            ret += x_ * g_to_temp_params[i];
            x_ *= x;
        }
        return ret;
            
    }


    static lm::Float const h = 6.62607004e-34;
    static lm::Float const kb = 1.380649e-23;
    static lm::Float const c = 299792458.0;
    static lm::Float const hc = h * c;
    static lm::Float blackbody(lm::Float wavelengthNM,lm::Float temperatureKelvin) {
        auto fr = c / (wavelengthNM *  1e-9);
        auto t1 = 2.0 * h *fr *fr*fr;
        auto t2 = 1.0 / (c*c* (glm::exp(h*fr/(kb * temperatureKelvin))-1.0));
        return  t1 * t2;
    }



    static lm::Float applySdssGBand(lm::Float blackBodyTemperature) {
        lm::Float ret = 0.0;
        for(int i = 0; i < 89; i++) {
            ret += blackbody(g_band_wavelengths_nm[i], blackBodyTemperature) * g_band_response[i];
        }
        ret /= g_band_weights_sum;
        return ret;
    }



    static lm::Float xFit_1931( lm::Float wave )
    {
        lm::Float t1 = (wave-442.0f)*((wave<442.0f)?0.0624f:0.0374f);
        lm::Float t2 = (wave-599.8f)*((wave<599.8f)?0.0264f:0.0323f);
        lm::Float t3 = (wave-501.1f)*((wave<501.1f)?0.0490f:0.0382f);
        return 0.362f*glm::exp(-0.5f*t1*t1) + 1.056f*glm::exp(-0.5f*t2*t2)
        - 0.065f*glm::exp(-0.5f*t3*t3);
    }
    static lm::Float yFit_1931( lm::Float wave )
    {
        lm::Float t1 = (wave-568.8f)*((wave<568.8f)?0.0213f:0.0247f);
        lm::Float t2 = (wave-530.9f)*((wave<530.9f)?0.0613f:0.0322f);
        return 0.821f*glm::exp(-0.5f*t1*t1) + 0.286f*glm::exp(-0.5f*t2*t2);
    }
    static lm::Float zFit_1931( lm::Float wave )
    {
        lm::Float t1 = (wave-437.0f)*((wave<437.0f)?0.0845f:0.0278f);
        lm::Float t2 = (wave-459.0f)*((wave<459.0f)?0.0385f:0.0725f);
        return 1.217f*glm::exp(-0.5f*t1*t1) + 0.681f*glm::exp(-0.5f*t2*t2);
    }

    static void tempToXYZ(lm::Float temp, lm::Float & X,lm::Float & Y,lm::Float & Z) {
        X = Y = Z = 0.0;
        for(int nmi = 340; nmi < 900; nmi++) {
            lm::Float nm = (lm::Float)nmi;
            X += xFit_1931(nm) * blackbody(nm,temp);
            Y += yFit_1931(nm) * blackbody(nm,temp);
            Z += zFit_1931(nm) * blackbody(nm,temp);
        }
    }

        



}



LM_NAMESPACE_BEGIN(LM_NAMESPACE)



typedef std::vector<std::vector<lm::Float>&> stdvec2d;
    

struct LightToCameraRaySegmentCDF {
    Vec3 weight;
    Vec3 p,d;
    lm::Float cdfSoFar;
    lm::Float localcdf;
    lm::Float t;
    lm::Float tSoFar;
    lm::Float a;
    lm::Float b;
    

};


void to_json(lm::Json& j, const LightToCameraRaySegmentCDF& p); 


void from_json(const lm::Json& j, LightToCameraRaySegmentCDF& p);


struct StarSource {
    lm::Vec3 intensity;
    lm::Vec3 position;
    int index;
};

void to_json(lm::Json& j, const StarSource& p); 


void from_json(const lm::Json& j, StarSource& p);


namespace stats {
    struct CachedSampleId {};
    struct TetraIdGuess{};

    struct SampleIdCacheHits {};
    struct SampleIdCacheMisses {};
    struct UsedCachedTetra {};
    struct UsedNeighborTetra {};
    struct InvalidNeighbor {};
    struct ResampleAccel {};
    struct TotalTetraTests {};

    struct MaxTransmittance {};
    struct FreePathTransmittance {};
    struct OpticalThickness {};

    struct RegularTrackingStrategyDistanceSample {};

    struct RegularTrackingStrategyTotalT{};//the distance until the last volume boundary
    struct RegularTrackingStrategyTotalEffT{};//the distance as above but only path segments that will participate in scatter are included
    struct RegularTrackingStrategyMinT{}; //the distance until the first volume boundary
    
    struct RegularTrackingStrategyTetraIndex{};//the tetra index the sample "landed" in

    struct RegularTrackingStrategyXi {};
    struct RegularTrackingStrategyTotalTau{}; //the optical thickness for the whole current ray
    struct RegularTrackingStrategyTauUntilScatter{}; //the optical thickness accumulated until the regular distance sample scattered
    struct RegularTrackingStrategyNormFac {};//the normalization factor to make all uniform samples scatter 
    struct RegularTrackingStrategyMuT {};//the extinction coefficient at scatter point

    struct EquiangularStrategyDistanceSample {};

    struct ScatteringAlbedo {};


    struct DistanceSampleRandomValues {};
    struct EquiDistanceSampleRandomValueVertexIndex {};
    struct RegularDistanceSampleRandomValueVertexIndex {};

    struct RegularContribution{};
    struct RegularRegularPDF{};
    struct RegularEquiPDF{};

    struct EquiContribution{};
    struct EquiEquiPDF{};
    struct EquiRegularPDF{};

    struct EmissiveContribution{};

    struct DistanceSamplesPDFs{};
    //2 distance samples, first one from equiangular, second one from regular 
    enum IJ {
        //equiangular sampling, distance sample 0
        _0_0,
        //equiangular sampling, distance sample 1
        _0_1,
        //regular sampling, distance sample 0
        _1_0,
        //regular sampling, distance sample 1
        _1_1
    };

    struct BoundaryVisitor{};
    struct LastBoundarySequence {};

    struct VRL{};
    typedef int TetraIndex;

    struct LightsInTetra{};

    struct DuplicateWatchdog {};

    




}

struct RaySegmentCDF {
    lm::Float localcdf;
    lm::Float t;
    lm::Float a;
    lm::Float b;
    lm::Float a_kelv;
    lm::Float b_kelv;
    int tetraI;
};




    
class ArepoLMMesh : public Mesh {

public:
    virtual void foreach_triangle(const ProcessTriangleFunc& process_triangle) const = 0;

    virtual lm::Mesh::Tri triangle_at(int face) const = 0;

    virtual int correspondingTetra(int face) const = 0;

    virtual const std::vector<int> & adjacentTs(int pointIndex) const = 0;


    virtual lm::Mesh::InterpolatedPoint surface_point(int face, lm::Vec2 uv) const = 0;

    virtual int num_triangles() const = 0;
};



class Volume_Arepo : public Volume {
public:
    virtual Bound bound() const = 0;
    virtual bool has_scalar() const = 0;
    virtual Float max_scalar() const = 0;
    virtual Float max_scalar(lm::Ray ray, lm::Float & out_t_forhowlong, lm::Float & out_aCoef, lm::Float & out_bCoef) const = 0;
    virtual Float eval_scalar(Vec3 p) const = 0;
    virtual Float eval_scalar(Vec3 p , Vec3 dir) const = 0;
    virtual bool has_color() const = 0;
    virtual Vec3 eval_color(Vec3 p) const = 0;
    virtual Float sample_distance(Ray ray,lm::Float tmin, lm::Float tmax, lm::Rng& rng,lm::Float & weight) const = 0;
    virtual Float eval_transmittance(lm::Ray ray, Float tmin, Float tmax) const = 0;
    virtual void visitBFS(lm::Vec3 startPos, std::function<bool(int tetraI,glm::tmat4x3<lm::Float> corners, int bfsLayer)> processor) const = 0;
    virtual int findTetra(lm::Vec3 pos) const = 0;

};


LM_NAMESPACE_END(LM_NAMESPACE)

/*
LM_NAMESPACE_BEGIN(LM_NAMESPACE::objloader)



class Mesh_Arepo : public OBJLoaderContext {
public:
    public:
    virtual bool load(
        const std::string& path,
        OBJSurfaceGeometry& geo,
        const ProcessMeshFunc& process_mesh,
        const ProcessMaterialFunc& process_material) override = 0;
      
};
LM_NAMESPACE_END(LM_NAMESPACE::objloader)

*/

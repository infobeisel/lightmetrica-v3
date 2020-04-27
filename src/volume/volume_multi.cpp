/*
	Lightmetrica - Copyright (c) 2019 Hisanari Otsu
	Distributed under MIT license. See LICENSE file for details.
*/

#include <pch.h>
#include <lm/core.h>
#include <lm/volume.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

/*
\rst
.. function:: volume::multigaussian

    Multiple Gaussian volumes.

    :param volumes: array of references to volume
\endrst
*/
class Volume_Multi : public Volume {
private:
	Bound bound_;
	std::vector< Volume* > volumes_den_;
    std::vector< Volume* > volumes_alb_;
    unsigned int size_;
    Float maxScalar_=0;
    Rng* rng;

public:
	LM_SERIALIZE_IMPL(ar) {
		ar(bound_);
	}

public:
	virtual void construct(const Json& prop) override {
		//bound_.min = json::value<Vec3>(prop, "bound_min", Vec3(Inf));
		//bound_.max = json::value<Vec3>(prop, "bound_max", Vec3(-Inf));
		auto volRefAlb = json::value< std::vector<std::string> >(prop,"volumes_alb");
        auto volRefDen = json::value< std::vector<std::string> >(prop,"volumes_den");
        if(!volRefAlb.size() || !volRefDen.size() || volRefDen.size()!=volRefAlb.size()){
            LM_THROW_EXCEPTION(Error::InvalidArgument,
			    "volumes_alb and/or volumes_den have an invalid size. They need to be of same size.");
        }
            
        size_=volRefAlb.size();
        
        for( int i=0 ; i<size_ ; i++ ){
            volumes_alb_.push_back( comp::get<Volume>(volRefAlb[i]) );
            volumes_den_.push_back( comp::get<Volume>(volRefDen[i]) );
            if(!volumes_alb_.back()->has_color())
                LM_THROW_EXCEPTION(Error::InvalidArgument,
			    "volumes_alb[{}] has no albedo/color",i);
            if(!volumes_den_.back()->has_scalar())
                LM_THROW_EXCEPTION(Error::InvalidArgument,
			    "volumes_den[{}] has no density",i);
        }
        
        // compute bounds
        Vec3 min=Vec3(Inf);
        Vec3 max=Vec3(-Inf);
        for( int i=0 ; i<size_ ; i++ ){
            const Bound b = volumes_den_[i]->bound();
            LM_INFO("min {}, {}, {}",b.min.x,b.min.y,b.min.z);
            min.x=(min.x>b.min.x)?b.min.x:min.x;
            min.y=(min.y>b.min.y)?b.min.y:min.y;
            min.z=(min.z>b.min.z)?b.min.z:min.z;

            max.x=(max.x<b.max.x)?b.max.x:max.x;
            max.y=(max.y<b.max.y)?b.max.y:max.y;
            max.z=(max.z<b.max.z)?b.max.z:max.z;
            
            maxScalar_+=volumes_den_[i]->max_scalar();
        }
        bound_.min=min;
        bound_.max=max;
    
        LM_INFO("min {}, {}, {}",bound().min.x,bound().min.y,bound().min.z);
        LM_INFO("max {}, {}, {}",bound().max.x,bound().max.y,bound().max.z);
        rng = new Rng(math::rng_seed());
	}

    Bound computeBound(const Vec3 sig, const Vec3 pos){
        auto compAC=[&](Float s) -> Float {
            return sqrt(-2._f*s*std::log(Eps));};
        Float ac = std::max(std::max(compAC(sig.x), compAC(sig.y)), compAC(sig.z));
        Bound b;
        b.max = pos+ac;
        b.min = pos-ac;
        return b;
        
    }
	
	virtual Bound bound() const override {
		return bound_;
	}

	virtual bool has_scalar() const override {
		return true;
	}

	virtual Float max_scalar() const override {
		return maxScalar_;
	}

	virtual Float eval_scalar(Vec3 p) const override {
        Float sum=0._f;
        for(auto* v : volumes_den_)
            sum+=v->eval_scalar(p);
        //for( int i=0 ; i<size_ ; i++ )
        //    sum+=volumes_[i]->eval_scalar(p);
        //return std::accumulate(volumes_.begin(),volumes_.end(),0._f,[p](Float a, Volume& b)->Float{
         //   return a+b->eval_scalar(p);});
       return sum;
	}

	virtual bool has_color() const override {
		return true;
	}

// Compute the color by randomly selecting within the different gaussians
// depending on scalar value at position p
	virtual Vec3 eval_color(Vec3 p) const override {
        std::vector<Float> scalars;
        Float sum=0;
        scalars.reserve(size_);
        // Perform a prefix sum of scalar evaluations
        // Cannot use std::accumulate since we need to keep all individual values
        for(auto* v : volumes_den_){
            sum+= v->eval_scalar(p);
            scalars.push_back(sum);
        }
        /*for( int i=0 ; i<size_ ; i++ ){
            sum += volumes_[i]->eval_scalar(p);
		    scalars.push_back( sum );
        }*/
        
        // Normalize to use as probabilities
        Float s = scalars.back();
        for(auto& sc : scalars)
            sc=sc/s;
        //for( int i=0 ; i<size_ ; i++ )
        //    scalars[i]=scalars[i]/s;
            
        // Ensures that the while loop exits
        scalars[scalars.size()-1]=1.0_f+Eps;    
        Float r = rng->u();
        Float current = 0._f;
        while(r >= scalars[current])
            current++;
        auto col = volumes_alb_[current]->eval_color(p);
        return col;
		
	}
};

LM_COMP_REG_IMPL(Volume_Multi, "volume::multi");

LM_NAMESPACE_END(LM_NAMESPACE)

#ifndef NLO_MT_VELOCIT_FUNCTOR_H_
#define NLO_MT_VELOCIT_FUNCTOR_H_

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"

#include "aliases.h"

using namespace walberla;

class VelocityFunctor
{

    const Vector3<real_t> velocity_;
    const real_t period_;
    const real_t height_;
 
    real_t constantTerm_; // part of the velocity that is constant in both time and space
    real_t amplitude_;
 public:
    VelocityFunctor(Vector3<real_t> const& velocity, real_t period, real_t height)
    : 
        velocity_(velocity),
        period_(period),
        height_(height)
    {
        constantTerm_ = real_t(4) * velocity_[0] / (height_ * height_);
    }
    
    void operator()(real_t time)
    {
        amplitude_ = constantTerm_ * real_t(0.5) * (real_t(1) - std::cos(real_t(2) * math::pi * time / period_));
    }
    
    Vector3< real_t > operator()(Vector3<real_t> const& pos, real_t)
    {
        return Vector3< real_t >(amplitude_ * pos[1] * (height_ - pos[1]), real_t(0), real_t(0));
    }
};

#endif
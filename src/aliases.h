#ifndef NLO_MT_ALIASES_H_
#define NLO_MT_ALIASES_H_

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"

using LatticeModel_T         = walberla::lbm::D3Q27<walberla::lbm::collision_model::D3Q27Cumulant>; 
using Stencil_T              = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;

using PdfField_T = walberla::lbm::PdfField< LatticeModel_T >;

using flag_t      = walberla::uint16_t;
using FlagField_T = walberla::FlagField< flag_t >;

// Communication Pack Info
using PackInfo_T =  walberla::pystencils::CumulantMRTPackInfo;

#endif
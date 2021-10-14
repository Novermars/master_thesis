#ifndef NLO_MT_ALIASES_H_
#define NLO_MT_ALIASES_H_

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"

#include "CumulantMRTNoSlip.h"
#include "CumulantMRTOutflow.h"
#include "CumulantMRTSimpleUBB.h"
#include "CumulantMRTPackInfo.h"
#include "CumulantMRTSweep.h"

using LatticeModel_T         = walberla::lbm::D3Q27<walberla::lbm::collision_model::D3Q27Cumulant>; 
using Stencil_T              = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;

using flag_t        = walberla::uint8_t;
using PdfField_T    = walberla::field::GhostLayerField<walberla::real_t, Stencil_T::Size>;
using VectorField_T = walberla::field::GhostLayerField<walberla::real_t, Stencil_T::D >;
using FlagField_T   = walberla::FlagField<flag_t>;

// Communication Pack Info
using PackInfo_T    = walberla::pystencils::CumulantMRTPackInfo;
using NoSlip_T      = walberla::lbm::CumulantMRTNoSlip;
using SimpleUBB_T   = walberla::lbm::CumulantMRTSimpleUBB;
using Outflow_T     = walberla::lbm::CumulantMRTOutflow;

using Mesh          = walberla::mesh::TriangleMesh;

#endif
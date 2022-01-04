#ifndef NLO_MT_ALIASES_H_
#define NLO_MT_ALIASES_H_

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"

#include "stencil/D3Q27.h"

#include "CumulantMRTNoSlip.h"
#include "CumulantMRTOutflow.h"
#include "CumulantMRTDynamicUBB.h"
#include "CumulantMRTPackInfo.h"
#include "CumulantMRTSweep.h"

namespace walberla 
{
    using Stencil_T = stencil::D3Q27;

    using flag_t = walberla::uint8_t;
    using PdfField_T = walberla::field::GhostLayerField<walberla::real_t, Stencil_T::Size>;
    using VectorField_T = walberla::field::GhostLayerField<walberla::real_t, Stencil_T::D>;
    using FlagField_T = walberla::FlagField<flag_t>;

    // Communication Pack Info
    using PackInfo_T = walberla::pystencils::CumulantMRTPackInfo;
    using CumulantMRTSweep = walberla::pystencils::CumulantMRTSweep;
    using NoSlip_T = walberla::lbm::CumulantMRTNoSlip;
    using DynamicUBB_T = walberla::lbm::CumulantMRTDynamicUBB;
    using Outflow_T = walberla::lbm::CumulantMRTOutflow;

    using Mesh = walberla::mesh::TriangleMesh;
}

// Flags for the boundary and fluid cells
const walberla::FlagUID FluidFlagUID("Fluid Flag");
const walberla::FlagUID NoSlipFlagUID("NoSlip Flag");
const walberla::FlagUID OutflowUID("Outflow Flag");
const walberla::FlagUID InflowUID("Inflow Flag");

// UIDs for the boundary and fluid cells
const walberla::BoundaryUID NoSlipBoundaryUID("NoSlip Boundary");
const walberla::BoundaryUID OutflowBoundaryUID("Outflow Boundary");
const walberla::BoundaryUID InflowBoundaryUID("Inflow Boundary");

// Standarized coloring for the mesh
const walberla::Mesh::Color noSlipColor{255, 255, 255}; // White
const walberla::Mesh::Color inflowColor{255, 0, 0};     // Red  
const walberla::Mesh::Color outflowColor{0, 0, 255};    // Blue

// Number of ghost layers
walberla::uint_t const numGhostLayers = walberla::uint_t(1);

#endif
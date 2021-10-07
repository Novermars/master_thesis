#ifndef NLO_MT_INFLOWBC_H_
#define NLO_MT_INFLOWBC_H_

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"

#include "VelocityFunctor.h"
#include "aliases.h"

using namespace walberla;

/*
 * Class to implement a dynamic inflow boundary for the aneurysm.
 * We assume that at the inflow, we have a fully developed periodic Poiseulle flow
 * simulating the pulsating blood flow
 * 
 * This code is very strongly influenced by tutorial 06 on boundary conditions
 * and the code written by Daniel Bauer for his bachelor thesis.
 */

// Number of ghost layers
uint_t const FieldGhostLayers = uint_t(1);

/*
 * We have four different kind of cells:
 * 1) Fluid cells
 * 2) Boundary cells with no-slip behavior
 * 3) Inflow boundary cells where we have a velocity boundary condition
 * 4) Outflow boundary cells where we have an extrapolated pressure boundary conditon
 */
FlagUID const FluidFlagUID("Fluid Flag");
FlagUID const NoSlipFlagUID("NoSlip Flag");
FlagUID const SimpleUBBFlagUID("SimpleUBB Flag");
FlagUID const DynamicUBBFlagUID("DynamicUBB Flag");
FlagUID const OutletFlagUID("Outlet Flag");

using NoSlip_T           = lbm::NoSlip<LatticeModel_T, flag_t>;
using SimpleUBB_T        = lbm::SimpleUBB< LatticeModel_T, flag_t >;
using DynamicUBB_T       = lbm::DynamicUBB<LatticeModel_T, flag_t, VelocityFunctor>;
using Outlet_T           = lbm::Outlet<LatticeModel_T, FlagField_T>;

using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, NoSlip_T, SimpleUBB_T, Outlet_T>;
//using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, NoSlip_T, DynamicUBB_T, Outlet_T>;

/*
 * Class that implements our custom boundary handling
 */

class MyBoundaryHandling
{
    BlockDataID const flagFieldID_;
    BlockDataID const pdfFieldID_;

    real_t period_;
    Vector3<real_t> maxInflowVelocity_;
    std::shared_ptr<lbm::TimeTracker> timeTracker_;
public:
    MyBoundaryHandling(const BlockDataID& flagFieldID, const BlockDataID& pdfFieldID, real_t period,
                        Vector3<real_t> maxInflowVelocity, const std::shared_ptr< lbm::TimeTracker >& timeTracker)
    : 
        flagFieldID_(flagFieldID), 
        pdfFieldID_(pdfFieldID), 
        period_(period), 
        maxInflowVelocity_(maxInflowVelocity),
        timeTracker_(timeTracker)
    {}

    BoundaryHandling_T* operator()(IBlock* const block, StructuredBlockStorage const* const storage) const;
};

BoundaryHandling_T* MyBoundaryHandling::operator()(IBlock* const block,
                                                   const StructuredBlockStorage* const storage) const
{
    Vector3<real_t> domainSize(real_c(storage->getNumberOfXCells()), real_c(storage->getNumberOfYCells()),
                                    real_c(storage->getNumberOfZCells()));

    real_t height = domainSize[1];

    //VelocityFunctor velocity(maxInflowVelocity_, period_, height);

    WALBERLA_ASSERT_NOT_NULLPTR(block)

    FlagField_T* flagField = block->getData< FlagField_T >(flagFieldID_);
    PdfField_T* pdfField   = block->getData< PdfField_T >(pdfFieldID_);

    auto const fluidFlag = flagField->getOrRegisterFlag(FluidFlagUID);

    BoundaryHandling_T* handling = new BoundaryHandling_T(
            "Boundary Handling", flagField, fluidFlag,
            NoSlip_T("NoSlip", NoSlipFlagUID, pdfField),
            SimpleUBB_T("SimpleUBB", SimpleUBBFlagUID, pdfField, maxInflowVelocity_),
            /*DynamicUBB_T("DynamicUBB", DynamicUBBFlagUID, pdfField, timeTracker_, storage->getLevel(*block), velocity,
                        block->getAABB()),*/
            Outlet_T("Outlet", OutletFlagUID, pdfField, flagField, fluidFlag)
        );

    CellInterval domainBB = storage->getDomainCellBB();
    storage->transformGlobalToBlockLocalCellInterval(domainBB, *block);

    /*
    handling->forceBoundary(DynamicUBBFlagUID, west);
    handling->forceBoundary(OutletFlagUID, east);
    handling->forceBoundary(NoSlipFlagUID, south);
    */
    handling->fillWithDomain(domainBB);

    return handling;
}

#endif
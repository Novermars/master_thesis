#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"

#include "mesh_common/DistanceComputations.h"
#include "mesh_common/DistanceFunction.h"
#include "mesh_common/MatrixVectorOperations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/vtk/CommonDataSources.h"
#include "mesh_common/vtk/VTKMeshWriter.h"
#include "mesh/blockforest/BlockExclusion.h"
#include "mesh/blockforest/BlockForestInitialization.h"
#include "mesh/boundary/BoundaryInfo.h"
#include "mesh/boundary/BoundaryLocation.h"
#include "mesh/boundary/BoundaryLocationFunction.h"
#include "mesh/boundary/BoundarySetup.h"
#include "mesh/boundary/BoundaryUIDFaceDataSource.h"
#include "mesh/boundary/ColorToBoundaryMapper.h"
#include "mesh_common/MeshOperations.h"


#include "InitialPDFsSetter.h"
#include "aliases.h"
#include "Parameters.h"
#include "circle.h"

#include "json.hpp"
#include "npy.hpp"
#include "spline.h"

#include <iostream>
#include <cmath>

using namespace walberla;

namespace NLO
{
// Hash function for std::unordered_map<NLO::Circle::Coords, T>
struct key_hash : public std::unary_function<NLO::Circle::Coords, std::size_t>
{
    std::size_t operator()(const NLO::Circle::Coords& coord) const
    {
        // Kindly stolen from waLBerla
        Cell cell{std::get<0>(coord), std::get<1>(coord), std::get<2>(coord)};
        std::size_t seed = 0;
        std::hash<cell_idx_t> hasher;

        seed ^= hasher(cell.x()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(cell.y()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(cell.z()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};

// Function that returns the a cublic spline interpolator to find the
// missing values from Andreas' data
tk::spline getInterpolater()
{
    nlohmann::json heartBeatData_json;
    std::ifstream heartBeatFile{"cleanHBData.json"};
    heartBeatFile >> heartBeatData_json;

    auto t_vals = heartBeatData_json["t"].get<std::vector<real_t>>();
    auto v_vals = heartBeatData_json["v"].get<std::vector<real_t>>();

    tk::spline interpolater(t_vals, v_vals);
    return interpolater;
}

auto determineInputCells(std::shared_ptr<StructuredBlockForest> blocks, BlockDataID flagFieldId)
{
        // We want to calculate the radius and the middle point of the circle that describes the inflow
        // For that we loop over all blocks in the SbF on the rank, and then in every block we loop over 
        // all the cells without the ghost nodes and check if the cell is part of the inflow boundary
        // If that's the case, we save its global coordinates 
        std::vector<Block const*> block_vec;
        int level = 0;
        blocks->getBlocks(block_vec, level);
        std::vector<int> inflowX;
        std::vector<int> inflowY;
        std::vector<int> inflowZ;
        for (auto const block: block_vec)
        {
            auto* flagField = block->getData<FlagField_T>(flagFieldId);
            flag_t inflowFlag = flagField->getFlag(InflowUID);
            
            for (cell_idx_t x = 0; x < static_cast<int>(flagField->xSize()); ++x)
            {
                for (cell_idx_t y = 0; y < static_cast<int>(flagField->ySize()); ++y)
                {
                    for (cell_idx_t z = 0; z < static_cast<int>(flagField->zSize()); ++z)
                    {
                        if (flagField->isFlagSet(x, y, z, inflowFlag))
                        {
                            Cell cell{x, y, z};
                            blocks->transformBlockLocalToGlobalCell(cell, *block);
                            inflowX.push_back(static_cast<int>(cell.x()));
                            inflowY.push_back(static_cast<int>(cell.y()));
                            inflowZ.push_back(static_cast<int>(cell.z()));
                        }
                    }
                }
            }
        }
        return std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>({inflowX, inflowY, inflowZ});
}

void writeMetaData(Parameters parameters, std::vector<int> const& inflowX0,
                   std::vector<int> const& inflowY0, std::vector<int> const& inflowZ0 )
{
    std::vector<NLO::Circle::Coords> inflowCells;
    NLO::Circle::Coords middle;
    uint_t radius;

    std::string outputFile{"cellOutput.json"};
    std::ofstream output_stream{outputFile};
    nlohmann::json output_json;

    output_json["timesteps"] = parameters.timeSteps_;
    output_json["numHeartBeats"] = parameters.numHeartBeats_;
    output_json["dx"] = parameters.dx_;
    output_json["omega"] = parameters.omega_;
    output_json["numConstNoises"] = parameters.numConstNoises_;
    output_json["numConstInflow"] = parameters.numConstInflow_;

    for (std::size_t idx = 0; idx != inflowX0.size(); ++idx)
    {
        inflowCells.push_back({inflowX0[idx], inflowY0[idx], inflowZ0[idx]});
    }

    NLO::Circle circle = NLO::Circle(inflowCells);
    middle = circle.middle();
    radius = circle.radius();
    auto radiusY = circle.radiusY();   
    auto radiusZ = circle.radiusZ();

    output_json["diameter"] = 2 * radius;
    output_json["middle"] = {std::get<0>(middle), std::get<1>(middle), std::get<2>(middle)};
    output_json["radiusY"] = radiusY;
    output_json["radiusZ"] = radiusZ;

    output_stream << output_json;    
}

} //end of namespace NLO

class InflowProfile
{
    NLO::Circle::Coords middle_;
    uint_t radius_;
    uint_t radiusSq_;
    uint_t radiusY_;
    uint_t radiusZ_;
    real_t dt_;
    real_t dx_;
    std::vector<real_t>* noise_values_;
    real_t gamma_;
    std::shared_ptr<lbm::TimeTracker> timeTracker_;
    real_t* centerVelocity_;
    tk::spline interpolated_;

public:
    InflowProfile(NLO::Circle::Coords const& middle, uint_t radiusY, uint_t radiusZ, 
                  real_t dt, real_t dx, real_t gamma,
                  std::shared_ptr<lbm::TimeTracker> timeTracker,
                  std::vector<real_t>* noise_values,
                  real_t* centerVelocity) 
    :
        middle_{middle},
        radiusY_{radiusY},
        radiusZ_{radiusZ},
        dt_{dt},
        dx_{dx},
        noise_values_{noise_values},
        gamma_{gamma},
        timeTracker_{timeTracker}
    {        
        radius_ = std::max(radiusY_, radiusZ_);
        radiusSq_ = radius_ * radius_;

        // load noise values from the first file
        std::string path = "data/noise_field_0.npy";
        centerVelocity_ = centerVelocity;
        std::vector<unsigned long> shape;
        bool fortran_order = false;
        npy::LoadArrayFromNumpy<double>(path, shape, fortran_order, *noise_values_);
    }

    Vector3< real_t > operator()( const Cell& pos, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block ) {
        auto inflowCell = pos;
        SbF->transformBlockLocalToGlobalCell(inflowCell, block);
        real_t distanceSq = std::pow(inflowCell.x() - std::get<0>(middle_), 2) +
                            std::pow(inflowCell.y() - std::get<1>(middle_), 2) +
                            std::pow(inflowCell.z() - std::get<2>(middle_), 2);
        real_t scaleFactor = (gamma_ + 2) / gamma_ * (1 - std::pow(distanceSq / static_cast<real_t>(radiusSq_), gamma_));
        scaleFactor = std::max(scaleFactor, 0.0);
        
        // Noise is scaled to have max_val = 1 for all time values
        auto noiseVector = findNoise(inflowCell);
        // Assume 10% noise of the center velocity in all directions
        Vector3<real_t> noiseScaling = Vector3<real_t>{0.1, 0.1, 0.1 };
        for (int idx = 0; idx != 3; ++idx)
            noiseVector[idx] = *centerVelocity_ * noiseScaling[idx] * noiseVector[idx];
        
        // Main flow is in the -x direction, other directions only have noise
        Vector3<real_t> velocity{-*centerVelocity_, 0, 0};
        return scaleFactor * (dt_ / dx_) * (velocity - noiseVector);
    }

private:
    // Returns the noise value corresponding to the cell's location in the inflow bdy
    // Assumes that we have a constant x value
    Vector3<real_t> findNoise(Cell const& globalCell)
    {
        // We want to find noise_values_[y, z, {0,1,2}]
        int yCoord = std::abs(globalCell.y() - std::get<1>(middle_));
        int zCoord = std::abs(globalCell.z() - std::get<2>(middle_));
        real_t x = (*noise_values_)[getIndex(yCoord, zCoord, 0)];
        real_t y = (*noise_values_)[getIndex(yCoord, zCoord, 1)];
        real_t z = (*noise_values_)[getIndex(yCoord, zCoord, 2)];
        return Vector3<real_t>{x, y, z};
    }

    // Converts 3D index to 1D index in flat array
    int getIndex(int yCoord, int zCoord, int direction)
    {
        int diameter = static_cast<int>(radius_ * 2) + 1;
        return 3 * yCoord * diameter + 3 * zCoord + direction;
    }

}; 

/*
 * waLBerla expects that the faces are colored, not the vertices
 * So in this function, we convert our vertex-colored mesh to a face colored mesh.
 * If not all vertices corresponding to a face have the same color,
 * then we give it the default color (white/NoSlipBoundary)
 */
void vertexToFaceColor(Mesh& mesh, const Mesh::Color& defaultColor)
{
    WALBERLA_CHECK(mesh.has_vertex_colors())
    mesh.request_face_colors();

    for (auto faceIt = mesh.faces_begin(); faceIt != mesh.faces_end(); ++faceIt)
    {
        Mesh::Color vertexColor;

        bool useVertexColor = true;

        auto vertexIt = mesh.fv_iter(*faceIt);
        WALBERLA_ASSERT(vertexIt.is_valid())

        // Take the color of the first vertex to be the color for this face
        vertexColor = mesh.color(*vertexIt);

        ++vertexIt;
        // If we have a vertex in this face which has a different color,
        // then we give the face the default color.
        while (vertexIt.is_valid() && useVertexColor)
        {
            if (vertexColor != mesh.color(*vertexIt)) useVertexColor = false;
            ++vertexIt;
        }

        mesh.set_color(*faceIt, useVertexColor ? vertexColor : defaultColor);
    }
}

/* Helper function to construct a Parameters object
 * All field have public visibility and hence can be adressed
 * directly
 */

Parameters constructParameters(std::shared_ptr<Config> config)
{
    auto parameters = config->getOneBlock("Parameters");

    real_t omega = parameters.getParameter< real_t >("omega", real_c(1.8));
    uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(1000));

    double remainingTimeLoggerFrequency =
        parameters.getParameter< double >("remainingTimeLoggerFrequency", 3.0); // in seconds
    const uint_t VTKwriteFrequency = parameters.getParameter< uint_t >("VTKwriteFrequency", 10);
    real_t numHeartBeats = parameters.getParameter<real_t>("numHeartBeatCycles", 1);
    uint_t numConstNoises = parameters.getParameter<uint_t>("numConstNoises", 100);
    bool generateInflowProfile = parameters.getParameter<bool>("generateInflowProfile", 0);
    uint_t numConstInflow = parameters.getParameter<uint_t>("numConstInflow", 10);
    std::string vtkPrefixPath = parameters.getParameter< std::string >("vtkPathPrefix");

    // read domain parameters
    auto domainParameters = config->getOneBlock("DomainSetup");
    std::string meshFile = domainParameters.getParameter< std::string >("meshFile");

    real_t const dx = domainParameters.getParameter< real_t >("dx", real_t(1));
    Vector3<uint_t> const cellsPerBlock = domainParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");

    auto stabilityCheckerParam = config->getOneBlock("StabilityChecker");
    uint_t stabilityChecker = stabilityCheckerParam.getParameter<uint_t>("StabilityChecker", 1000);
    return Parameters(dx, timesteps, VTKwriteFrequency, stabilityChecker, numConstNoises, cellsPerBlock,
                      remainingTimeLoggerFrequency, numHeartBeats, meshFile, omega, generateInflowProfile,
                      numConstInflow, vtkPrefixPath);
}

/* Loads the mesh from the file including the nonstandard vertex colors
 * Then it is broadcasted over the possible MPI processes
 * Finally, we convert the vertex colors to face colors as waLBerla
 * expects it
 */
std::shared_ptr<Mesh> readAndColorMesh(const std::string& fileName)
{
    // read in mesh with vertex colors on a single process and broadcast it
    auto mesh = make_shared<Mesh>();
    mesh->request_vertex_colors();
    mesh::readAndBroadcast(fileName, *mesh);

    // color faces according to vertices
    vertexToFaceColor(*mesh, noSlipColor);

    return mesh;
}

shared_ptr<mesh::DistanceOctree<Mesh>> buildDistanceOctree(shared_ptr<Mesh> mesh)
{
    // add information to mesh that is required for computing signed distances from a point to a triangle
    auto triDist = make_shared<mesh::TriangleDistance<Mesh>>(mesh);
    return make_shared<mesh::DistanceOctree<Mesh>>(triDist);
}


std::shared_ptr<StructuredBlockForest> createBlockForest(Parameters const& params,
                                                         AABB const& aabb,
                                                         std::shared_ptr<mesh::DistanceOctree<Mesh>> distanceOctree
                                                        )
{
  mesh::ComplexGeometryStructuredBlockforestCreator bfc(aabb, Vector3<real_t>(real_c(1)),
                                                        mesh::makeExcludeMeshExterior(distanceOctree, real_c(1))
                                                       );
    bfc.setPeriodicity(Vector3<bool>(false));

    return bfc.createStructuredBlockForest(params.cellsPerBlock_);
}

int main(int argc, char* argv[])
{
    walberla::Environment walberlaEnv(argc, argv);
    // We get weird MPI errors without this line
    mpi::MPIManager::instance()->useWorldComm();
    auto config = walberlaEnv.config();
    auto parameters = constructParameters(config);
    std::string vtkPathPrefix = parameters.vtkPathPrefix_;
    auto mesh = readAndColorMesh(parameters.meshFileName_);

    // REMARK: scale mesh with 1/dx
    typename Mesh::Scalar meshScaleFactor = real_c(1) / parameters.dx_;
    mesh::scale(*mesh, Vector3<typename Mesh::Scalar>(meshScaleFactor));

    auto distanceOctree = buildDistanceOctree(mesh);
    auto aabb = computeAABB(*mesh);
    aabb.scale(real_c(1.2)); // increase aabb to make sure that boundary conditions are far away from domain boundaries
    auto blocks = createBlockForest(parameters, aabb, distanceOctree);
    WALBERLA_LOG_DEVEL_VAR(blocks->getNumberOfBlocks());

    // Let the root process write the calculated octree to a file
    WALBERLA_ROOT_SECTION()
    {
        distanceOctree->writeVTKOutput(vtkPathPrefix + "distanceOctree");
    }

    // In this LBM simulation, we need three fields, which are very similar to (multi)dimensional arrays
    // One for the macroscopic velocity values
    BlockDataID velocityFieldId = field::addToStorage< VectorField_T >(blocks, "velocity", real_c(0), field::fzyx);
    // One for the flags which determine if a cell is part of the fluid or boundary
    BlockDataID flagFieldId     = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");
    // And one for the mesoscopic particule distribution function values;
    BlockDataID pdfFieldId      = field::addToStorage< PdfField_T >(blocks, "pdf field", real_c(0.0), field::fzyx);

    // Set the values of both the tmp and real field to a default value for the density rho
    real_t rho = real_c(1);
    pystencils::InitialPDFsSetter pdfSetter(pdfFieldId, velocityFieldId, rho);

    // Does one sweep of the lbm algorithm for the fluid cells
    pystencils::CumulantMRTSweep cumulantMRTSweep(pdfFieldId, velocityFieldId, parameters.omega_);
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
    {
        pdfSetter(&(*blockIt));
    }

    // register flags at the flag field
    for (auto blockIterator = blocks->begin(); blockIterator != blocks->end(); ++blockIterator) {
        FlagField_T *flagField = blockIterator->getData<FlagField_T>(flagFieldId);
        flagField->registerFlag(FluidFlagUID);
        flagField->registerFlag(NoSlipFlagUID);
        flagField->registerFlag(OutflowUID);
        flagField->registerFlag(InflowUID);
    }

    // map boundaryUIDs to colored mesh faces
    mesh::ColorToBoundaryMapper<mesh::TriangleMesh> colorToBoundaryMapper =
            mesh::ColorToBoundaryMapper<mesh::TriangleMesh>(mesh::BoundaryInfo(NoSlipBoundaryUID));

    colorToBoundaryMapper.set(noSlipColor, mesh::BoundaryInfo(NoSlipBoundaryUID));
    colorToBoundaryMapper.set(outflowColor, mesh::BoundaryInfo(OutflowBoundaryUID));
    colorToBoundaryMapper.set(inflowColor, mesh::BoundaryInfo(InflowBoundaryUID));

    auto boundaryLocations = colorToBoundaryMapper.addBoundaryInfoToMesh(*mesh);

    // print mesh with mapped boundaryUIDs
    mesh::VTKMeshWriter<mesh::TriangleMesh> meshWriter(mesh, "meshBoundaries", 1);
    meshWriter.addDataSource(make_shared<mesh::BoundaryUIDFaceDataSource<mesh::TriangleMesh> >(boundaryLocations));
    meshWriter.addDataSource(make_shared<mesh::ColorFaceDataSource<mesh::TriangleMesh> >());
    meshWriter.addDataSource(make_shared<mesh::ColorVertexDataSource<mesh::TriangleMesh> >());
    meshWriter();

    mesh::BoundarySetup boundarySetup(blocks, makeMeshDistanceFunction(distanceOctree), numGhostLayers);

    // set whole region inside the mesh to fluid
    boundarySetup.setFlag<FlagField_T>(flagFieldId, FluidFlagUID, mesh::BoundarySetup::INSIDE);

    // set whole region outside the mesh to no-slip
    boundarySetup.setFlag<FlagField_T>(flagFieldId, NoSlipFlagUID, mesh::BoundarySetup::OUTSIDE);

    // set outflow flag to outflow boundary
    boundarySetup.setBoundaryFlag<FlagField_T>(flagFieldId, OutflowUID, OutflowBoundaryUID,
                                               makeBoundaryLocationFunction(distanceOctree, boundaryLocations),
                                               mesh::BoundarySetup::OUTSIDE);

    // set inflow flag to inflow boundary
    boundarySetup.setBoundaryFlag<FlagField_T>(flagFieldId, InflowUID, InflowBoundaryUID,
                                               makeBoundaryLocationFunction(distanceOctree, boundaryLocations),
                                               mesh::BoundarySetup::OUTSIDE);

    // Only generate new inflow profile if specified
    // Warning: Must have proper amount of files and the dx and mesh should be the same!
    // Only useful for debugging
    if (parameters.generateInflowProfile_)
    {
        WALBERLA_ROOT_SECTION()
        {
            WALBERLA_LOG_DEVEL("Generating new inflow profile!");
        }

        // fancy structured binding 
        auto [inflowX, inflowY, inflowZ] = NLO::determineInputCells(blocks, flagFieldId);

        // We then send the global coordinates of the cells belonging to the inflow boundary
        // to rank zero
        std::vector<int> inflowX0 = mpi::gatherv(inflowX, 0, MPI_COMM_WORLD);
        std::vector<int> inflowY0 = mpi::gatherv(inflowY, 0, MPI_COMM_WORLD);
        std::vector<int> inflowZ0 = mpi::gatherv(inflowZ, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD); // make sure that all the ranks are synchronised here
        if (mpi::MPIManager::instance()->rank() == 0)
        {
            // 
            NLO::writeMetaData(parameters, inflowX0, inflowY0, inflowZ0);

            // Now call the Python program which generates the noise array
            if (std::system("python3 generate_noise_array.py") == -1)
            {
                WALBERLA_LOG_DEVEL_VAR("Error in Python script, exiting!")
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
        }
    }    

    // Make sure that we finish the preprocessing before we continue
    MPI_Barrier(MPI_COMM_WORLD);

    nlohmann::json metaData_json;
    std::ifstream metaDataFile{"metaData.json"};
    metaDataFile >> metaData_json;
    auto timeSteps = metaData_json["timesteps"];
    real_t dt = metaData_json["dt"];
    real_t dx = metaData_json["dx"];
    uint_t radiusY = metaData_json["radiusY"];
    uint_t radiusZ = metaData_json["radiusZ"];
    NLO::Circle::Coords middle = {metaData_json["middle"][0], metaData_json["middle"][1], metaData_json["middle"][2]};
    real_t gamma = 9;

    auto interpolater = NLO::getInterpolater();
    real_t centerVelocity = interpolater(0);

    std::shared_ptr< lbm::TimeTracker > timeTracker = std::make_shared< lbm::TimeTracker >();
    std::vector<real_t> noiseValues;
    std::function< Vector3< real_t >(const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
        inflowProfile = InflowProfile(middle, radiusY, radiusZ, dt, dx, gamma, timeTracker, &noiseValues, &centerVelocity);

    NoSlip_T noSlip(blocks, pdfFieldId);
    DynamicUBB_T dynamicUBB(blocks, pdfFieldId, inflowProfile);
    Outflow_T outflow(blocks, pdfFieldId, real_c(1.0));

    noSlip.fillFromFlagField<FlagField_T>(blocks, flagFieldId, NoSlipFlagUID, FluidFlagUID);
    dynamicUBB.fillFromFlagField<FlagField_T>(blocks, flagFieldId, InflowUID, FluidFlagUID);
    outflow.fillFromFlagField<FlagField_T>(blocks, flagFieldId, OutflowUID, FluidFlagUID);

    WALBERLA_ROOT_SECTION()
    {
        WALBERLA_LOG_DEVEL_VAR(timeSteps);
    }
    
    SweepTimeloop timeloop(blocks->getBlockStorage(), timeSteps);

    blockforest::communication::UniformBufferedScheme<Stencil_T> communication(blocks);
    communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));

    // Add the partial sweeps to the timeloop
    timeloop.add() << BeforeFunction(communication, "communication") << Sweep(noSlip);
    timeloop.add() << Sweep(dynamicUBB);
    timeloop.add() << Sweep(outflow);
    timeloop.add() << Sweep(cumulantMRTSweep);

    // Time logger
    timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), parameters.remainingTimeLoggerFrequency_),
                                 "remaining time logger");
    timeloop.addFuncAfterTimeStep(makeSharedFunctor(timeTracker), "time tracking");

    // This is a hack which forces the inflow boundary values to be determined at every timestep
    auto updateBdyValues = [&](){
        dynamicUBB.fillFromFlagField<FlagField_T>(blocks, flagFieldId, InflowUID, FluidFlagUID);
    };

    auto updateNoiseValues = [&](){
        // Update only every numConstNoises steps
        auto time = static_cast<int>(timeTracker->getTime());
        if (time % parameters.numConstNoises_ == 0)
        {
            int noiseIdx = static_cast<int>( time / parameters.numConstNoises_);
            std::string path = "data/noise_field_" + std::to_string(noiseIdx) + ".npy";
            std::vector<unsigned long> shape;
            bool fortran_order = false;
            npy::LoadArrayFromNumpy<double>(path, shape, fortran_order, noiseValues);
        }
        // Calculate the velocity at the center
        // We map the time values to [0, 1) due to periodicity 
        centerVelocity = interpolater(std::fmod(time * dt, 1.0));
    };

    timeloop.addFuncAfterTimeStep(updateNoiseValues, "update noise values");
    timeloop.addFuncAfterTimeStep(updateBdyValues, "update bdy values");

    // Set velocity in boundary cells to zero for better output
    auto zeroSetterFunction = [&](IBlock* block) {
        VectorField_T* velocityField   = block->getData< VectorField_T >(velocityFieldId);
        FlagField_T* flagField = block->getData< FlagField_T >(flagFieldId);
        flag_t fluidFlag = flagField->getFlag(FluidFlagUID);

        WALBERLA_FOR_ALL_CELLS(
                velIt, velocityField, flagIt, flagField,
                if (*flagIt != fluidFlag) {
                    velIt[0] = real_c(0);
                    velIt[1] = real_c(0);
                    velIt[2] = real_c(0);
                }
                ) // WALBERLA_FOR_ALL_CELLS
        };
    timeloop.add() << Sweep(zeroSetterFunction, "Zerosetter");
    
    timeloop.addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker<VectorField_T>(
	                                blocks, velocityFieldId, parameters.stabilityCheckFrequency_)));

    if (parameters.vtkWriteFrequency_ > 0)
    {
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "cumulant_mrt_velocity_field", parameters.vtkWriteFrequency_, 0,
                                                        true, vtkPathPrefix, "simulation_step", false, true, true, false, 0);
	    // Set the sampling size for the output
	    // The output size is divided by resolution in all directions
	    // So the final size is resolution^3 times smaller
	    constexpr double resolution = 2;
	    static_assert(resolution >= 1, "You probably want to have this number be at least 1");
        vtkOutput->setSamplingResolution(resolution);

        auto velWriter = make_shared< field::VTKWriter< VectorField_T > >(velocityFieldId, "Velocity");
        vtkOutput->addCellDataWriter(velWriter);

        timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
    }

    // write flag field
    auto vtkFlagField = field::createVTKOutput<FlagField_T>(flagFieldId, *blocks, "flag_field", uint_c(1), uint_c(0));
    vtkFlagField();

    // write domain decomposition
    vtk::writeDomainDecomposition(blocks);
    timeloop.run();
}

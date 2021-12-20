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

#include <iostream>
#include <cmath>

//#include <Python.h>

using namespace walberla;

FlagUID const FluidFlagUID("Fluid Flag");
FlagUID const NoSlipFlagUID("NoSlip Flag");
FlagUID const OutflowUID("Outflow Flag");
FlagUID const InflowUID("Inflow Flag");

BoundaryUID const NoSlipBoundaryUID("NoSlip Boundary");
BoundaryUID const OutflowBoundaryUID("Outflow Boundary");
BoundaryUID const InflowBoundaryUID("Inflow Boundary");

namespace NLO
{
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
} //end of namespace NLO

class InflowProfile
{
public:

    InflowProfile(std::shared_ptr<std::unordered_map<NLO::Circle::Coords, 
                                                     Vector3<real_t>, NLO::key_hash>> const& values_map) 
    :
        values_map_(values_map)
    {        
        std::string fileName{"data/inflow_profile_0.json"};
        std::ifstream file_streamNoise(fileName);
        nlohmann::json values;
        file_streamNoise >> values;
        for (std::size_t idx = 0; idx != values.size(); ++idx)
        {
            values_map_->insert({{static_cast<int>(values[idx][0]), 
                                  static_cast<int>(values[idx][1]), 
                                  static_cast<int>(values[idx][2])}, 
                                  Vector3<real_t>{values[idx][3], 
                                                  values[idx][4], 
                                                  values[idx][5]}});
        }
    }

    Vector3< real_t > operator()( const Cell& pos, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block ) {
        auto inflowCell = pos;
        SbF->transformBlockLocalToGlobalCell(inflowCell, block);
        return (*values_map_)[NLO::Circle::Coords{inflowCell.x(), inflowCell.y(), inflowCell.z()}];
    }
private:
    std::shared_ptr<std::unordered_map<NLO::Circle::Coords, Vector3<real_t>, NLO::key_hash>> values_map_;
}; 

/*
// AN181
std::string modelName = "AN181";
constexpr Mesh::Color noSlipColor{255, 255, 255};
constexpr Mesh::Color inflowColor{255, 0, 12};
constexpr Mesh::Color outflowColor{0, 42, 255};
*/

// AN182
//std::string modelName = "AN182";
constexpr Mesh::Color noSlipColor{255, 255, 255};
constexpr Mesh::Color inflowColor{255, 0, 42};
constexpr Mesh::Color outflowColor{0, 98, 255};

// Number of ghost layers
uint_t const numGhostLayers = uint_t(1);

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
    bool generateInflowProfile = parameters.getParameter<uint_t>("generateInflowProfile", 0) == 0;

    // read domain parameters
    auto domainParameters = config->getOneBlock("DomainSetup");
    std::string meshFile = domainParameters.getParameter< std::string >("meshFile");

    real_t const dx = domainParameters.getParameter< real_t >("dx", real_t(1));
    Vector3<uint_t> const cellsPerBlock = domainParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");

    auto stabilityCheckerParam = config->getOneBlock("StabilityChecker");
    uint_t stabilityChecker = stabilityCheckerParam.getParameter<uint_t>("StabilityChecker", 1000);
    return Parameters(dx, timesteps, VTKwriteFrequency, stabilityChecker, numConstNoises, cellsPerBlock,
                      remainingTimeLoggerFrequency, numHeartBeats, meshFile, omega, generateInflowProfile);
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
  /// REMARK: we usually leave dx here at 1 and scale the mesh instead (see below)
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
        distanceOctree->writeVTKOutput("vtk_out/distanceOctree");
    }

    // In this LBM simulation, we need three fields, which are very similar to (multi)dimensional arrays
    // One for the macroscopic velocity values
    BlockDataID velocityFieldId = field::addToStorage< VectorField_T >(blocks, "velocity", real_c(0), field::fzyx);
    // One for the flags which determine if a cell is part of the fluid or boundary
    BlockDataID flagFieldId     = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");
    // And one for the mesoscopic particule distribution function values;
    BlockDataID pdfFieldId      = field::addToStorage< PdfField_T >(blocks, "pdf field", real_c(0.0), field::fzyx);

    /// REMARK: this must not be skipped; there are two pdfFields (incl. a temporary one for pointer swapping) and both
    /// must have meaningful initial values
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
        // We then send the global coordinates of the cells belonging to the inflow boundary
        // to rank zero
        std::vector<int> inflowX0 = mpi::gatherv(inflowX, 0, MPI_COMM_WORLD);
        std::vector<int> inflowY0 = mpi::gatherv(inflowY, 0, MPI_COMM_WORLD);
        std::vector<int> inflowZ0 = mpi::gatherv(inflowZ, 0, MPI_COMM_WORLD);

        std::vector<NLO::Circle::Coords> inflowCells;
        NLO::Circle::Coords middle;
        real_t radius;

        MPI_Barrier(MPI_COMM_WORLD); // make sure that all the ranks are synchronised here
        if (mpi::MPIManager::instance()->rank() == 0)
        {
            { // Makes sure file is closed before running python script
            std::string outputFile{"cellOutput.json"};
            std::ofstream output_stream{outputFile};
            nlohmann::json output_json;
            output_json["timesteps"] = parameters.timeSteps_;
            output_json["numHeartBeats"] = parameters.numHeartBeats_;
            output_json["dx"] = parameters.dx_;
            output_json["omega"] = parameters.omega_;
            output_json["numConstNoises"] = parameters.numConstNoises_;
            for (std::size_t idx = 0; idx != inflowX0.size(); ++idx)
            {
                inflowCells.push_back({inflowX0[idx], inflowY0[idx], inflowZ0[idx]});
                // Put the cell's coordinates in the json file for the python program to read it
                output_json["cells"][idx] = {inflowX0[idx], inflowY0[idx], inflowZ0[idx]};
            }

            NLO::Circle circle = NLO::Circle(inflowCells);
            middle = circle.middle();
            radius = circle.radius();   
            output_json["diameter"] = 2 * radius;
            output_json["middle"] = {std::get<0>(middle), std::get<1>(middle), std::get<2>(middle)};
            output_stream << output_json;
            }

            // Now call the Python program which generates the inflow profile including the noise
            if (std::system("python3 generate_noise.py") == -1)
            {
                WALBERLA_LOG_DEVEL_VAR("Error in Python script, exiting!")
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
        }
    }    
    // Make sure that we finish the preprocessing before we continue
    MPI_Barrier(MPI_COMM_WORLD);

    std::shared_ptr< lbm::TimeTracker > timeTracker = std::make_shared< lbm::TimeTracker >();
    std::shared_ptr<std::unordered_map<NLO::Circle::Coords, Vector3<real_t>, NLO::key_hash>> noise 
        = std::make_shared<std::unordered_map<NLO::Circle::Coords, Vector3<real_t>, NLO::key_hash>>();
    std::function< Vector3< real_t >(const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
        inflowProfile = InflowProfile(noise);

    NoSlip_T noSlip(blocks, pdfFieldId);
    DynamicUBB_T dynamicUBB(blocks, pdfFieldId, inflowProfile);
    Outflow_T outflow(blocks, pdfFieldId, real_c(1.0));

    noSlip.fillFromFlagField<FlagField_T>(blocks, flagFieldId, NoSlipFlagUID, FluidFlagUID);
    dynamicUBB.fillFromFlagField<FlagField_T>(blocks, flagFieldId, InflowUID, FluidFlagUID);
    outflow.fillFromFlagField<FlagField_T>(blocks, flagFieldId, OutflowUID, FluidFlagUID);

    nlohmann::json metaData_json;
    std::ifstream metaDataFile{"metaData.json"};
    metaDataFile >> metaData_json;
    auto timeSteps = metaData_json["timesteps"];

    WALBERLA_ROOT_SECTION()
    {
        WALBERLA_LOG_DEVEL_VAR(timeSteps);
    }
    
    SweepTimeloop timeloop(blocks->getBlockStorage(), timeSteps);

    blockforest::communication::UniformBufferedScheme<Stencil_T> communication(blocks);
    communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));

    /// REMARK: never use write "timeloop.add() << Sweep(A) << Sweep(B);" (see merge request !486)
    timeloop.add() << BeforeFunction(communication, "communication") << Sweep(noSlip);
    //timeloop.add() << Sweep(simpleUBB);
    timeloop.add() << Sweep(dynamicUBB);
    timeloop.add() << Sweep(outflow);
    timeloop.add() << Sweep(cumulantMRTSweep);

    // Time logger
    timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), parameters.remainingTimeLoggerFrequency_),
                                 "remaining time logger");
    timeloop.addFuncAfterTimeStep(makeSharedFunctor(timeTracker), "time tracking");

    // 
    auto updateBdyValues = [&](){
        dynamicUBB.fillFromFlagField<FlagField_T>(blocks, flagFieldId, InflowUID, FluidFlagUID);
    };

    auto updateNoiseValues = [&](){
        noise->clear();
        std::string fileName = "data/inflow_profile_" + std::to_string(static_cast<int>(timeTracker->getTime())) + ".json";
        std::ifstream file_streamNoise(fileName);
        nlohmann::json inflow_profile;
        file_streamNoise >> inflow_profile;
        for (std::size_t idx = 0; idx != inflow_profile.size(); ++idx)
        {
            noise->insert({{static_cast<int>(inflow_profile[idx][0]),
                            static_cast<int>(inflow_profile[idx][1]),
                            static_cast<int>(inflow_profile[idx][2])}, 
                            Vector3<real_t>{inflow_profile[idx][3], 
                                            inflow_profile[idx][4], 
                                            inflow_profile[idx][5]}});
        }
    };

    timeloop.addFuncAfterTimeStep(updateNoiseValues, "update noise values");
    timeloop.addFuncAfterTimeStep(updateBdyValues, "update bdy values");

    // set velocity in boundary cells to zero
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

    if (parameters.vtkWriteFrequency_ > 0)
    {
        std::string const path = "vtk_out/";
        /// REMARK: force_pvtu must be set to true; otherwise the "pvti" format is used and paraview can only display
        /// pvti files if the domain is partitioned in blocks with a block for EVERY part of the domain; this is likely
        /// to be not the case when using smaller dx
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "cumulant_mrt_velocity_field", parameters.vtkWriteFrequency_, 0,
                                                        true, path, "simulation_step", false, true, true, false, 0);

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
    //std::cout << "Test succesful\n";
}
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

#include "CumulantMRTPackInfo.h"
#include "CumulantMRTSweep.h"
#include "InitialPDFsSetter.h"
//#include "BoundaryHandling.h"
#include "aliases.h"
#include "Parameters.h"

#include <iostream>

using namespace walberla;

FlagUID const FluidFlagUID("Fluid Flag");
FlagUID const NoSlipFlagUID("NoSlip Flag");
FlagUID const OutflowUID("Outflow Flag");
FlagUID const InflowUID("Inflow Flag");

constexpr Mesh::Color noSlipColor{255, 255, 255};

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

Parameters constructParameters(std::shared_ptr<Config> config)
{
    auto parameters = config->getOneBlock("Parameters");

    real_t omega = parameters.getParameter< real_t >("omega", real_c(1.8));
    Vector3< real_t > initialVelocity =
            parameters.getParameter< Vector3< real_t > >("initialVelocity", Vector3< real_t >());
    uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(1000));

    double remainingTimeLoggerFrequency =
        parameters.getParameter< double >("remainingTimeLoggerFrequency", 3.0); // in seconds
    const uint_t VTKwriteFrequency = parameters.getParameter< uint_t >("VTKwriteFrequency", 1000);

    // read domain parameters
    auto domainParameters = config->getOneBlock("DomainSetup");
    std::string meshFile = domainParameters.getParameter< std::string >("meshFile");

    real_t const dx = domainParameters.getParameter< real_t >("dx", real_t(1));
    Vector3<uint_t> const cellsPerBlock = domainParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");

    auto stabilityCheckerParam = config->getOneBlock("StabilityChecker");
    uint_t stabilityChecker = stabilityCheckerParam.getParameter<uint_t>("StabilityChecker", 1000);
    return Parameters(dx, initialVelocity, timesteps, VTKwriteFrequency, stabilityChecker, cellsPerBlock,
                      remainingTimeLoggerFrequency, meshFile, omega);
}

std::shared_ptr<Mesh> readMesh(const std::string& fileName)
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
  mesh::ComplexGeometryStructuredBlockforestCreator bfc(aabb, Vector3<real_t>(params.dx_),
                                                        mesh::makeExcludeMeshInterior(distanceOctree, params.dx_)
                                                       );
    bfc.setPeriodicity(Vector3<bool>(false));

    return bfc.createStructuredBlockForest(params.cellsPerBlock_);
}


int main(int argc, char* argv[])
{
    walberla::Environment walberlaEnv(argc, argv);
    mpi::MPIManager::instance()->useWorldComm();
    auto config = walberlaEnv.config();
    auto parameters = constructParameters(config);
    auto mesh = readMesh(parameters.meshFileName_);
    auto distanceOctree = buildDistanceOctree(mesh);
    auto aabb = computeAABB(*mesh);
    auto blocks = createBlockForest(parameters, aabb, distanceOctree);

    mesh::BoundarySetup boundarySetup(blocks, makeMeshDistanceFunction(distanceOctree), numGhostLayers);

    BlockDataID velocityFieldId = field::addToStorage< VectorField_T >(blocks, "velocity", real_c(0.0), field::fzyx);
    BlockDataID flagFieldId     = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");
    BlockDataID pdfFieldId      = field::addToStorage< PdfField_T >(blocks, "pdf field", real_c(0.0), field::fzyx);

    pystencils::CumulantMRTSweep cumulantMRTSweep(pdfFieldId, velocityFieldId, parameters.omega_);

    // Register which cells have a NoSlip boundary
    NoSlip_T noSlip(blocks, pdfFieldId);
    noSlip.fillFromFlagField<FlagField_T>(blocks, flagFieldId, FlagUID("NoSlip"), FluidFlagUID);
    mesh::ColorToBoundaryMapper<mesh::TriangleMesh> colorToBoundaryMapperNoSlip(
        mesh::BoundaryInfo(mesh::BoundaryInfo(BoundaryUID("NoSlip")))
    );
    colorToBoundaryMapperNoSlip.set(
        mesh::TriangleMesh::Color(255, 255, 255),
        mesh::BoundaryInfo(BoundaryUID("NoSlip"))
    );

    // Register which cells are the inflow
    SimpleUBB_T simpleUBB(blocks, pdfFieldId, 1.0, 1.0, 1.0);
    simpleUBB.fillFromFlagField<FlagField_T>(blocks, flagFieldId, FlagUID("SimpleUBB"), FluidFlagUID);
    mesh::ColorToBoundaryMapper<mesh::TriangleMesh> colorToBoundaryMapperSimpleUBB(
        mesh::BoundaryInfo(mesh::BoundaryInfo(BoundaryUID("SimpleUBB")))
    );
    colorToBoundaryMapperSimpleUBB.set(
        mesh::TriangleMesh::Color(255, 0, 12),
        mesh::BoundaryInfo(BoundaryUID("SimpleUBB"))
    );
    // Register which cells are the outflow
    Outflow_T outflow(blocks, pdfFieldId);
    outflow.fillFromFlagField<FlagField_T>(blocks, flagFieldId, FlagUID("Outflow"), FluidFlagUID);
    mesh::ColorToBoundaryMapper<mesh::TriangleMesh> colorToBoundaryMapperOutflow(
        mesh::BoundaryInfo(mesh::BoundaryInfo(BoundaryUID("Outflow")))
    );
    colorToBoundaryMapperOutflow.set(
        mesh::TriangleMesh::Color(0, 42, 255),
        mesh::BoundaryInfo(BoundaryUID("Outflow"))
    );

    SweepTimeloop timeloop(blocks->getBlockStorage(), parameters.timeSteps_);

    blockforest::communication::UniformBufferedScheme<Stencil_T> communication(blocks);
    communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));

    timeloop.add() << BeforeFunction(communication, "communication") << Sweep(noSlip) << Sweep(simpleUBB) << Sweep(outflow);
    timeloop.add() << Sweep(cumulantMRTSweep);

    // Time logger
    timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), parameters.remainingTimeLoggerFrequency_),
                                 "remaining time logger");
    timeloop.run();
    std::cout << "Test succesful\n";
}
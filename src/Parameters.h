#ifndef NLO_MT_PARAMETERS_H_
#define NLO_MT_PARAMETERS_H_

#include "aliases.h"

using walberla::real_t;
using walberla::uint_t;
using walberla::Vector3;

struct Parameters
{
    real_t dx_;
    real_t omega_;
    Vector3<real_t> initialVelocity_;
    Vector3<uint_t> cellsPerBlock_;

    uint_t timeSteps_;
    uint_t vtkWriteFrequency_;
    uint_t stabilityCheckFrequency_;
    uint_t numConstNoises_;
    double remainingTimeLoggerFrequency_;

    real_t numHeartBeats_;

    std::string meshFileName_;
    bool generateInflowProfile_;

    Parameters(real_t dx, uint_t timeSteps,
               uint_t vtkWriteFrequency, uint_t stabilityCheckFrequency,
               uint_t numConstNoises,
               Vector3<uint_t> cellsPerBlock, double remainingTimeLoggerFrequency,
               real_t numHeartBeats, std::string const& meshFileName, real_t omega,
               bool generateInflowProfile)
    :
        dx_(dx),
        omega_(omega),
        cellsPerBlock_(cellsPerBlock),
        timeSteps_(timeSteps),
        vtkWriteFrequency_(vtkWriteFrequency),
        stabilityCheckFrequency_(stabilityCheckFrequency),
        numConstNoises_(numConstNoises),
        remainingTimeLoggerFrequency_(remainingTimeLoggerFrequency),
        numHeartBeats_(numHeartBeats),
        meshFileName_(meshFileName),
        generateInflowProfile_(generateInflowProfile)
    {}
};

#endif
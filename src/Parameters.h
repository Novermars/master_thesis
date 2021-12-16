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
    double remainingTimeLoggerFrequency_;

    real_t numHeartBeats_;

    std::string meshFileName_;

    Parameters(real_t dx, uint_t timeSteps,
               uint_t vtkWriteFrequency, uint_t stabilityCheckFrequency,
               Vector3<uint_t> cellsPerBlock, double remainingTimeLoggerFrequency,
               real_t numHeartBeats, std::string const& meshFileName, real_t omega)
    :
        dx_(dx),
        omega_(omega),
        cellsPerBlock_(cellsPerBlock),
        timeSteps_(timeSteps),
        vtkWriteFrequency_(vtkWriteFrequency),
        stabilityCheckFrequency_(stabilityCheckFrequency),
        remainingTimeLoggerFrequency_(remainingTimeLoggerFrequency),
        numHeartBeats_(numHeartBeats),
        meshFileName_(meshFileName)
    {}
};

#endif
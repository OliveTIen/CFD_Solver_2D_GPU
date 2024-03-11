#include "EvolveGPU.h"
#include "../gpu/../SignDefine.h"
#include "../space/FluxGPU.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"

void GPU::Time::EvolveDevice(REAL dt, int flag, ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, std::map<int, int>& edge_periodic_pair) {
    if (flag == _EVO_explicit) {
        EvolveExplicitDevice(dt, element_device, elementField_device, edge_device, edge_periodic_pair);
    }
    else {
        LogWriter::writeLogAndCout("Error: invalid evolve method.\n", LogWriter::Error, LogWriter::Error);
        exit(1);
    }
}

void GPU::Time::EvolveExplicitDevice(REAL dt, ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, std::map<int, int>& edge_periodic_pair) {
    // 数值通量
    U2NITS::Space::Flux::calculateFluxHost(element_device, elementField_device, edge_device, edge_periodic_pair);
    // 时间推进
    for (int ie = 0; ie < element_device.num_element; ie++) {
        REAL omega = element_device.volume[ie];
        for (int j = 0; j < 4; j++) {
            elementField_device.U[j][ie] -= dt / omega * elementField_device.Flux[j][ie];
        }
    }
}

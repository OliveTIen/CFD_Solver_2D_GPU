#include "Evolve.h"
#include "../global/GlobalPara.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "../math/Math.h"

void U2NITS::Time::EvolveHost(REAL dt, int flag, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    if (flag == _EVO_explicit) {
        EvolveExplicitHost(dt, element_host, elementField_host, edge_host, edge_periodic_pair);
    }
    else {
        LogWriter::writeLogAndCout("Error: invalid evolve method.\n", LogWriter::Error, LogWriter::Error);
        exit(1);
    }
}

void U2NITS::Time::EvolveExplicitHost(REAL dt, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    // 数值通量
    U2NITS::Space::Flux::calculateFluxHost(element_host, elementField_host, edge_host, edge_periodic_pair);
    // 时间推进
    for (int ie = 0; ie < element_host.num_element; ie++) {
        REAL omega = element_host.volume[ie];
        for (int j = 0; j < 4; j++) {
            elementField_host.U[j][ie] -= dt / omega * elementField_host.Flux[j][ie];
        }

        //Math::restrictU2d(
        //    elementField_host.U[0][ie],
        //    elementField_host.U[1][ie],
        //    elementField_host.U[2][ie],
        //    elementField_host.U[3][ie]
        //);
    }
}

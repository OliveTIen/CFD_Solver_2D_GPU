#include "Evolve.h"
#include "../../SignDefine.h"
#include "../space/CalculateFlux2.cuh"

void GPU::Time::EvolveHost(REAL dt, int flag, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    if (flag == _EVO_explicit) {
        EvolveExplicitHost(dt, element_host,elementField_host,edge_host, edge_periodic_pair);
    }
    else if (flag == _EVO_rk3) {
        //TODO: RK3
    }
    else {
        EvolveExplicitHost(dt, element_host, elementField_host, edge_host, edge_periodic_pair);
    }
}

void GPU::Time::EvolveExplicitHost(REAL dt, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    // 数值通量
    GPU::calculateFluxHost(element_host, elementField_host, edge_host, edge_periodic_pair);
    // 时间推进
    for (int ie = 0; ie < element_host.num_element; ie++) {
        REAL omega = element_host.volume[ie];
        for (int j = 0; j < 4; j++) {
            elementField_host.U[j][ie] -= dt / omega * elementField_host.Flux[j][ie];
        }
    }
}

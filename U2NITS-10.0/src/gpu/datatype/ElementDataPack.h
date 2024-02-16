#ifndef __ELEMENT_DATA_PACK_H__
#define __ELEMENT_DATA_PACK_H__

#include "Define.h"

namespace GPU {
    // µ±Ç°
    class Element {
    public:
        bool* isNull;

        REAL* x;
        REAL* y;

        REAL* U1;
        REAL* U2;
        REAL* U3;
        REAL* U4;

        REAL* Ux1;
        REAL* Ux2;
        REAL* Ux3;
        REAL* Ux4;
        REAL* Uy1;
        REAL* Uy2;
        REAL* Uy3;
        REAL* Uy4;

    public:
        void alloc(int num_element) {
            this->isNull = new bool[num_element] {};
            this->x = new REAL[num_element]{};
            this->y = new REAL[num_element]{};

            this->U1 = new REAL[num_element]{};
            this->U2 = new REAL[num_element]{};
            this->U3 = new REAL[num_element]{};
            this->U4 = new REAL[num_element]{};

            this->Ux1 = new REAL[num_element]{};
            this->Ux2 = new REAL[num_element]{};
            this->Ux3 = new REAL[num_element]{};
            this->Ux4 = new REAL[num_element]{};
            this->Uy1 = new REAL[num_element]{};
            this->Uy2 = new REAL[num_element]{};
            this->Uy3 = new REAL[num_element]{};
            this->Uy4 = new REAL[num_element]{};
        }
        void free() {
            delete[] isNull;
            delete[] x;
            delete[] y;

            delete[] U1;
            delete[] U2;
            delete[] U3;
            delete[] U4;

            delete[] Ux1;
            delete[] Ux2;
            delete[] Ux3;
            delete[] Ux4;
            delete[] Uy1;
            delete[] Uy2;
            delete[] Uy3;
            delete[] Uy4;
        }
        void cuda_alloc(int num_element) {
            cudaMalloc((void**)&this->isNull, sizeof(bool) * num_element);
            cudaMalloc((void**)&this->x, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->y, sizeof(REAL) * num_element);

            cudaMalloc((void**)&this->U1, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->U2, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->U3, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->U4, sizeof(REAL) * num_element);

            cudaMalloc((void**)&this->Ux1, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->Ux2, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->Ux3, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->Ux4, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->Uy1, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->Uy2, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->Uy3, sizeof(REAL) * num_element);
            cudaMalloc((void**)&this->Uy4, sizeof(REAL) * num_element);

            if (isNull == nullptr || x == nullptr || y == nullptr
                || U1 == nullptr || U2 == nullptr || U3 == nullptr || U4 == nullptr
                || Ux1 == nullptr || Ux2 == nullptr || Ux3 == nullptr || Ux4 == nullptr
                || Uy1 == nullptr || Uy2 == nullptr || Uy3 == nullptr || Uy4 == nullptr) {
                throw "cuda malloc error";
            }
        }
        void cuda_free() {
            cudaFree(this->isNull);
            cudaFree(this->x);
            cudaFree(this->y);

            cudaFree(this->U1);
            cudaFree(this->U2);
            cudaFree(this->U3);
            cudaFree(this->U4);

            cudaFree(this->Ux1);
            cudaFree(this->Ux2);
            cudaFree(this->Ux3);
            cudaFree(this->Ux4);
            cudaFree(this->Uy1);
            cudaFree(this->Uy2);
            cudaFree(this->Uy3);
            cudaFree(this->Uy4);
        }
        void cuda_copy_to_device(Element* element_device, int num_element) {
            cudaMemcpy(element_device->isNull, this->isNull, sizeof(bool) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->x, this->x, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->y, this->y, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);

            cudaMemcpy(element_device->U1, this->U1, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->U2, this->U2, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->U3, this->U3, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->U4, this->U4, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);

            cudaMemcpy(element_device->Ux1, this->Ux1, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->Ux2, this->Ux2, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->Ux3, this->Ux3, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->Ux4, this->Ux4, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->Uy1, this->Uy1, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->Uy2, this->Uy2, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->Uy3, this->Uy3, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);
            cudaMemcpy(element_device->Uy4, this->Uy4, sizeof(REAL) * num_element, cudaMemcpyHostToDevice);

        }
        void cuda_copy_to_host(Element* element_host, int num_element) {
            cudaMemcpy(element_host->isNull, this->isNull, sizeof(bool) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->x, this->x, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->y, this->y, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);

            cudaMemcpy(element_host->U1, this->U1, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->U2, this->U2, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->U3, this->U3, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->U4, this->U4, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);

            cudaMemcpy(element_host->Ux1, this->Ux1, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->Ux2, this->Ux2, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->Ux3, this->Ux3, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->Ux4, this->Ux4, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->Uy1, this->Uy1, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->Uy2, this->Uy2, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->Uy3, this->Uy3, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
            cudaMemcpy(element_host->Uy4, this->Uy4, sizeof(REAL) * num_element, cudaMemcpyDeviceToHost);
        }

    };

    class ElementDataPack {
    public:
        GPU::Element self;
        GPU::Element neighbors[3];

        int num_element = 0;

    public:
        void alloc(int num_element) {
            this->num_element = num_element;
            self.alloc(num_element);
            neighbors[0].alloc(num_element);
            neighbors[1].alloc(num_element);
            neighbors[2].alloc(num_element);
        }
        void free() {
            self.free();
            neighbors[0].free();
            neighbors[1].free();
            neighbors[2].free();
            num_element = 0;
        }
        void cuda_alloc(int num_element) {
            this->num_element = num_element;
            self.cuda_alloc(num_element);
            neighbors[0].cuda_alloc(num_element);
            neighbors[1].cuda_alloc(num_element);
            neighbors[2].cuda_alloc(num_element);
        }
		void cuda_free() {
            self.cuda_free();
            neighbors[0].cuda_free();
            neighbors[1].cuda_free();
            neighbors[2].cuda_free();
            num_element = 0;
        }
        void cuda_copy_to_device(ElementDataPack* device) {
            self.cuda_copy_to_device(&(device->self), num_element);
			neighbors[0].cuda_copy_to_device(&(device->neighbors[0]), num_element);
			neighbors[1].cuda_copy_to_device(&(device->neighbors[1]), num_element);
			neighbors[2].cuda_copy_to_device(&(device->neighbors[2]), num_element);
        }
		void cuda_copy_to_host(ElementDataPack* host) {
            self.cuda_copy_to_host(&(host->self), num_element);
			neighbors[0].cuda_copy_to_host(&(host->neighbors[0]), num_element);
			neighbors[1].cuda_copy_to_host(&(host->neighbors[1]), num_element);
			neighbors[2].cuda_copy_to_host(&(host->neighbors[2]), num_element);
        }   
    };
    /*
    inline void cuda_memcpy(ElementDataPack* dist, const ElementDataPack* src, cudaMemcpyKind kind) {

    }
    */
}


#endif
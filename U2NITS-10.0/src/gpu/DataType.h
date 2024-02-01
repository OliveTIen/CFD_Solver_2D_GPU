#ifndef GPU_DATATYPE_H
#define GPU_DATATYPE_H

#include <cuda_runtime.h>

// define data type
#define REAL double

namespace GPU {
    // 未使用
    class Coordinate {
    public:
        REAL* x;
        REAL* y;
    public:
		void alloc(int num_element) {
            this->x = new REAL[num_element]{};
            this->y = new REAL[num_element]{};
		}
		void free() {
			delete[] x;
			delete[] y;
		}
    };
    // 未使用 守恒量
    class Value {
    public:
        REAL* U1;
        REAL* U2;
        REAL* U3;
        REAL* U4;
    public:
        void alloc(int num_element) {
            this->U1 = new REAL[num_element]{};
            this->U2 = new REAL[num_element]{};
            this->U3 = new REAL[num_element]{};
            this->U4 = new REAL[num_element]{};
        }
        void free() {
            delete[] U1;
            delete[] U2;
            delete[] U3;
            delete[] U4;
        }
    };
    // 未使用
    class ValueWithGradient: public Value {
    public:
        bool* isNull;
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
            Value::alloc(num_element);
            this->isNull = new bool[num_element]{};
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
            Value::free();
			delete[] isNull;
            delete[] Ux1;
            delete[] Ux2;
            delete[] Ux3;
            delete[] Ux4;
            delete[] Uy1;
            delete[] Uy2;
            delete[] Uy3;
            delete[] Uy4;
        }
    };

    // 当前
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
            this->isNull = new bool[num_element]{};
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
    };

    //class ElementPlus : public Element {
    //public:
    //    REAL* flux;
    //};

    class Edge {
    public:
        int* element_index1;
        int* element_index2;
        int* node_index1;
        int* node_index2;
        int num_edge = 0;
    };

}

#endif // !GPU_DATATYPE_H

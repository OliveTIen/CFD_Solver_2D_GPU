#ifndef __ELEMENT_DATA_PACK_H__
#define __ELEMENT_DATA_PACK_H__

#include "DataType.h"

namespace GPU {
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
}


#endif
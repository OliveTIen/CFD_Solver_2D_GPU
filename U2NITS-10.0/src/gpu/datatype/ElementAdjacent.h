#ifndef __ELEMENTADJACENT_H__
#define __ELEMENTADJACENT_H__

// 类中函数默认是inline，因此可以不放进cpp文件
// 全部写在头文件中即可

namespace GPU {
    class ElementAdjacent {

    private:
        class PrivateNeighbor {
        public:
            int* index;
            void alloc(int num_element) {
                index = new int[num_element];
                for (int i = 0; i < num_element; i++) {
                    index[i] = -1;
                }
            }
            void free() {
                delete[] index;
            }
        };

    public:
        int m_num_element = 0;
        PrivateNeighbor* neighbors;
        void alloc(int num_neighbor, int num_element) {
            m_num_element = num_element;
            neighbors = new PrivateNeighbor[num_neighbor];
            for (int i = 0; i < num_neighbor; i++) {
                neighbors[i].alloc(num_element);
            }
        }
        void free() {
            for (int i = 0; i < m_num_element; i++) {
                neighbors[i].free();
            }
            m_num_element = 0;
        }
    };

}


#endif
#ifndef GPU_SOLVER_H
#define GPU_SOLVER_H

#include "DataType.h"
#include "ElementDataPack.h"

class GPUSolver {
public:
	GPU::ElementDataPack element_host;
	GPU::ElementDataPack element_device;

public:
	void initialze();

	void iteration();

    void free(){this->element_host.free();}
	//void run();

};



#endif // GPU_SOLVER_H

/*
����1����Ԫ�����������ݣ�
nodeCoordinates[3][3]
values[4]
neighborValues[3][3]

��ȡ�ļ��󣬽�FVM_2D�����ݽṹת��ΪElement�ṹ��
    1. ��ʼ��selfxy, node1, node2, node3 �������������⣬self����CPU�ϼ���ģ��Ժ��ת�Ƶ�GPU�ϣ�
    2. ��ʼ��selfValue.U1~U4, neighbor1.U1~U4, neighbor2.U1~U4, neighbor3.U1~U4
    3. �ݶȳ�ʼ��Ϊ0

�����ݶȣ���ÿ��ѭ������Ҫ���㣩
    ��Element����GPU
    ���㲢���µ�Ԫ�ݶ�(self.Ux1~Uy4)��Ȼ��ͬ��(barrier)���ȴ��ݶȼ������
    ������CPU
    ��CPU�н���Ԫ�ݶȴ��ݸ��ھӵ�Ԫ������neighbor1.Ux1~Uy4, neighbor2.Ux1~Uy4, neighbor3.Ux1~Uy4

����߽�ֵ��
���ݵ�Ԫ�ڵ����ꡢ�ھӵ�Ԫ�ݶȡ��ھӵ�Ԫֵ����߽�ֵ��ͬ�����ȴ��߽�ֵ�������
���ݱ߽�ֵ�����������������߽���ֵͨ��
������ֵͨ�������㲢���µ�Ԫֵself.U1~U4
����Ԫֵ���ݸ��ھӵ�Ԫ������neighbor1.U1~U4, neighbor2.U1~U4, neighbor3.U1~U4

ע�⣺
1. Element�ĳ�Ա��ָ�룬����ֱ�Ӵ����õ���Щ�ʹ���Щ

*/
#ifndef GMESH_MESH_READER
#define GMESH_MESH_READER
#include <string>
// �������ڰ�װ��ȡ����Ĳ���������FVM_2D�����ӷ�׳̶�
class GmeshMeshReader {

public:
	// ��ȡinp�������ݣ���ʼ���������
	static int readMeshFile(std::string suffix);
};
#endif // !GMESH_MESH_READER

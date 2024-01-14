#ifndef GMESH_MESH_READER
#define GMESH_MESH_READER
#include <string>
// 该类用于包装读取网格的操作，减少FVM_2D代码的臃肿程度
class GmeshMeshReader {

public:
	// 读取inp网格数据，初始化网格变量
	static int readMeshFile(std::string suffix);
};
#endif // !GMESH_MESH_READER

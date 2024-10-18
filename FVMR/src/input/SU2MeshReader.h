#ifndef SU2_MESH_READER
#define SU2_MESH_READER
#include <string>
#include <vector>
class SU2MeshReader {
public:

	// 状态。读取网格时，读取到某一位置，会有特定状态
	enum State {
		state_START,
		state_NDIME,
		state_NPOIN,
		state_NPOIN_data,
		state_NELEM,
		state_NELEM_data,
		state_NMARK,
		state_NMARK_data,
	};

	// 存储两个节点
	class SimpleEdge {
	public:
		int nodeIDs[2]{-1,-1};
		SimpleEdge() {};
		SimpleEdge(int nodeID_1, int nodeID_2) { nodeIDs[0] = nodeID_1; nodeIDs[1] = nodeID_2; }
	};
	// 存储名称，以及边的节点
	class SimpleBoundary {
	public:
		std::string name;
		std::vector<SimpleEdge> edges;

		SimpleBoundary() {}
		SimpleBoundary(std::string name_) { name = name_; }
		SimpleBoundary(const std::string& name, const std::vector<SimpleEdge>& edges)
			: name(name), edges(edges) {
		}
	};

public:
	// readFile的拆分版，拆成两个子函数
	static void read_mesh_and_process(std::string filePath, bool convertRectToTriangle);
	// 读取网格文件
	static void readMesh(std::string filePath, bool convertRectToTriangle, int& maxNodeID, int& maxElementID, std::vector<SimpleBoundary>& tmp_boudaries);
	// 缩放网格坐标。应在process前调用
	static void scaleMesh(double factor);
	static void process(int maxNodeID, int maxElementID, std::vector<SimpleBoundary>&tmp_boudaries);
	// 
};

#endif // !SU2_MESH_READER

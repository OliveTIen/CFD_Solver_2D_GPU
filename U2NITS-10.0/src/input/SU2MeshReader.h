#ifndef SU2_MESH_READER
#define SU2_MESH_READER
#include <string>
#include <vector>
class SU2MeshReader {
public:
	enum ErrorCode {
		E_NO_ERROR,
		ERROR_READ_FILE,
		ERROR_PARSE
	};
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

	
	class SimpleEdge {
	public:
		int nodeIDs[2]{-1,-1};
		SimpleEdge() {};
		SimpleEdge(int nodeID_1, int nodeID_2) { nodeIDs[0] = nodeID_1; nodeIDs[1] = nodeID_2; }
	};
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
	static int readFile_2(std::string filePath, bool convertRectToTriangle);
	static int readMesh(std::string filePath, bool convertRectToTriangle, int& maxNodeID, int& maxElementID, std::vector<SimpleBoundary>& tmp_boudaries);
	static void process(int maxNodeID, int maxElementID, std::vector<SimpleBoundary>&tmp_boudaries);
	// 
};

#endif // !SU2_MESH_READER

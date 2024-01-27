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
	static int readFile(std::string filePath, bool convertRectToTriangle);
};

#endif // !SU2_MESH_READER

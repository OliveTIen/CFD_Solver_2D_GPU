#ifndef SU2_MESH_READER
#define SU2_MESH_READER
#include <string>
#include <vector>

/**
 * \class SU2MeshReader
 * \brief
 * Read su2 mesh and stores data into FVM_2D
 * \details public member function @ref SU2MeshReader::read_mesh_and_process 
 * is the interface.<br>
 * TODO:<br>
 * The problem is it's coupling with legacy FVM_2D, and after reading mesh the
 * data stored in FVM_2D is unused but not freed, leading to memory waste.
 */
class SU2MeshReader {
public:

	/**
	 * Status when reading mesh file.
	 * It's like a state machine.
	*/
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

	/**
	 * Stores two nodes
	*/
	class SimpleEdge {
	public:
		int nodeIDs[2]{-1,-1};
		SimpleEdge() {};
		SimpleEdge(int nodeID_1, int nodeID_2) { nodeIDs[0] = nodeID_1; nodeIDs[1] = nodeID_2; }
	};

	/**
	 * Stores name and nodes of edges
	*/
	class SimpleBoundary {
	public:
		std::string name;
		std::vector<SimpleEdge> edges;

		SimpleBoundary() {}
		SimpleBoundary(std::string name_) { name = name_; }
		SimpleBoundary(const std::string& name, const std::vector<SimpleEdge>& edges)
			: name(name), edges(edges) {}
	};

public:
	/**
	 * Read mesh, scale mesh, and process
	 * @param filePath - file path.
	 * @param convertRectToTriangle(bool) -  convert to triangle.
	 * @return none
	*/
	static void read_mesh_and_process(std::string filePath, bool convertRectToTriangle);

private:
	/**
	 * Read su2 mesh file, including node indices, node coordinates, element indices, and element-node relations
	 * @param filePath - file path.
	 * @param tmp_boudaries - out param, storing boundaries.
	 * @return none
	*/
	static void readMesh(std::string filePath, bool convertRectToTriangle, int& maxNodeID, int& maxElementID, std::vector<SimpleBoundary>& tmp_boudaries);
	/**
	 * Scale mesh node coordinates. Should called before process
	 * @param factor - scale factor.
	*/
	static void scaleMesh(double factor);
	/**
	 * Initialize Edges, pNodeTable and so on, and initialize boundaries.
	 * @param maxNodeID - max node ID.
	 * @param maxElementID - max element ID.
	*/
	static void process(int maxNodeID, int maxElementID, std::vector<SimpleBoundary>&tmp_boudaries);
};

#endif // !SU2_MESH_READER

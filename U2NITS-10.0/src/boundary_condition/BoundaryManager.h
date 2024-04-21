#ifndef BOUNDARY_MANAGER_H
#define BOUNDARY_MANAGER_H
#include <string>
namespace U2NITS {
	class BoundaryManager {
	public:
		static int boundaryNameToType(std::string name);
	};
}
#endif // !BOUNDARY_MANAGER_H

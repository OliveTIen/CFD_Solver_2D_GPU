#ifndef BOUNDARY_MANAGER_H
#define BOUNDARY_MANAGER_H
#include <string>
namespace U2NITS {
	/*
	ÐÂBoundaryManager
	*/
	class BoundaryManager {
	private:
		static BoundaryManager* p_instance;
		BoundaryManager() {};

	public:
		static BoundaryManager* getInstance();
		static int boundaryNameToType(std::string name);

	private:
		
	};
}
#endif // !BOUNDARY_MANAGER_H

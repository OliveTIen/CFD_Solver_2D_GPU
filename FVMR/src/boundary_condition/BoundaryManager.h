#ifndef BOUNDARY_MANAGER_H
#define BOUNDARY_MANAGER_H
#include <string>
namespace U2NITS {
	/*
	* \class BoundaryManager
	* \brief
	* Convert boundary name to type
	*/
	class BoundaryManager {
	private:
		static BoundaryManager* p_instance;
		BoundaryManager() {};

	public:
		static BoundaryManager* getInstance();
		/*
		* Convert boundary name to type
		* @param name - boundary name
		*/
		static int boundaryNameToType(std::string name);

	private:
		
	};
}
#endif // !BOUNDARY_MANAGER_H

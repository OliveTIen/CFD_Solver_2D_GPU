#pragma once

typedef double su2double;
typedef int SU2_Comm;
#include <iostream>

/*!
* \class CDriverBase
* \ingroup Drivers
* \brief Base class for all drivers.
* \author <NAME>
* \date 2024-03-09
*/

class CDriverBase {
protected:
	int rank;
	int size;
	char* config_file_name;
	unsigned short nZone;

	//CConfig* driver_config = nullptr; /*!< \brief Definition of the driver configuration. */
	//COutput* driver_output = nullptr; /*!< \brief Definition of the driver output. */
	//CConfig** config_container;

	//CConfig* main_config = nullptr;     /*!< \brief Reference to base (i.e. ZONE 0) configuration (used in driver API). */
	//CGeometry* main_geometry = nullptr; /*!< \brief Reference to base (i.e. ZONE, INST, MESH 0) geometry (used in driver API). */

public:
	CDriverBase(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator);
	virtual ~CDriverBase(void);
	virtual void Run() {};
	virtual void Finalize() {};

protected:
	void InitializeContainers();
};
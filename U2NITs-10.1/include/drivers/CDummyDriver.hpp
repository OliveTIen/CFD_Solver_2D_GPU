
#pragma once
#include "CDriver.hpp"

/*!
 * \brief CDummyDriver class that constructs the driver without running a solver.
 * \ingroup Drivers
 */
class CDummyDriver : public CDriver {

public:

    /*!
     * \brief Constructor of the class.
     * \param[in] confFile - Configuration file name.
     * \param[in] val_nZone - Total number of zones.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CDummyDriver(char* confFile,
        unsigned short val_nZone
    );

    /*!
     * \brief Destructor of the class.
     */
    ~CDummyDriver() override {}

    /*!
     * \brief Does nothing except printing the information that no solver is running.
     */
    void StartSolver() override;

};

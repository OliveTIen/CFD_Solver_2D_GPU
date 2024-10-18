/*!
 * \file main.cpp
 * \brief main file
 * \details
 * This is the entry point that calls @ref U2NITS::CDriver::start()
 * \author Olive Tien
 */

#include "drivers/CDriver.h"

int main() {
	U2NITS::CDriver::getInstance()->start();
	return 0;
}

// Doxygen标记参考: https://blog.csdn.net/qq_36631379/article/details/121980455
/**
 * \mainpage FVMR Documentation
 * #### Description
 * FVMR (Finite Volume Method Realtime) is a 2D fluid solver written in C++/CUDA,
 * based on the finite volume method and unstructured grids, capable of real-time
 * visualization of solution results. This repository contains the project source
 * code and example grids. The project runs on the Windows platform. <br>
 * #### Quick Tour
 * To quickly familiarize yourself with the code, 
 * please start with the entry point @ref main.cpp .
 * 
 * 
 */
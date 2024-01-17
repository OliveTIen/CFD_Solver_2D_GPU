#ifndef HEAD_H
#define HEAD_H

#include <omp.h>

#include <conio.h> 
#include <iostream>
#include <io.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#elif defined __linux__
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <stdlib.h> 
#include "./include/Eigen/Core"
#include "./include/Eigen/Dense"
#include <time.h>
#include <ctime>
#include <thread>
//#include "Math.h"
#include "./include/rapidjson/document.h"
#include "./include/rapidjson/prettyWriter.h"
#include "./include/rapidjson/stringbuffer.h"

#include "SignDefine.h"
#include "globalPara.h"

#ifdef _WIN32

#elif defined __linux__
typedef struct _COORD {
	short X;
	short Y;
} COORD, * PCOORD;
#endif




#endif

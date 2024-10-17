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
#include <time.h>
#include <ctime>
#include <thread>


#include "../global/globalPara.h"
#include "../gpu/datatype/DefineType.h"

#ifdef _WIN32

#elif defined __linux__
typedef struct _COORD {
	short X;
	short Y;
} COORD, * PCOORD;
#endif




#endif

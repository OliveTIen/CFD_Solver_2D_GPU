#include "CExit.h"
#include "../drivers/CDriver.h"
#include <stdio.h>

void CExit::saveAndExit(int _Code) {
	U2NITS::CDriver::saveAndExit(_Code);
}

void CExit::pressAnyKeyToExit() {

	rewind(stdin);// 清空缓冲区，否则会立刻退出
	exit(pressAnyKeyToContinue());
}

int CExit::pressEnterToContinue() {

	int c = 0;
	printf("Press ENTER to continue... ");
	fflush(stdout);
	do c = getchar(); while ((c != '\n') && (c != EOF));
	return c;
}

int CExit::pressAnyKeyToContinue() {

	//// 已弃用。测试失败。不行。只有输入enter才起作用
	//int c = 256;// char最大不超过127
	//printf("Press any key to continue... ");
	//fflush(stdout);
	//do c = getchar(); while ((c == 256));
	//return c;


	printf("Press any key to continue... ");
	while (1) {
		if (_kbhit()) {
			return _getch();
		}
	}

}

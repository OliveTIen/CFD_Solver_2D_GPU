#include "CExit.h"
#include "../drivers/CDriver.h"
#include <stdio.h>

void CExit::saveAndExit(int _Code) {
	U2NITS::CDriver::saveAndExit(_Code);
}

void CExit::pressAnyKeyToExit() {

	rewind(stdin);// ��ջ�����������������˳�
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

	//// �����á�����ʧ�ܡ����С�ֻ������enter��������
	//int c = 256;// char��󲻳���127
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

#ifndef CEXIT_H
#define CEXIT_H

class CExit {
public:
	// ���ڵ���������ʹ�á�������һ�������ļ����˳�
	static void saveAndExit(int _Code);
	// ��������˳���Ŀǰֻ�ܰ�enter�˳�
	static void pressAnyKeyToExit();
	static int pressEnterToContinue();
	static int pressAnyKeyToContinue();
};

#endif // !CEXIT_H

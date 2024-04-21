#ifndef CEXIT_H
#define CEXIT_H

class CExit {
public:
	// 仅在迭代过程中使用。保存上一步续算文件并退出
	static void saveAndExit(int _Code);
	// 按任意键退出。目前只能按enter退出
	static void pressAnyKeyToExit();
	static int pressEnterToContinue();
	static int pressAnyKeyToContinue();
};

#endif // !CEXIT_H

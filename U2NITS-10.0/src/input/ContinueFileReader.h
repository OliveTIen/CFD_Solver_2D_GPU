#ifndef _CONTINUEFILEREADER_H_
#define _CONTINUEFILEREADER_H_

#include <string>
#include <vector>


/*
由于FVM_2D::readContinueFile()函数包含对指针表初始化等操作，很复杂，
目前先不忙完成该类
*/
class ContinueFileReader {
	// 建议将reader writer合并 放进fileio文件夹 查询一下其他软件的目录组织

public:
	static int readContinueFile_1();
	// 添加Ux Uy 目前用不到
	static int readContinueFile_2_unused_addUxUy();

	static int readContinueFile_1_1();

private:
	static void process(int maxNodeID, int maxElementID, std::vector<std::vector<int>>& edges_of_all_sets);

};

#endif // _CONTINUEFILEREADER_H_
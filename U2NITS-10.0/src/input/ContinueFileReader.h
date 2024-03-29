#ifndef _CONTINUEFILEREADER_H_
#define _CONTINUEFILEREADER_H_

#include <string>

/// <summary>
///  由于FVM_2D::readContinueFile()函数包含对指针表初始化等操作，很复杂，
/// 目前先不忙完成该类
/// </summary>
class ContinueFileReader {
	// 建议将reader writer合并 放进fileio文件夹 查询一下其他软件的目录组织

public:
	static int readContinueFile_1();
	// [todo]不经过FVM_2D
	static int readContinueFile_2();

};

#endif // _CONTINUEFILEREADER_H_
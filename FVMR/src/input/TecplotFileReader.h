#ifndef TECPLOT_FILE_READER
#define TECPLOT_FILE_READER
#include <string>

/*
2024-05-27 待完成
在计算Ma0.77算例时，中途出现NaN，没有pause续算文件，且recovery文件被意外覆盖
一种解决办法是从Euler方程算例开始续算
另一种长久之计是编写从Tecplot输出文件续算的代码。Tecplot是基于格点的，读取后需要转换到格心
并且tecplot文件所采用的编号是GPUID而不是ID
*/

class TecplotFileReader {
public:
	static TecplotFileReader* getInstance();

private:
	static TecplotFileReader* p_instance;

	TecplotFileReader() {};
};

#endif
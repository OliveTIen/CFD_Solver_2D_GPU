/*
具有正交属性的类，适合用桥接模式
*/

#ifndef _CSOLVER_H_
#define _CSOLVER_H_
#include "equation/CEquation.h"
#include "platform/CPlatform.h"

namespace U2NITS {
	class CSolver {
	private:
		//CEquation* m_equation;
		//CPlatform* m_platform;
	public:
		virtual void initialize() = 0;
		virtual void iterate() = 0;
		virtual void updateResidual() = 0;
		virtual void finalize() = 0;
	};
}
#endif // !_CSOLVER_H_

/*
步骤：（优化后）
读取配置（控制参数） FileReader.readConfig()
初始化流场：读取续算文件 / 读取网格文件 + 初始化 FileReader.readField()
输出初始信息（程序信息、边界参数等）ConsolePrinter.updateScreen()
输出残差文件头 FileIOManager.writeHistFileHead()
计时器开始
初始化资源（GPU内存）Solver.initialize()
迭代开始
更新部分数据（文件名）
计算（流场） Solver.iterate()
计算（残差）仅需输出时计算 Solver.updateResidual()
更新输出信息（包括记录下一次光标位置、输出进度）ConsolePrinter.updateScreen()
输出文件（流场、残差）FileWriter.write()
判断是否跳出循环（终止）this->checkStop()
输出结束信息（提示词）ConsolePrinter.updateScreen()
释放资源 Solver.finalize()
*/

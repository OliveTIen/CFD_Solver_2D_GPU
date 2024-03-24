﻿---
date:   2024-01-12
author: tgl

---

# 日志
## 代办
3. 尝试课件的算例？？？可是它是非结构网格不可压，不知道是否需要另外写求解器
4. 尝试添加粘性(N-S)处理复杂边界
## 问题
1. 原来的项目中，输入square_per时出现内存问题
## 已解决
### ConsoleWindow问题
以前能够实现光标移动，可是更改代码后发现失效了。经排查，原来是开头输出input.toml时，由于输出内容超过1页，导致屏幕滚动，
然而目前的ConsolePrinter只支持一行滚动的

2024-03-19
在计算梯度的代码中，发现了LeastSquare的内存错误：nValidNeighbor提前+1。
以前直接改成GPU程序，导致该错误没发现，因此引发cuda_error 700异常。
看来在改为GPU程序前，一定要在CPU上测试一下

仍然出现异常。检查文件读取程序，发现element的指标管理混乱

2024-03-16
1.出现bug，debug模式下无法运行，release模式下可以运行；显示栈被破坏；用堆则会陷入死循环
后来发现是程序的问题，原来是U2ruvwp函数没有修改，二维数组用三维的函数造成越界。这说明三维程序改为二维程序时要小心
2.ptxas fatal : Unresolved extern function. 
原因解释：nvcc旧版本不允许__device__函数在别处定义。为了兼容旧版本，-rdc选项是默认关闭的，需要手动打开
(https://blog.csdn.net/weixin_47744790/article/details/134275227)
解决措施：首先，Visual Studio中右击项目属性-CUDA C/C++-Common-Generate Relocatable Device Code设置为true
然后，链接器-输入-附加依赖项-添加cudadevrt.lib
(https://blog.csdn.net/zhangzhe_0305/article/details/79898850)
可是还是报错，于是在FluxGPU.cu中把Roe函数注释掉了
即使我在CUDA Linker中添加了附加依赖项cudadevrt.lib，还是报错。于是撤销该操作。
后来发现是因为我设置的是Release的配置。应该将Debug也设置了。
事实上，只需要开启-rdc选项，无需添加cudadevrt.lib

2024-03-14
性能分析工具
Visual Studio Profiler https://blog.csdn.net/sujunzy666/article/details/19963077


2024-03-08
为什么收敛这么慢？虽然算了接近1万步，但流场变化缓慢
GPU利用率低的原因：
1. 计算量太小，GPU并不能发挥作用
2. GPU内存带宽不足，限制了计算速率
试图模仿SU2进行重构，但是难度很大。这就像我学Android Socket编程一样，如果上来就加多线程、设计模式，会
让人摸不着头脑。但是如果从最简单的案例慢慢往上加功能，就很好理解

2024-02-26
TODO:
- 在GPUSolver2中添加存储旧ID的结构，实现从GPUID到ID的转换。(采取了另一种方式)
- 在以上基础上，改写write file函数 √
- 计算节点值 √

2024-02-20
添加GPU。
	TODO:
	重写输出函数，直接用ElementSoA计算node值，而不是先复制到Element_2D
	更具体地，由于事先无法确定需要输出哪些内容，计算过程应放在文件输出中
	计算Flux
测量内存带宽：参见书《并行计算与高性能计算》P274，文件位于D:\tgl\Local\HPC\GPU
进行至：void FVM_2D::calculateNodeValue_GPU

2024-02-09 ~ 
发现Limiter.cpp中Barth限制器的错误。
首先，限制器的除法好像除反了。并没有体现限制器的作用。
其次，限制器只需用于二阶精度格式，它用于消除色散引起的伪振荡，跟人工粘性的作用类似。

我现在采用的是AoS数据结构，即若干个数组并列。它参照CudaVS项目中的向量加法函数。不知道这样读取会不会影响性能。

2024-02-06 ~ 2024-02-08
完善FVM2D中solve_GPU函数。
接下来要编写GPU_calculateGradient()函数，参见GPU_space.h
修复代码错误：原CPU代码中“nNeighbor == 3”部分遗漏代码“dUm(2, 0) = dU[2][iVar]”(Reconstructor.cpp第110行)
进行至写__global__ void GPU::calculateGradientKernel(GPU::ElementDataPack& deviceDataPack)函数（2024-02-08）

2024-02-01
感觉科学计算不太适合用太多的面向对象。用纯过程和函数式编程更利于并行。
面向对象(OO)更适合桌面应用、GUI。OO是一种编程范式，它不能一招鲜吃遍天。
阅读OpenCFD-SCU，发现其cuda并行代码位于solver内计算通量的部分，这意味着我需要将Edge_2D的数据与方法分离，
将处理Edge_2D数据的方法放进新类中。
把所有cuda kelnel放进kernel.cu中，总共不超过200行，因为kernel都是简单的加减乘除运算

GPU编程的难点在于，它的所有kernel类似于函数式编程，仅知道当前编号。而非结构网格要操作的数据位置是不固定的。
例如进行规约操作，对于CPU编程，只需简单求和即可，但在GPU上需要二分法求和。

可能的解决办法：
在计算前，将所需数据(相邻单元数据、梯度、)整理成数组，将数组指针传给kernel。
虽然考虑到复制会耗费时间，但不试一试怎么知道最终效果如何呢。

鉴于GPU核心数上千（是CPU的数百倍），带宽大，可以将计算全部放进GPU

方案1：单元仅存最少数据：
nodeCoordinates[3][3]
values[4]
neighborValues[3][3]
先计算梯度
然后计算边界值
再计算边界数值通量(要用到黎曼求解器)
改变单元值

经过该步骤后，传输到主机，基于更新后的values，将所有单元的neighborValues更新，再传输进GPU
然后进行下一步循环

方案2：
nodeCoordinates[3][3]
fieldValues[4]
fieldValue


2024-01-27
pEdges向量在用push添加pEdge时，仅仅添加了
 
2024-01-26
今天编写SU2文件读取。需要重新理解以前的代码，进行至vBoundary2D
！新发现，SU2文件与ReadContinueFile的读取方式很像

2024-01-20
之前测试线性重构发现和常量重构精度差不多，原来是忘记将计算结果存进单元了
今天写了些限制器，但是线性重构发散的问题还没有解决。。。或许要添加防发散机制 是网格的问题？
问题描述：
  线性重构 isNan()检测到rho<0
  文件第26步出现-nan(ind)，共12处，为同一个单元的三个顶点的值。
  -nan(ind)出现的原因是分母为0或者负数开平方

  在isNan()中添加对rho<0的检测，未检测到
  在梯度计算器中添加对Ux Uy==nan的检测，未检测到
  更换为网格cylinder_2_41后，AMatrix报错：Error: the coefficient matrix has 0 as a pivot. Please fix the input and try again.
    共有两个异常点，分别位于左端和右端
    发生于element ID = 1549, x=1.98499, y=-0.917184
    如果忽视该警告，能够算下去。说明可能是网格的问题

2024-01-17
添加线性重构。
发现OpenMP，以后可以研究如何保证数据的原子性：https://blog.csdn.net/xuyiqiang87/article/details/52763653
粘性求解器进行至Solver_2D.cppP866，

2024-01-15
复杂的设计模式、不添加注释、奇怪的命名都会增加阅读代码的难度
先在TestCode项目中测试新建的TomlFileManager类
我想创建一种文件输出类，输入参数有：文件路径、文件的打开方式（覆盖还是添加）、内容

遇到的问题：
计算残差时，显示无穷。这是因为我将所有单元的绝对误差相加，超出了double的表示范围。正确做法是用相对误差表示，
或者取无穷范数（即元素最大值）

2024-01-14
新增toml格式文件读取（以前还自己造了个轮子，看来白造了），放在Include中
> 转换工具 https://tooltt.com/json2toml/
发现书籍：重构:改善现有代码的设计
不到万不得已不需要重构。或者先把最常用的重构了，一次重构一点
统计代码行数的技巧：在src目录下，右键git bash here，然后输入
```
find . -name "*.cpp" -or -name "*.h" | grep -v "include"  | xargs wc -l
```
> 参考资料：https://blog.csdn.net/jimojianghu/article/details/129792250

2024-01-13
之前创建了一些项目，可是时间太长，重新阅读很费劲。
于是仿照UNITS目录结构创建了本项目，方便后面的移植
UNITS项目`D:\tgl\Local\THUnits-Interface-UNITs\codes_UNITS_old\UNITs-v3.4-2021-12\UNITs_Solver`

当前目标：3维欧拉求解器
首先把FVM_2D的代码移植了过来，然后添加了git

难点：Solver_2D::getEdgeFlux_wallNonViscous，位于Solver_2D
Solver_2D::getEdgeFlux_farfield
Solver_2D::cal_ruvp_farfield_new cpp第825行
首先要理解Euler方程的数值边界条件处理。参见课堂笔记md文档
但我看王乾的代码，我确实已经实现了，为什么圆柱绕流的结果不理想？原来是边界条件问题，input.json速度设置为0
接下来是找标准算例验证自己程序的正确性。希望远场边界是正确的
此外添加四边形网格，修改程序框架


# 目录结构
include - 外部库

# 代码规范
文件用UTF-8格式保存

## 命名
目录名称用小写+下划线
文件名用驼峰命名法

# 注意事项
.h和.cpp分开是为了避免编译耗时太长
多用git

# Visual Studio 配置宏介绍
来源：https://blog.csdn.net/qianniulaoren/article/details/133160383
也可以在项目属性中点击任意项-编辑，展开宏，以进行浏览。
如果在常规的浏览中找不到，可以在调试的浏览中找


$(TargetDir)                                           # 目标输出文件所在的目录
$(TargetName) = U2NITS-10.0                            # 目标的名称
$(TargetExt) = .exe                                    # 目标的扩展名
$(TargetFileName) = $(TargetName)$(TargetExt)          # 目标输出文件名，包括扩展名
$(TargetPath) = $(TargetDir)$(TargetFileName)                      

$(SolutionDir) = D:\tgl\Local\HPC\U2NITS-10.0\         # 解决方案目录，即主目录
$(SolutionName) = U2NITS-10.0
$(SolutionExt) = .sln
$(SolutionFileName) = $(SolutionName)$(SolutionExt)    # .sln文件全名
$(SolutionPath) = $(SolutionDir)$(SolutionFileName)    # .sln文件全路径

$(Platform) = x64                                      # 解决方案平台名称，如x86、x64
$(Configuration) = Debug                               # 当前的编译配置名称，如Release、Debug
$(IntDir) = $(Platform)\$(Configuration) = x64\Debug\  # 编译器使用的中间目录，产出obj文件

$(ProjectName) = U2NITS-10.0                           # 当前工程名称
$(ProjectDir) = $(SolutionDir)$(ProjectName)\
              = D:\tgl\Local\HPC\U2NITS-10.0\U2NITS-10.0\ # 项目目录

$(OutDir) = $(SolutionDir)$(Platform)\$(Configuration)\
          = D:\tgl\Local\HPC\U2NITS-10.0\x64\Debug\    # 输出目录 该目录下生成exe pdb等文件

例如我将调试的工作目录修改为$(ProjectDir)WorkingDirectory\。设置后注意点“应用”

关于堆内存二维数组

void testHeap2DArray() {
	const int nVar = 4;// 守恒量个数，二维为4
	const int nX = 2;// 坐标个数，二维为2
	const int nVN = 3;// 有效邻居个数



	//real dUdX[nX][nVar]{};// 已初始化为0
	// 法一 指向数组的指针。delete时卡住，
	//auto dX = new real[nVN][nX];
	//auto dU = new real[nVN][nVar];
	// real (*dX)[2] 表示一个指向包含 2 个元素的 real 数组的指针，而 real * dX[2] 则是一个包含两个 real 指针的数组。
	/*
	释放时应使用delete[] dX。因为此时dX相当于一个一维数组，有nVN个元素，每个元素大小为2*sizeof(real)

	*/

	// 法一 动态二维数组 报错
	//real** dX = new real * [nVN];
	//real** dU = new real * [nVN];
	//for (int i = 0; i < nVN; i++) {
	//	dX[i] = new real[nX];
	//	dU[i] = new real[nVar];
	//}

	//// 法三 栈内存，能运行，但是大小固定
	//real dXs[nVN][nX];
	//real dUs[nVN][nVar];
	//auto dX = dXs;
	//auto dU = dUs;

	// 法四 全部当作一维数组
	
	real* dX = new real[nVN * nX];
	real* dUdX = new real[nX * nVar];
	real* dU = new real[nVN * nVar];

	real dXtrans[nX][nVN];
	

	for (int i = 0; i < nVN; i++) {
		for (int j = 0; j < nX; j++) {
			for (int k = 0; k < nVar; k++) {
				//dX[i][j] = i + j;
				//dUdX[j][k] = j + k;
				//dU[i][k] = i + k;

				dX[i * nX + j] = i + j;
				dUdX[j * nVar + k] = j + k;
				dU[i * nVar + k] = i + k;
			}
		}
	}
	GPU::Matrix::printMatrix(nVN, nX, (real*)dX);
	GPU::Matrix::printMatrix(nX, nVar, (real*)dUdX);
	GPU::Matrix::mul_ixj_jxk(nVN, nX, nVar, (real*)dX, (real*)dUdX, (real*)dU);
	GPU::Matrix::printMatrix(nVN, nVar, (real*)dU);

	//for (int i = 0; i < nVN; i++) {
	//	delete[] dX[i];
	//	delete[] dU[i];
	//}
	delete[] dX;
	delete[] dUdX;
	delete[] dU;
}


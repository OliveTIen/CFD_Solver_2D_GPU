---
date:   2024-01-12
author: tgl

---

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

# 日志

2024-01-13
之前创建了一些项目，可是时间太长，重新阅读很费劲。
UNITS项目参见 `D:\tgl\Local\THUnits-Interface-UNITs\codes_UNITS_old\UNITs-v3.4-2021-12\UNITs_Solver`

当前目标：3维欧拉求解器
首先把FVM_2D的代码移植了过来，然后添加了git

2024-01-14
新增toml格式文件读取（以前还自己造了个轮子，看来白造了），放在Include中
> 转换工具 https://tooltt.com/json2toml/
发现书籍：重构:改善现有代码的设计
不到万不得已不需要重构。或者先把最常用的重构了，一次重构一点
统计代码行数的技巧：在src目录下，右键git bash here，然后输入
``` shell
// 
find . -name "*.cpp" -or -name "*.h" | grep -v "include"  | xargs wc -l
// 统计fun3d共73万行代码
find . -name "*.f90" -or -name "*.cpp" -or -name "*.h" | grep -v "include"  | xargs wc -l

```

> 参考资料：https://blog.csdn.net/jimojianghu/article/details/129792250

2024-01-15
复杂的设计模式、不添加注释、奇怪的命名都会增加阅读代码的难度

2024-01-17
添加线性重构。
发现OpenMP，以后可以研究如何保证数据的原子性：https://blog.csdn.net/xuyiqiang87/article/details/52763653
粘性求解器进行至Solver_2D.cppP866，

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
     
2024-01-26
今天编写SU2文件读取。需要重新理解以前的代码，进行至vBoundary2D
[to do]readSU2File与ReadContinueFile很多代码高度重合，可以合并

2024-01-27
pEdges向量在用push添加pEdge时，仅仅添加了

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

2024-02-06
完善FVM2D中solve_GPU函数。
接下来要编写GPU_calculateGradient()函数，参见GPU_space.h
修复代码错误：原CPU代码中“nNeighbor == 3”部分遗漏代码“dUm(2, 0) = dU[2][iVar]”(Reconstructor.cpp第110行)
进行至写__global__ void GPU::calculateGradientKernel(GPU::ElementDataPack& deviceDataPack)函数（2024-02-08）

2024-02-09
发现Limiter.cpp中Barth限制器的错误。
首先，限制器的除法好像除反了。并没有体现限制器的作用。
其次，限制器只需用于二阶精度格式，它用于消除色散引起的伪振荡，跟人工粘性的作用类似。

2024-02-20
添加GPU。
	TODO:
	重写输出函数，直接用ElementSoA计算node值，而不是先复制到Element_2D
	更具体地，由于事先无法确定需要输出哪些内容，计算过程应放在文件输出中
	计算Flux
测量内存带宽：参见书《并行计算与高性能计算》P274，文件位于D:\tgl\Local\HPC\GPU
进行至：void FVM_2D::calculateNodeValue_GPU
我现在采用的是AoS数据结构，即若干个数组并列。它参照CudaVS项目中的向量加法函数。不知道这样读取会不会影响性能。

2024-02-26
TODO:
- 在GPUSolver2中添加存储旧ID的结构，实现从GPUID到ID的转换。(采取了另一种方式)
- 在以上基础上，改写write file函数 √
- 计算节点值 √


2024-03-08
为什么收敛这么慢？虽然算了接近1万步，但流场变化缓慢
GPU利用率低的原因：
1. 计算量太小，GPU并不能发挥作用
2. GPU内存带宽不足，限制了计算速率
试图模仿SU2进行重构，但是难度很大。这就像我学Android Socket编程一样，如果上来就加多线程、设计模式，会
让人摸不着头脑。但是如果从最简单的案例慢慢往上加功能，就很好理解

2024-03-14
性能分析工具
Visual Studio Profiler https://blog.csdn.net/sujunzy666/article/details/19963077

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

2024-03-19
在计算梯度的代码中，发现了LeastSquare的内存错误：nValidNeighbor提前+1。
以前直接改成GPU程序，导致该错误没发现，因此引发cuda_error 700异常。
看来在改为GPU程序前，一定要在CPU上测试一下

仍然出现异常。检查文件读取程序，发现element的指标管理混乱

2024-03-25
续算时，输出物理时间有误。可能是因为未读取


2024-04-01
关于枚举
在添加emitSignal()函数后，发现无法暂停，原来是SignalPack sp = emitSignalPack(istep, maxIteration, t, T, gpuSolver2.residualVector[0]);
这一语句中，sp.pauseSignal并没有被赋值为_timeReached。
最开始以为是enum类型的特性。最后发现是因为函数内部定义了一个局部变量PauseSignal，而成员变量没有被赋值。是盲目copy代码的锅。
enum的特性指的是如果要用int赋给它，需要强制类型转换，其它没区别。

2024-04-01
GPU报错 cuda error 701
输入` compute-sanitizer --launch-timeout=0 --tool=memcheck ./U2NITS-10.0.exe > opt.txt 2>&1`
生成的opt.txt中，
 Invalid __global__ read of size 4 bytes
原来是核函数kernel中读取了host的内存edge_device.num_edge。措施：添加成员int* _num_edge_，用cudaMemcpy申请内存
还是不行。我怀疑是不能用host定义的结构体

看来GPU部分要重来。
搜索CUDA封装，发现MUDA是CUDA的不错的封装，但是不是我想要的
https://zhuanlan.zhihu.com/p/659664377
它过于封装，看不见memcpy，以至于无法了解底层。可以用来应急(如果效率至上，能抓到老鼠都是好猫)
还是不行，报错。即使用C++20，还是报错：
1>D:\tgl\Local\HPC\U2NITS-10.0\TestCode\src\include\muda\muda\launch\kernel.h(92,1): fatal  error C1001: 内部编译器错误。

查看OpenCFD-SCU发现，其CUDA都是存储的全局变量

cuder也是一种封装

在TestCode中新建了CTestGPUSoA，用于测试CUDA编程，成功。看来核函数传参时只能传具体数据的指针，不能传ElementSOA的引用。引用也是指针，而global不允许传host指针。

cudaError701 too many resources requested for launch

2024-04-14
编程的原则：
尽量隐藏内部细节，用接口实现各层的通信。如果不隐藏，如果底层被修改，高层也要一并修改
目前的问题：
时间积分，显式太慢，用局部时间步又发散。隐式
求压力系数。我看UNITs中直接输出的表面压力 Cpp(I,J,K) = 2.0_8*(PP(I,J,K)-PPF)  
一种方式是积分表面压力，除以参考面积和参考动压
参考面积(参考弦长)可以直接输入，不用想太复杂。

发现了stdgpu，GPU上的类STL库 https://stotko.github.io/stdgpu/getting_started/integrating_into_your_project.html

2024-04-24 按2024-01-14的方法统计代码行数，为13825行

2024-05-03 待修复bug：续算readContinueFile时，若修改pause_的文件名中的数字，会找不到续算文件。不知道是因为小数位数的问题还是什么

2024-05-07 
上周完成了粘性（层流）。下面开始GPU
出现cuda error 701(too much resources required to launch)。假设是寄存器数量不够，将通量计算部分的block size调整为128，错误消失
但出现cuda error 700。用compute-sanitizer进行内存检查，具体代码为
```
compute-sanitizer --launch-timeout=0 --tool=memcheck ./U2NITS-10.0.exe > opt.txt 2>&1
```
最后发现，虽然getCudaLastError是在梯度计算被触发的，但实际上错误位于通量求解部分的代码，这也说明了GPU代码和CPU代码是错开执行的。
根据opt.txt，是因为FluxGPU.cu第298行，double* inf_ruvp = GlobalPara::boundaryCondition::_2D::inf::ruvp使用了指向host内存的指针，引发了内存的非法访问
添加cudaMalloc和cudaFree后，解决了该问题

PS：在排查error 701时，系统学习了CUDA编程的基础知识，包括内存模型、硬件架构等

2024-05-09
为了测试float，需要把所有double改为myfloat，但是cpptoml不支持float，会报错：static_assert failed "invalid value type"
于是进入cpptoml.h第282行，添加了float，但不知道有什么隐患

2024-05-13
GPU程序性能瓶颈为内存拷贝和释放，占用50%以上时间。为避免每次都拷贝，需要把calculate dt放在GPU上完成，这就涉及到规约运算
cuda-samples有两个项目是关于规约的，其中reductionMultiBlockCG是单文件项目，可以先看一看

(05-15)已经学习了GPU规约算法。计划：
在gpuSolver中添加数据结构，命名为elementFieldDt，存储与dt相关的物理量，例如每个单元的dt，计算dt所用的中间变量
改造成GPU代码，

(05-16)出现cuda error 13: invalid device symbol。跟规约程序对比后，将GPU::Math::reduce_device的最后一个参数取消使用，而直接
在函数体内从GPU::Math::p_operator_min复制到目标p_func_host

(05-17)无法将float *[4]转换为const float *[] (无法将 float ** 转化为 (const float) * *)
const float * * p 表示该指针的指向可以修改，指向的内容可以修改，指向再指向的内容不可修改
const float * const * p 表示该指针的指向不可修改，指向的内容可以修改，指向再指向的内容不可修改
const float const * * p 两个const位于float两侧，只相当于一个const，因此等价于const float * * p
这会导致const失效，因此编译器在编译时就报错阻止该危险行为
https://isocpp.org/wiki/faq/const-correctness#constptrptr-conversion

" * "也可以修饰const，https://blog.csdn.net/m0_37806112/article/details/81252151
https://zhuanlan.zhihu.com/p/90720012
读法上，从右向左结合。关键是const与 * 的相对位置，类型说明符char不重要
const char* p - p的指向可以修改，所指向的内容不可修改。(const (char(* p)))，首先读到(* p)，表示指针，然后读到char，即指向的内容为char，然后读到const，即指向的内容char不可修改
char const* p - 同上。(char (const(* p))，首先读到(* p)，表示是指针，然后读到const (* p)，即指向的内容(此时还不知道什么类型)不可修改，然后读到char，即指向的内容的类型为char
char* const p - p的指向不可修改(且必须初始化)，所指向的内容可以修改。首先读到const p，表示是常量(不可修改)，然后读到 * ，表示是指针，然后读到(char(* (const p)))，即指向的内容是char类型
const char* const p - p不可修改; 其指向的内容不可修改。从右向左读，即(const (char(* (const p))))。首先读到const，表示是常量，然后读到 * ，表示是指针，然后读到char，表示类型，然后读到const，表示指向的内容不可修改

对于(const int)类型变量例如const int a=10，则(&a)类型为 (const int* )

free
   关于free如何知道堆内存的大小：https://www.zhihu.com/question/23196195
   用alloc返回的指针，其前面有一块区域存储内存大小等信息
   因此堆指针可以互相交换，交换后各自free即可

(05-17)屡次出现out of bounds错误。用COMPUTE-SANITIZER检查发现是内存访问错误，根据以往经验应该是在device中使用了指向host内容的指针
于是尝试修改device代码，传参就传整个的elementField_device，而不是传elementField_device.ruvp等。问题解决
总结：
如果函数参数为myfloat* ruvp_device[4]，则会报错 cuda error，因为实际上你只传了一个二重指针myfloat** ruvp_device，函数认为ruvp_device[4]是指针数组，存储的数据都是连续的
指向分析：
ruvp_device在device内存中指向一个大小为4的数组，数组保存着4个地址，这些地址分别指向device内存中4个长数组
传参后，

疑问：为什么可以往核函数中传一维数组指针element_var.alphaC(参见CalculateDtGPU.cu: line 115: update_alphaC_device)，而不能传二维数组指针？
或许可以搜搜“给cuda核函数传递二维数组的方法”https://blog.csdn.net/lingerlanlan/article/details/24399909

(05-17)虽然程序能运行，但是仍然存在问题
device端计算dt最小值时，有20%的概率最小值略大于host函数的结果。例如host一直计算的是9.038e-5，device有时候是9.038e-5，有时候是1.073e-4等
可能是规约函数的问题。当n为奇数时，最后一个元素没有参与运算。
规约函数要求元素是2的幂次方的整数倍，否则在某一步会出现奇数个元素，此时末尾元素没有参与规约

该网址(https://www.cnblogs.com/5long/p/algorithms-on-cuda-reduction.html)说，“内核块的线程数量必须是2的幂次，否则计算结果是不正确的”

2024-05-19 迁移项目问题
准备迁移到另外一台机器运行，但是
cuda error the provided PTX was compiled with an unsupported toolchain
应该是机器的驱动版本太老。更新目标计算机驱动后解决问题

[https://blog.csdn.net/sophicchen/article/details/120782209]
命令行输入 nvcc -V，显示的是运行时版本(Runtime API)，因为编译器编译出的程序要跟着编译器版本走
命令行输入 nvidia-smi，显示的是驱动版本(Driver API)
要求驱动版本大于等于运行时版本，程序才能正常运行。

驱动下载网址：https://www.nvidia.cn/geforce/drivers/
驱动有两种 studio和game ready，

解决规约问题后，计算速度仍然很慢，看来要测量单独计算速度，而不要把输出的时长算进去

2024-05-28
如果发现绿色波浪线所描述的错误与实际不一致(明明错误已经解决但波浪线不消失)，可以删除"BROWSE.VC-1166caf"(位于`D:\Documents\VisualStudio_Temp`)

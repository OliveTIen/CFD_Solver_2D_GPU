---
date:   2024-01-12
author: tgl

---

# 日志
## 代办
1. 当务之急，输出残差曲线，查看收敛情况 ResidualCalculator中继续 √
2. 将一些文件放进目录层级中。例如把输入文件的类放进input中 √
3. 尝试课件的算例？？？可是它是非结构网格不可压，不知道是否需要另外写求解器
4. 尝试添加粘性(N-S)处理复杂边界
## 问题
1. 原来的项目中，输入square_per时出现内存问题
## 已解决
### ConsoleWindow问题
以前能够实现光标移动，可是更改代码后发现失效了。经排查，原来是开头输出input.toml时，由于输出内容超过1页，导致屏幕滚动，
然而目前的ConsolePrinter只支持一行滚动的
 
2024-01-17
添加线性重构。
发现OpenMP，以后可以研究如何保证数据的原子性：https://blog.csdn.net/xuyiqiang87/article/details/52763653

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
但我看王乾的代码，我确实已经实现了，为什么圆柱绕流的结果不理想？原来是我input.json速度设置为0
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

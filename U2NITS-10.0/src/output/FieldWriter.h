#ifndef FIELD_WRITER_H
#define FIELD_WRITER_H
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "../global/GlobalPara.h"
#include "../Node_2D.h"
#include "../Element_2D.h"
#include "../gpu/datatype/NodeSoA.h"
#include "../gpu/datatype/ElementSoA.h"
#include "../gpu/datatype/FieldSoA.h"
/*
一个没有数据的类，如果全部是静态成员函数，那么它与命名空间的效果相同。
*/

class FieldWriter {
private:
	static int numTecplotFileWritten;

public:
	FieldWriter() {};

	static void writeTecplotFile_1(
		double t_current,
		std::string filePath,
		std::string title,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::OutputNodeFieldSoA& output_node_field	);

	static void writeContinueFile_1(
		int i_step,
		double t_current,
		std::string filePath,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		real* elementField_U[4]);

	//static void writeContinueFile_1(
	//	int i_step,
	//	double t_current,
	//	std::string filePath,
	//	GPU::NodeSoA& nodes,
	//	GPU::ElementSoA& elements,
	//	real* element_U_old[4]) {
	//};


	// 添加了Ux Uy
	static void writeContinueFile_addUxUy(
		int i_step,
		double t_current,
		std::string filePath,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::FieldSoA& elementField);

	// 已输出的文件数量
	static int getNumTecplotFileWritten();
};

#endif
#ifndef FIELD_WRITER_H
#define FIELD_WRITER_H
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "../GlobalPara.h"
#include "../Node_2D.h"
#include "../Element_2D.h"
#include "../gpu/datatype/NodeSoA.h"
#include "../gpu/datatype/ElementSoA.h"
#include "../gpu/datatype/FieldSoA.h"
/*
一个没有数据的类，如果全部是静态成员函数，那么它与命名空间的效果相同。
*/

class FieldWriter {

public:
	FieldWriter() {};

	// 适合FVM_2D自带的数据结构
	static void writeTecplotFile(
		double t_current,
		std::string filePath,
		std::string title,
		const std::vector<Node_2D>& nodes, 
		const std::vector<Element_2D>& elements,
		const std::vector<double>& rho_nodes,
		const std::vector<double>& u_nodes,
		const std::vector<double>& v_nodes,
		const std::vector<double>& p_nodes		);
	static void writeTecplotFile_GPU(
		double t_current,
		std::string filePath,
		std::string title,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::OutputNodeFieldSoA& output_node_field	);
	static void writeContinueFile(
		int i_step,
		double t_current,
		std::string filePath,
		const std::vector<Node_2D>& nodes,
		const std::vector<Element_2D>& elements,
		void* pBoundaryManager);
	static void writeContinueFile_GPU(
		int i_step,
		double t_current,
		std::string filePath,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::FieldSoA& elementField,
		void* pBoundaryManager);

};

#endif
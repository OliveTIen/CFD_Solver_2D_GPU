#ifndef FIELD_WRITER_H
#define FIELD_WRITER_H
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "../GlobalPara.h"
#include "../Node_2D.h"
#include "../Element_2D.h"
/*
一个没有数据的类，如果全部是静态成员函数，那么它与命名空间的效果相同。


*/

class FieldWriter {
private:
	bool filePathHasInitialized = false;
	std::ofstream m_outfile;
	std::string m_filePath;

public:
	FieldWriter() {};

	FieldWriter(std::string f_name) {
		m_filePath = f_name;
		filePathHasInitialized = true;
	}

	void setFilePath(std::string f_name) {
		m_filePath = f_name;
		filePathHasInitialized = true;
	}
	// 适合FVM_2D自带的数据结构
	void writeTecplotFile(
		double t_current,
		std::string title,
		const std::vector<Node_2D>& nodes, 
		const std::vector<Element_2D>& elements,
		const std::vector<double>& rho_nodes,
		const std::vector<double>& u_nodes,
		const std::vector<double>& v_nodes,
		const std::vector<double>& p_nodes		);
	//// 计算节点场变量值ruvp，输出到文件
	//void writeTecplotFile(
	//);
	void writeContinueFile(
		int i_step,
		double t_current,
		std::string filePath,
		const std::vector<Node_2D>& nodes,
		const std::vector<Element_2D>& elements,
		void* pBoundaryManager);

private:
	//// 计算节点场变量值ruvp
	//void calculateNodeFieldValues(
	//	int num_node,
	//	double* rho,
	//	double* u,
	//	double* v,
	//	double* p);
};

/*
参考函数：
void FVM_2D::writeTecplotFile(std::string f_name, double t_current) {
	//std::cout << "wirtefile:" << f_name << std::endl;
	calculateNodeValue();
	std::ofstream f(f_name);
	if (!f.is_open())std::cout << "Error: Fail to open file " << f_name << std::endl;

	//header
	f << R"(TITLE = "Example: 2D Finite Element Data")" << "\n";
	f << R"(VARIABLES = "X", "Y")";
	if (GlobalPara::output::output_var_ruvp[0])		f << R"(, "rho")";
	if (GlobalPara::output::output_var_ruvp[1])		f << R"(, "u")";
	if (GlobalPara::output::output_var_ruvp[2])		f << R"(, "v")";
	if (GlobalPara::output::output_var_ruvp[3])		f << R"(, "p")";
	f << "\n";

	f << "ZONE NODES="<< nodes.size()
		<<", ELEMENTS = "<< elements.size()
		<<", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL"
		<<", SOLUTIONTIME="<<std::scientific<< t_current<<std::fixed
		<<std::endl;

	//nodes
	for (int i = 0; i < nodes.size(); i++) {
		f << nodes[i].x << ' ' << nodes[i].y;
		if (GlobalPara::output::output_var_ruvp[0])		f << ' ' << rho_nodes[i];
		if (GlobalPara::output::output_var_ruvp[1])		f << ' ' << u_nodes[i];
		if (GlobalPara::output::output_var_ruvp[2])		f << ' ' << v_nodes[i];
		if (GlobalPara::output::output_var_ruvp[3])		f << ' ' << p_nodes[i];
		f << std::endl;
	}

	//elements
	for (int ie = 0; ie < elements.size(); ie++) {
		f << elements[ie].nodes[0] << " ";
		f << elements[ie].nodes[1] << " ";
		f << elements[ie].nodes[2] << " ";
		f << elements[ie].nodes[2] << std::endl;
	}

	f.close();
}

*/

#endif
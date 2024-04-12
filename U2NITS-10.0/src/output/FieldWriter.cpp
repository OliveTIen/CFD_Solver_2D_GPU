#include "FieldWriter.h"
#include "../BoundaryManager_2D.h"
#include "../gpu/datatype/NodeSoA.h"
#include "../gpu/datatype/ElementSoA.h"
#include "../gpu/datatype/FieldSoA.h"
#include "../FVM_2D.h"

int FieldWriter::numTecplotFileWritten = 0;

void FieldWriter::writeTecplotFile_GPU(double t_current, std::string filePath, std::string title, GPU::NodeSoA& nodes, GPU::ElementSoA& elements, GPU::OutputNodeFieldSoA& output_node_field) {
	// Tecplot输出非结构网格时，节点编号从1开始
	// 目前可能存在问题的地方：
	// Element是三角形单元，因此输出0 1 2 2号节点，没有输出3号节点

	//// 判断路径是否初始化
	//if (!filePathHasInitialized) {
	//	std::cout << "Error: File path has not been initialized." << std::endl;
	//	return;
	//}
	//// 打开文件
	//std::ofstream& outfile = m_outfile;
	//outfile.open(m_filePath);
	//if (!outfile.is_open()) {
	//	std::cout << "Error: Fail to open file " << m_filePath << std::endl;
	//	return;
	//}
	// 打开文件
	std::ofstream outfile;
	outfile.open(filePath);
	if (!outfile.is_open())std::cout << "Error: Fail to open file " << filePath << std::endl;

	// 输出标题、变量名
	outfile << R"(TITLE = ")" << title << R"(")" << "\n";
	outfile << R"(VARIABLES = "X", "Y")";
	if (GlobalPara::output::output_var_ruvp[0])		outfile << R"(, "rho")";
	if (GlobalPara::output::output_var_ruvp[1])		outfile << R"(, "u")";
	if (GlobalPara::output::output_var_ruvp[2])		outfile << R"(, "v")";
	if (GlobalPara::output::output_var_ruvp[3])		outfile << R"(, "p")";
	outfile << "\n";
	// 域定义
	outfile << "ZONE NODES=" << nodes.num_node
		<< ", ELEMENTS = " << elements.num_element
		<< ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL"
		<< ", SOLUTIONTIME=" << std::scientific << t_current << std::fixed
		<< std::endl;

	// 节点
	for (int i = 0; i < nodes.num_node; i++) {
		outfile << nodes.xy[0][i] << ' ' << nodes.xy[1][i];
		REAL& rho = output_node_field.ruvp[0][i];
		REAL& p = output_node_field.ruvp[3][i];
		if (rho < 0) {
			//std::cout<<"rho<0"<<std::endl;
		}
		if (GlobalPara::output::output_var_ruvp[0])		outfile << ' ' << output_node_field.ruvp[0][i];
		if (GlobalPara::output::output_var_ruvp[1])		outfile << ' ' << output_node_field.ruvp[1][i];
		if (GlobalPara::output::output_var_ruvp[2])		outfile << ' ' << output_node_field.ruvp[2][i];
		if (GlobalPara::output::output_var_ruvp[3])		outfile << ' ' << output_node_field.ruvp[3][i];
		outfile << std::endl;
	}

	// 单元
	for (int ie = 0; ie < elements.num_element; ie++) {
		outfile << elements.nodes[0][ie] + 1 << " ";
		outfile << elements.nodes[1][ie] + 1 << " ";
		outfile << elements.nodes[2][ie] + 1 << " ";
		outfile << elements.nodes[2][ie] + 1 << std::endl;
	}

	outfile.close();
	numTecplotFileWritten++;
}

void FieldWriter::writeContinueFile_GPU(int i_step, double t_current, std::string filePath, GPU::NodeSoA& nodes, GPU::ElementSoA& elements, GPU::FieldSoA& elementField) {
	// 输出暂存文件，用于下次续算
	// 目前还不敢直接用GPUID

    // 变量定义
	std::ofstream outfile;
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	BoundaryManager_2D& boundaryManager = pFVM2D->boundaryManager;

	// 打开文件
	outfile.open(filePath);
	if (!outfile.is_open())std::cout << "Error: Fail to open file " << filePath << std::endl;
	// 输出时间步数
	outfile << "t, istep" << std::endl;
	outfile << t_current << " " << i_step << std::endl;
	// 输出节点
	outfile << "nodes: ID, x, y" << std::endl;
	for (int in = 0; in < nodes.num_node; in++) {
		outfile << pFVM2D->nodes[in].ID << " ";// ID
		//outfile << nodes.ID[in] << " ";// GPUID
		outfile << nodes.xy[0][in] << " ";
		outfile << nodes.xy[1][in] << std::endl;
	}

	outfile << "elements: ID, nodes, U" << std::endl;
	for (int ie = 0; ie < elements.num_element; ie++) {
		outfile << pFVM2D->elements[ie].ID << " ";// ID
		outfile << pFVM2D->elements[ie].nodes[0] << " ";
		outfile << pFVM2D->elements[ie].nodes[1] << " ";
		outfile << pFVM2D->elements[ie].nodes[2] << " ";
		//outfile << elements.ID[ie] << " ";//GPUID
		//outfile << elements.nodes[0][ie] << " ";
		//outfile << elements.nodes[1][ie] << " ";
		//outfile << elements.nodes[2][ie] << " ";
		outfile << elementField.U[0][ie] << " ";
		outfile << elementField.U[1][ie] << " ";
		outfile << elementField.U[2][ie] << " ";
		outfile << elementField.U[3][ie] << std::endl;
	}

	outfile << "boundaries: ID, name; edgeIDs" << std::endl;
	for (int ib = 0; ib < boundaryManager.boundaries.size(); ib++) {
		outfile << boundaryManager.boundaries[ib].ID << " ";
		outfile << boundaryManager.boundaries[ib].name << std::endl;
		//一行太长会导致读取时buffer长度不够，从而导致难以预料的结果(死循环、读取失败等)
		int nEdge = (int)boundaryManager.boundaries[ib].pEdges.size();//控制输出' '还是'\n'
		bool check = 0;//控制
		for (int iEdge = 0; iEdge < nEdge; iEdge++) {
			outfile << boundaryManager.boundaries[ib].pEdges[iEdge]->ID;
			if (iEdge % 10 == 9) {
				outfile << std::endl;
				check = 1;
			}
			else {
				outfile << ' ';
				check = 0;
			}
		}
		if (!check)outfile << std::endl;//若上次输出的' '
	}
	outfile.close();

}

void FieldWriter::writeContinueFile_2(int i_step, double t_current, std::string filePath, GPU::NodeSoA& nodes, GPU::ElementSoA& elements, GPU::FieldSoA& elementField) {
	// 添加Ux Uy

	// 变量定义
	std::ofstream outfile;
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	BoundaryManager_2D& boundaryManager = pFVM2D->boundaryManager;

	// 打开文件
	outfile.open(filePath);
	if (!outfile.is_open())std::cout << "Error: Fail to open file " << filePath << std::endl;
	// 输出时间步数
	outfile << "t, istep" << std::endl;
	outfile << t_current << " " << i_step << std::endl;
	// 输出节点
	outfile << "nodes: ID, x, y" << std::endl;
	for (int in = 0; in < nodes.num_node; in++) {
		outfile << pFVM2D->nodes[in].ID << " ";// ID
		//outfile << nodes.ID[in] << " ";// GPUID
		outfile << nodes.xy[0][in] << " ";
		outfile << nodes.xy[1][in] << std::endl;
	}

	outfile << "elements: ID, nodes[3], U[4], Ux[4], Uy[4]" << std::endl;
	for (int ie = 0; ie < elements.num_element; ie++) {
		outfile << pFVM2D->elements[ie].ID << " ";// ID
		outfile << pFVM2D->elements[ie].nodes[0] << " ";
		outfile << pFVM2D->elements[ie].nodes[1] << " ";
		outfile << pFVM2D->elements[ie].nodes[2] << " ";
		outfile << elementField.U[0][ie] << " ";
		outfile << elementField.U[1][ie] << " ";
		outfile << elementField.U[2][ie] << " ";
		outfile << elementField.U[3][ie] << " ";
		outfile << elementField.Ux[0][ie] << " ";
		outfile << elementField.Ux[1][ie] << " ";
		outfile << elementField.Ux[2][ie] << " ";
		outfile << elementField.Ux[3][ie] << " ";
		outfile << elementField.Uy[0][ie] << " ";
		outfile << elementField.Uy[1][ie] << " ";
		outfile << elementField.Uy[2][ie] << " ";
		outfile << elementField.Uy[3][ie];
		outfile << std::endl;
	}

	outfile << "boundaries: ID, name; edgeIDs" << std::endl;
	for (int ib = 0; ib < boundaryManager.boundaries.size(); ib++) {
		outfile << boundaryManager.boundaries[ib].ID << " ";
		outfile << boundaryManager.boundaries[ib].name << std::endl;
		//一行太长会导致读取时buffer长度不够，从而导致难以预料的结果(死循环、读取失败等)
		int nEdge = (int)boundaryManager.boundaries[ib].pEdges.size();//控制输出' '还是'\n'
		bool check = 0;//控制
		for (int iEdge = 0; iEdge < nEdge; iEdge++) {
			outfile << boundaryManager.boundaries[ib].pEdges[iEdge]->ID;
			if (iEdge % 10 == 9) {
				outfile << std::endl;
				check = 1;
			}
			else {
				outfile << ' ';
				check = 0;
			}
		}
		if (!check)outfile << std::endl;//若上次输出的' '
	}
	outfile.close();
}

int FieldWriter::getNumTecplotFileWritten() {
	return numTecplotFileWritten; 
}



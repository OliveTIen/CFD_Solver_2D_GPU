#include "FieldWriter.h"
#include "../BoundaryManager_2D.h"
#include "../gpu/datatype/NodeSoA.h"
#include "../gpu/datatype/ElementSoA.h"
#include "../gpu/datatype/FieldSoA.h"
#include "../FVM_2D.h"

void FieldWriter::writeTecplotFile(double t_current, std::string title, const std::vector<Node_2D>& nodes, const std::vector<Element_2D>& elements, const std::vector<double>& rho_nodes, const std::vector<double>& u_nodes, const std::vector<double>& v_nodes, const std::vector<double>& p_nodes) {
	// 这里读写时用的是旧ID，因此需要旧的nodes和elements
	if (!filePathHasInitialized) {
		std::cout << "Error: File path has not been initialized." << std::endl;
		return;
	}
	std::ofstream& f = m_outfile;
	f.open(m_filePath);
	if (!f.is_open()) {
		std::cout << "Error: Fail to open file " << m_filePath << std::endl;
		return;
	}

	//header
	f << R"(TITLE = ")" << title << R"(")" << "\n";
	f << R"(VARIABLES = "X", "Y")";
	if (GlobalPara::output::output_var_ruvp[0])		f << R"(, "rho")";
	if (GlobalPara::output::output_var_ruvp[1])		f << R"(, "u")";
	if (GlobalPara::output::output_var_ruvp[2])		f << R"(, "v")";
	if (GlobalPara::output::output_var_ruvp[3])		f << R"(, "p")";
	f << "\n";

	f << "ZONE NODES=" << nodes.size()
		<< ", ELEMENTS = " << elements.size()
		<< ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL"
		<< ", SOLUTIONTIME=" << std::scientific << t_current << std::fixed
		<< std::endl;

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

void FieldWriter::writeTecplotFile_GPU(double t_current, std::string title, GPU::NodeSoA& nodes, GPU::ElementSoA& elements, GPU::OutputNodeFieldSoA& output_node_field) {
	// Tecplot输出非结构网格时，节点编号从1开始
	// 目前可能存在问题的地方：
	// Element是三角形单元，因此输出0 1 2 2号节点，没有输出3号节点

	// 判断路径是否初始化
	if (!filePathHasInitialized) {
		std::cout << "Error: File path has not been initialized." << std::endl;
		return;
	}
	// 打开文件
	std::ofstream& f = m_outfile;
	f.open(m_filePath);
	if (!f.is_open()) {
		std::cout << "Error: Fail to open file " << m_filePath << std::endl;
		return;
	}

	// 输出标题、变量名
	f << R"(TITLE = ")" << title << R"(")" << "\n";
	f << R"(VARIABLES = "X", "Y")";
	if (GlobalPara::output::output_var_ruvp[0])		f << R"(, "rho")";
	if (GlobalPara::output::output_var_ruvp[1])		f << R"(, "u")";
	if (GlobalPara::output::output_var_ruvp[2])		f << R"(, "v")";
	if (GlobalPara::output::output_var_ruvp[3])		f << R"(, "p")";
	f << "\n";
	// 域定义
	f << "ZONE NODES=" << nodes.num_node
		<< ", ELEMENTS = " << elements.num_element
		<< ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL"
		<< ", SOLUTIONTIME=" << std::scientific << t_current << std::fixed
		<< std::endl;

	// 节点
	for (int i = 0; i < nodes.num_node; i++) {
		f << nodes.xy[0][i] << ' ' << nodes.xy[1][i];
		if (GlobalPara::output::output_var_ruvp[0])		f << ' ' << output_node_field.ruvp[0][i];
		if (GlobalPara::output::output_var_ruvp[1])		f << ' ' << output_node_field.ruvp[1][i];
		if (GlobalPara::output::output_var_ruvp[2])		f << ' ' << output_node_field.ruvp[2][i];
		if (GlobalPara::output::output_var_ruvp[3])		f << ' ' << output_node_field.ruvp[3][i];
		f << std::endl;
	}

	// 单元
	for (int ie = 0; ie < elements.num_element; ie++) {
		f << elements.nodes[0][ie] + 1 << " ";
		f << elements.nodes[1][ie] + 1 << " ";
		f << elements.nodes[2][ie] + 1 << " ";
		f << elements.nodes[2][ie] + 1 << std::endl;
	}

	f.close();
}

void FieldWriter::writeContinueFile(
	int i_step,
	double t_current,
	std::string filePath,
	const std::vector<Node_2D>& nodes,
	const std::vector<Element_2D>& elements,
	void* pBoundaryManager) {

	// 输出暂存文件，用于下次续算
	// 变量定义
	std::ofstream& outfile = m_outfile;
	BoundaryManager_2D& boundaryManager = *(BoundaryManager_2D*)pBoundaryManager;
	// 打开文件
	outfile.open(filePath);
	if (!outfile.is_open())std::cout << "Error: Fail to open file " << filePath << std::endl;
	// 输出时间步数
	outfile << "t, istep" << std::endl;
	outfile << t_current << " " << i_step << std::endl;
	// 输出节点
	outfile << "nodes: ID, x, y" << std::endl;
	for (int in = 0; in < nodes.size(); in++) {
		//// 该函数需要修改，因为是根据单元值和梯度计算节点值。事实上我们已经计算了节点值
		//std::vector<double> uu_i = nodes[in].calNodeValue();
		outfile << nodes[in].ID << " ";
		outfile << nodes[in].x << " ";
		outfile << nodes[in].y << std::endl;
	}

	outfile << "elements: ID, nodes, U" << std::endl;
	for (int ie = 0; ie < elements.size(); ie++) {
		outfile << elements[ie].ID << " ";
		outfile << elements[ie].nodes[0] << " ";
		outfile << elements[ie].nodes[1] << " ";
		outfile << elements[ie].nodes[2] << " ";
		outfile << elements[ie].U[0] << " ";
		outfile << elements[ie].U[1] << " ";
		outfile << elements[ie].U[2] << " ";
		outfile << elements[ie].U[3] << std::endl;
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

void FieldWriter::writeContinueFile_GPU(int i_step, double t_current, std::string filePath, GPU::NodeSoA& nodes, GPU::ElementSoA& elements, GPU::FieldSoA& elementField, void* pBoundaryManager) {
	// 输出暂存文件，用于下次续算
	// 目前还不敢直接用GPUID

    // 变量定义
	std::ofstream& outfile = m_outfile;
	BoundaryManager_2D& boundaryManager = *(BoundaryManager_2D*)pBoundaryManager;
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;

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



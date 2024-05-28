#include "FieldWriter.h"
#include "../BoundaryManager_2D.h"
#include "../gpu/datatype/NodeSoA.h"
#include "../gpu/datatype/ElementSoA.h"
#include "../gpu/datatype/FieldSoA.h"
#include "../FVM_2D.h"
#include "LogWriter.h"
#include "../global/CExit.h"
#include "../math/Common.h"
#include "../input/TomlFileManager.h"
#include "../solvers/GPUSolver2.h"
#include "../solvers/SolverDataGetter.h"

FieldWriter* FieldWriter::p_instance = nullptr;

FieldWriter* FieldWriter::getInstance() {
	if (p_instance == nullptr) {
		p_instance = new FieldWriter();
	}
	return p_instance;
}



// 加引号
inline std::string Quote(const std::string s) {
	return "\"" + s + "\"";
}

// 将stringVector转换为string，用", "连接，且各元素加上双引号
std::string stringVector_to_string_withQuote(const std::vector<std::string>& v) {
	std::string out;
	size_t size = v.size();
	for (int i = 0; i < size; i++) {
		out = out + Quote(v[i]);

		if (i != size - 1) {
			out = out + ", ";
		}
	}
	return out;
}

void FieldWriter::write_tecplot_volume_file(myfloat t_current, std::string filePath, std::string title, GPU::NodeSoA& nodes, GPU::ElementSoA& elements, GPU::OutputNodeFieldSoA& output_node_field) {

/*
参照Tecplot安装目录 E:\Tecplot\Tecplot\Tecplot 360 EX 2023 R1\doc\360_data_format_guide.pdf

zonetype:
  fetriangle: 三角元 fe triangle
  fequadrilateral: 四边元
  fetetrahedron: 四面体
  febrick: 六面体
datapacking:
  point: 格点
  block: 格心
*/

	// 打开文件
	std::ofstream outfile;
	outfile.open(filePath);
	if (!outfile.is_open()) {
		LogWriter::logAndPrintError("failed to open file " + filePath +" @write_tecplot_volume_file\n");
		CExit::saveAndExit(-1);
	}

	std::string zone_title = "volume field";
	outfile 
	// 标题、变量名
		<< "title= " << Quote(title) << "\n"
		<< "variables= " << m_outputScheme_volume.get_variable_names() << "\n"
	// 域定义
		<< "zone "
		<< "t= " << Quote(zone_title) << ", "
		<< "nodes= " << nodes.num_node << ", "
		<< "elements= " << elements.num_element << ", "
		<< "datapacking= " << "point" << ", "
		<< "zonetype= " << "fequadrilateral" << ", "
		<< "solutiontime= " << std::scientific << t_current << std::fixed << std::endl;

	// 节点
	for (int i = 0; i < nodes.num_node; i++) {
		outfile << m_outputScheme_volume.get_variable_value_string(i, nodes, output_node_field) << "\n";
	}

	// 单元
	for (int ie = 0; ie < elements.num_element; ie++) {
		outfile << elements.nodes[0][ie] + 1 << " ";
		outfile << elements.nodes[1][ie] + 1 << " ";
		outfile << elements.nodes[2][ie] + 1 << " ";
		outfile << elements.nodes[2][ie] + 1 << std::endl;
	}

	outfile.close();
	m_numTecplotFileWritten++;
}

void FieldWriter::write_tecplot_boundary_file(myfloat t_current,std::string filePath,GPU::NodeSoA& node_host,GPU::EdgeSoA& edge_host,GPU::OutputNodeFieldSoA& output_node_field,GPU::BoundaryV2& boundary_host_new) {
	/*
	由于流场是二维的，因此边界是一维的，要绘制一维线图
	*/
	std::ofstream outfile;
	outfile.open(filePath);
	if (outfile.fail()) {
		LogWriter::logAndPrintError("failed to open file " + filePath + " @write_tecplot_boundary_file\n");
		return;
	}

	std::string title = "tecplot boundary file";
	outfile
		// 标题、变量名
		<< "title= " << Quote(title) << "\n"
		<< "variables= " << m_outputScheme_boundary.get_variable_names() << "\n";
	
	const auto& edgeSets = boundary_host_new.edgeSets;
	//m_referenceData.CD_vector = std::vector<myfloat>(edgeSets.size(), 0.0);// 所有edgeSet的阻力。值初始化为0
	//m_referenceData.CL_vector = std::vector<myfloat>(edgeSets.size(), 0.0);
	const myfloat S_ref = GlobalPara::constant::referenceArea;// 参考面积
	const myfloat one_on_Sref = 1 / S_ref;
	for (size_t iSet = 0; iSet < edgeSets.size(); iSet++) {
		const auto& edgeSet = edgeSets[iSet];// 当前边集合
		std::string zone_title = "boundary " + std::to_string(edgeSet.ID) + " " + edgeSet.name;// boundary 1 far
		myint boundary_edge_num = edgeSet.edge_vector.size();// 当前edge set的edge数量
		outfile
			<< std::scientific
			<< "zone "
			<< "t= " << Quote(zone_title) << ", "
			<< "solutiontime= " << t_current << ", "
			<< "strandid= " << 0 << ", "
			<< "i= " << boundary_edge_num << ", "
			<< "j= " << 1 << ", "
			<< "f= " << "point" << "\n";
		myint node0_old = -1;
		myint node1_old = -1;
		for (myint iEdge = 0; iEdge < boundary_edge_num; iEdge++) {
			/*
			输出当前edgeSet的边界点数据
			通常，同一个edgeSet中相邻两个edge是连续的，因此前者的node1等于后者的node0
			但也会有特殊情况，导致输出混乱
			*/
			myint edgeID = edgeSet.edge_vector[iEdge];
			myint node0 = edge_host.nodes[0][edgeID];
			myint node1 = edge_host.nodes[1][edgeID];

			if (node0 != node1_old) {
				outfile << m_outputScheme_boundary.get_variable_value_string(node0, node_host, output_node_field) << "\n";
			}
			outfile << m_outputScheme_boundary.get_variable_value_string(node1, node_host, output_node_field) << "\n";

			node0_old = node0;
			node1_old = node1;

			///*
			//计算升阻力，由边界压力积分得到
			//公式：
			//dCD = - Cp / S_ref * dS * nx
			//dCL = - Cp / S_ref * dS * ny
			//然后对边界进行环路积分
			//*/
			//
			//myint nx = edge_host.normal[0][edgeID];
			//myint ny = edge_host.normal[1][edgeID];
			//myint dS = edge_host.length[edgeID];// 边长度
			//myint Cp = 0.5 * (output_node_field.Cp[node0] + output_node_field.Cp[node1]);
			//myint dC_tmp = (-Cp) * one_on_Sref * dS;// 中间变量，为了减少计算量而设置，无物理含义
			//m_referenceData.CD_vector[iSet] += dC_tmp * nx;
			//m_referenceData.CL_vector[iSet] += dC_tmp * ny;
		}
	}

	outfile.close();
}

void FieldWriter::ReferenceData::calculate_edgeSet_force(GPU::EdgeSoA& edge_host, GPU::OutputNodeFieldSoA& output_node_field, GPU::BoundaryV2& boundary_host_new) {
	/*
	遍历边界边，根据压力系数积分得到升阻力
	公式：
	dCD = Cp / S_ref * dS * nx
	dCL = Cp / S_ref * dS * ny
	*/
	const auto& edgeSets = boundary_host_new.edgeSets;
	const myfloat S_ref = GlobalPara::constant::referenceArea;// 参考面积
	const myfloat one_on_Sref = 1 / S_ref;

	// 初始化向量size。初始化值为0
	this->CD_vector = std::vector<myfloat>(edgeSets.size(), 0.0);
	this->CL_vector = std::vector<myfloat>(edgeSets.size(), 0.0);

	for (size_t iSet = 0; iSet < edgeSets.size(); iSet++) {
		const auto& edgeSet = edgeSets[iSet];// 当前edgeSet
		myint boundary_edge_num = edgeSet.edge_vector.size();// 当前edgeSet的edge数量
		for (myint iEdge = 0; iEdge < boundary_edge_num; iEdge++) {
			myint edgeID = edgeSet.edge_vector[iEdge];
			myint node0 = edge_host.nodes[0][edgeID];
			myint node1 = edge_host.nodes[1][edgeID];
			myfloat nx = edge_host.normal[0][edgeID];
			myfloat ny = edge_host.normal[1][edgeID];
			myfloat dS = edge_host.length[edgeID];// 边长度
			myfloat Cp = 0.5 * (output_node_field.Cp[node0] + output_node_field.Cp[node1]);
			myfloat dC_tmp = Cp * one_on_Sref * dS;// 中间变量，为了减少计算量而设置，无物理含义
			this->CD_vector[iSet] += dC_tmp * nx;
			this->CL_vector[iSet] += dC_tmp * ny;
		}
	}
}

void FieldWriter::write_tecplot_hist_file(std::string filePath, int iteration, myfloat t_physics, myfloat residual[4], GPU::EdgeSoA& edge_host, GPU::OutputNodeFieldSoA& output_node_field, GPU::BoundaryV2& boundary_host_new) {
	/*
	判断是否应该写文件头
	在第一次调用该函数的情况下：
	若hist文件不存在，则需要写文件头
	若GlobalPara::basic::_continue==false，则需要写文件头
	*/
	bool b_should_write_head = false;
	if (!m_called_write_tecplot_hist_file) {
		std::ifstream infile(filePath, std::ios_base::in);// 只读
		if (!infile.is_open()) b_should_write_head = true;
		if (GlobalPara::basic::_continue == false)b_should_write_head = true;
	}
	m_called_write_tecplot_hist_file = true;

	std::ofstream outfile;
	if (b_should_write_head) {
		outfile.open(filePath, std::ios_base::app);
		if (outfile.fail()) {
			LogWriter::logAndPrintWarning("Write hist file head failed.\n");
		}

		std::string title = "tecplot hist file";
		outfile
			// 标题、变量名
			<< "title= " << Quote(title) << "\n"
			<< "variables= " << m_outputScheme_hist.get_variable_names() << "\n";

		std::string zone_title = "zone title";
		outfile
			<< std::scientific
			<< "zone "
			<< "t= " << Quote(zone_title) << ", "
			<< "f= " << "point" << "\n";
		/*
		outfile << R"(TITLE=")" << GlobalPara::basic::filename << R"(")" << "\n";
		outfile << R"(VARIABLES="Iteration" "Residual_rho" "Residual_rhou" "Residual_rhov" "Residual_rhop")" << "\n";
		outfile << R"(ZONE  F=POINT)" << "\n";		
		*/

		outfile.close();
	}

	outfile.open(filePath, std::ios_base::app);
	if (outfile.fail()) {
		LogWriter::logAndPrintWarning("Write hist file failed.\n");
	}

	m_referenceData.calculate_edgeSet_force(edge_host, output_node_field, boundary_host_new);
	m_histData.update(iteration, t_physics, residual, m_referenceData, SolverDataGetter::getSolverInstance()->boundary_host_new);
	outfile << m_outputScheme_hist.get_variable_value_string(m_histData) << "\n";

	outfile.close();

}


void FieldWriter::writeContinueFile_1(int i_step, myfloat t_current, std::string filePath, GPU::NodeSoA& nodes, GPU::ElementSoA& elements, myfloat* elementField_U[4]) {
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
	outfile << "t, istep, CFL" << std::endl;
	outfile << t_current << " " << i_step << " "<< GlobalPara::time::CFL << std::endl;// CFL是自定义CFL策略中需要读取的，没有该策略读取了也没用
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
		outfile << elementField_U[0][ie] << " ";
		outfile << elementField_U[1][ie] << " ";
		outfile << elementField_U[2][ie] << " ";
		outfile << elementField_U[3][ie] << std::endl;
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

void FieldWriter::allocNodeFieldDataUsingOutputScheme(GPU::OutputNodeFieldSoA& nodeField, myint num_node) {
	/*
	根据输出方案申请内存。需要在读取完Config和Mesh后调用，因为要知道输出方案和节点个数
	*/
	// 判断是否已读取Config
	if (!m_initialized) {
		LogWriter::logAndPrintError("outputScheme has not been initialized. @FieldWriter::allocNodeFieldDataUsingOutputScheme\n");
		exit(-1);
	}
	// 判断是否读取网格Mesh。若num_node==0表示未读取网格
	if (num_node == 0) {
		LogWriter::logAndPrintError("num_node == 0. @FieldWriter::allocNodeFieldDataUsingOutputScheme\n");
		exit(-1);
	}
	//判断是否重复申请内存
	if (nodeField.b_ruvp_allocated) {
		LogWriter::logAndPrintError("nodeField.ruvp has been allocated. @FieldWriter::allocNodeFieldDataUsingOutputScheme\n");
		exit(-1);
	}

	nodeField.alloc_ruvp(num_node);
	
	if (m_outputScheme_boundary.Cp || m_outputScheme_volume.Cp) {
		nodeField.Cp.resize(num_node);
	}
	if (m_outputScheme_boundary.Ma || m_outputScheme_volume.Ma) {
		nodeField.Ma.resize(num_node);
	}

	nodeField.b_all_allocated = true;
}

void FieldWriter::freeNodeFieldData(GPU::OutputNodeFieldSoA& nodeField) {
	nodeField.free_ruvp();
}

void FieldWriter::update_nodeField() {
	GPU::GPUSolver2* solver = GPU::GPUSolver2::getInstance();
	//update_nodeField_ruvp(solver->element_host, solver->outputNodeField, solver->element_vruvp);
	update_nodeField_ruvp(solver->element_host, solver->outputNodeField, solver->elementField_host.ruvp);
	update_nodeField_other_variables(solver->outputNodeField);
}

void FieldWriter::update_nodeField_ruvp(GPU::ElementSoA& element, GPU::OutputNodeFieldSoA& nodeField, myfloat* element_vruvp[4]) {
	/*
	将格心数据element_vruvp转化为格点数据nodeField的ruvp

	可能的隐患：
	三角单元用四个点存储，一般令nodeID[3]=-1。但若nodeID[3]==nodeID[2]，会导致nodeID[2]权重偏大
	*/
	auto& node_ruvp = nodeField.ruvp;
	myint num_node = nodeField.num_node;
	myint num_element = element.num_element;
	const int nNodePerElement = 4;
	const int nValuePerNode = 4;
	int* node_neighborElement_num;// 记录每个节点的邻居单元数量
	// 申请资源 初始化
	node_neighborElement_num = new int[num_node]();// 小括号自动初始化为0
	for (int i = 0; i < 4; i++) {
		memset(node_ruvp[i], 0, num_node * sizeof(myfloat));
	}

	// 将单元值加到其子节点值上
	for (int iElement = 0; iElement < num_element; iElement++) {
		// 统计每个节点的邻居单元数量，并修改节点值
		for (int jNode = 0; jNode < nNodePerElement; jNode++) {
			// 获取单元节点ID
			int GPUID_of_node = element.nodes[jNode][iElement];// GPUID of node 0
			// nodeID[3]=-1时，跳过该节点
			if (GPUID_of_node < 0 || GPUID_of_node >= num_node)continue;// 跳过本次循环，但不跳出循环体
			// 该ID对应的邻居单元数量+1
			node_neighborElement_num[GPUID_of_node]++;
			// 该ID对应的节点所有值加上单元值
			for (int kValue = 0; kValue < nValuePerNode; kValue++) {
				node_ruvp[kValue][GPUID_of_node] += element_vruvp[kValue][iElement];
			}
		}
	}

	// 节点值除以邻居单元数，得到平均值，作为节点ruvp值
	for (int iNode = 0; iNode < num_node; iNode++) {
		// 为了避免除以0，分母等于0则跳过
		if (node_neighborElement_num[iNode] == 0)continue;
		// node_ruvp除以邻居单元数，得到平均值
		for (int kValue = 0; kValue < nValuePerNode; kValue++) {
			node_ruvp[kValue][iNode] /= node_neighborElement_num[iNode];
		}
	}

	// 释放资源
	delete[] node_neighborElement_num;
}

void FieldWriter::update_nodeField_other_variables(GPU::OutputNodeFieldSoA& nodeField) {
	// 检查是否申请内存
	if (!nodeField.b_all_allocated) {
		LogWriter::logAndPrintError("nodeField has not been fully allocated.\n");
		exit(-1);
	}

	bool b_Cp = m_outputScheme_boundary.Cp || m_outputScheme_volume.Cp;
	bool b_Ma = m_outputScheme_boundary.Ma || m_outputScheme_volume.Ma;
	myint num_node = nodeField.num_node;
	myfloat CD = 0.0;
	myfloat CL = 0.0;
	for (myint i = 0; i < num_node; i++) {
		myfloat rho = nodeField.ruvp[0][i];
		myfloat u = nodeField.ruvp[1][i];
		myfloat v = nodeField.ruvp[2][i];
		myfloat p = nodeField.ruvp[3][i];
		if (b_Cp) {
			// 升力系数 Cp = (p - p_inf)/p_dynamic_inf
			m_referenceData.update_farfield();
			nodeField.Cp[i] = (p - m_referenceData.p_inf) / m_referenceData.p_dynamic_inf;
			// 升阻力 见write_tecplot_boundary_file
		}
		if (b_Ma) {
			myfloat V2 = u * u + v * v;
			myfloat a2 = GlobalPara::constant::gamma * p / rho;
			myfloat Ma2 = V2 / a2;
			nodeField.Ma[i] = sqrt(Ma2);
		}
	}
	
}


void FieldWriter::initialize_outputScheme_usingConfig() {
	// 初始化输出方案，只在readConfig时调用
	TomlFileManager* t = TomlFileManager::getInstance();
	t->getValueIfExists("output.outputScheme_volume_preset", m_outputScheme_volume.preset);
	t->getValueIfExists("output.outputScheme_boundary_preset", m_outputScheme_boundary.preset);
	m_outputScheme_volume.updateBoolDataByPreset();
	m_outputScheme_boundary.updateBoolDataByPreset();
	m_initialized = true;
}

void FieldWriter::FieldOutputScheme::updateBoolDataByPreset() {
	rho = false;
	u = false;
	v = false;
	p = false;
	Cp = false;
	Ma = false;

	switch (preset) {
	case type_ruvp:
		rho = true;
		u = true;
		v = true;
		p = true;
		break;

	case type_ruvp_Cp:
		rho = true;
		u = true;
		v = true;
		p = true;
		Cp = true;
		break;

	case type_ruvp_Ma:
		rho = true;
		u = true;
		v = true;
		p = true;
		Ma = true;
		break;

	default:
		LogWriter::logAndPrintError("invalid preset @FieldWriter::FieldOutputScheme::updateBoolDataByPreset.\n");
		exit(-1);
		break;
	}
}

std::string FieldWriter::FieldOutputScheme::get_variable_names() {
	// 根据output_type，获取tecplot文件头的variables的内容
	std::vector<std::string> variables_vector_string;

	if (true)variables_vector_string.push_back("x");
	if (true)variables_vector_string.push_back("y");
	if (rho)variables_vector_string.push_back("rho");
	if (u)variables_vector_string.push_back("u");
	if (v)variables_vector_string.push_back("v");
	if (p)variables_vector_string.push_back("p");
	if (Cp)variables_vector_string.push_back("Cp");
	if (Ma)variables_vector_string.push_back("Ma");

	return stringVector_to_string_withQuote(variables_vector_string);
}

std::string FieldWriter::FieldOutputScheme::get_variable_value_string(myint nodeID, GPU::NodeSoA& node_host, GPU::OutputNodeFieldSoA& output_node_field) {
	std::vector<myfloat> float_vector;

	if (true)float_vector.push_back(node_host.xy[0][nodeID]);
	if (true)float_vector.push_back(node_host.xy[1][nodeID]);
	if (rho)float_vector.push_back(output_node_field.ruvp[0][nodeID]);
	if (u)float_vector.push_back(output_node_field.ruvp[1][nodeID]);
	if (v)float_vector.push_back(output_node_field.ruvp[2][nodeID]);
	if (p)float_vector.push_back(output_node_field.ruvp[3][nodeID]);
	if (Cp)float_vector.push_back(output_node_field.Cp[nodeID]);
	if (Ma)float_vector.push_back(output_node_field.Ma[nodeID]);

	// float vector转string，用空格分隔。不带换行符
	std::stringstream ss;
	ss << std::scientific;
	for (size_t i = 0; i < float_vector.size(); i++) {
		ss << float_vector[i];// 这部分代码占用了18.75%的时间
		if (i != float_vector.size() - 1) {
			ss << " ";
		}
	}
	return ss.str();
}

void FieldWriter::ReferenceData::update_farfield() {

	myfloat rho_far = GlobalPara::boundaryCondition::_2D::inf::ruvp[0];
	myfloat u_far = GlobalPara::boundaryCondition::_2D::inf::ruvp[1];
	myfloat v_far = GlobalPara::boundaryCondition::_2D::inf::ruvp[2];
	myfloat p_far = GlobalPara::boundaryCondition::_2D::inf::ruvp[3];
	myfloat V2_far = u_far * u_far + v_far * v_far;
	this->p_inf = p_far;
	this->p_dynamic_inf = 0.5 * rho_far * V2_far;
}

std::string FieldWriter::HistOutputScheme::get_variable_names() {
	std::vector<std::string> variables_vector_string;

	variables_vector_string.push_back("Iteration");
	variables_vector_string.push_back("physical_time");
	variables_vector_string.push_back("R_1");
	variables_vector_string.push_back("R_2");
	variables_vector_string.push_back("R_3");
	variables_vector_string.push_back("R_4");
	variables_vector_string.push_back("C_L");
	variables_vector_string.push_back("C_D");
	// "Iteration" "R_1" "R_2" "R_3" "R_4" "R_5" "R_6" "C_L" "C_D" "C_M_x" "C_M_y"

	return stringVector_to_string_withQuote(variables_vector_string);
}

std::string FieldWriter::HistOutputScheme::get_variable_value_string(HistData& histData) {
	std::vector<myfloat> float_vector;
	// iteration是int，单独处理
	float_vector.push_back(histData.physical_time);
	float_vector.push_back(histData.residual[0]);
	float_vector.push_back(histData.residual[1]);
	float_vector.push_back(histData.residual[2]);
	float_vector.push_back(histData.residual[3]);
	float_vector.push_back(histData.CL_wall);
	float_vector.push_back(histData.CD_wall);

	std::stringstream ss;
	ss << std::scientific;
	// iteration是int，单独处理
	ss << histData.iteration << " ";
	// float vector转string，用空格分隔。不带换行符
	for (size_t i = 0; i < float_vector.size(); i++) {
		ss << float_vector[i];
		if (i != float_vector.size() - 1) {
			ss << " ";
		}
	}
	return ss.str();
}

void FieldWriter::HistData::update(int _iteration, myfloat _physical_time, const myfloat _res[4], ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new) {
	std::vector<myfloat> _res_vector(4);
	for (int i = 0; i < 4; i++) {
		_res_vector[i] = _res[i];
	}
	update(_iteration, _physical_time, _res_vector, referenceData, boundary_host_new);
}

void FieldWriter::HistData::update(int _iteration, myfloat _physical_time, const std::vector<myfloat>& _res, ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new) {
	this->iteration = _iteration;
	this->physical_time = _physical_time;
	this->residual = _res;// 只要没用std::move，都是复制赋值 https://runebook.dev/zh/docs/cpp/container/vector/operator=
	set_CDCL_wall(referenceData, boundary_host_new);
}

void FieldWriter::HistData::set_CDCL_wall(ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new) {
	/*
	1 仅输出翼型气动力：搜索类型为wall的edgeSet，对气动力求和，作为histData的气动力
	2 输出所有边界气动力：需要存储边界名称
	目前先仅完成1，一步一步慢慢来
	*/
	const auto& edgeSets = boundary_host_new.edgeSets;
	this->CD_wall = 0.0;
	this->CL_wall = 0.0;
	if (referenceData.CD_vector.empty() || referenceData.CL_vector.empty()) {
		return;
	}
	for (size_t iSet = 0; iSet < edgeSets.size(); iSet++) {
		int setType = edgeSets[iSet].type;
		if (setType == _BC_wall_adiabat ||
			setType == _BC_wall_isothermal ||
			setType == _BC_wall_nonViscous ) {
			this->CD_wall += referenceData.CD_vector[iSet];
			this->CL_wall += referenceData.CL_vector[iSet];
		}
	}
}

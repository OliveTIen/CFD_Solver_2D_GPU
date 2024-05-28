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



// ������
inline std::string Quote(const std::string s) {
	return "\"" + s + "\"";
}

// ��stringVectorת��Ϊstring����", "���ӣ��Ҹ�Ԫ�ؼ���˫����
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
����Tecplot��װĿ¼ E:\Tecplot\Tecplot\Tecplot 360 EX 2023 R1\doc\360_data_format_guide.pdf

zonetype:
  fetriangle: ����Ԫ fe triangle
  fequadrilateral: �ı�Ԫ
  fetetrahedron: ������
  febrick: ������
datapacking:
  point: ���
  block: ����
*/

	// ���ļ�
	std::ofstream outfile;
	outfile.open(filePath);
	if (!outfile.is_open()) {
		LogWriter::logAndPrintError("failed to open file " + filePath +" @write_tecplot_volume_file\n");
		CExit::saveAndExit(-1);
	}

	std::string zone_title = "volume field";
	outfile 
	// ���⡢������
		<< "title= " << Quote(title) << "\n"
		<< "variables= " << m_outputScheme_volume.get_variable_names() << "\n"
	// ����
		<< "zone "
		<< "t= " << Quote(zone_title) << ", "
		<< "nodes= " << nodes.num_node << ", "
		<< "elements= " << elements.num_element << ", "
		<< "datapacking= " << "point" << ", "
		<< "zonetype= " << "fequadrilateral" << ", "
		<< "solutiontime= " << std::scientific << t_current << std::fixed << std::endl;

	// �ڵ�
	for (int i = 0; i < nodes.num_node; i++) {
		outfile << m_outputScheme_volume.get_variable_value_string(i, nodes, output_node_field) << "\n";
	}

	// ��Ԫ
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
	���������Ƕ�ά�ģ���˱߽���һά�ģ�Ҫ����һά��ͼ
	*/
	std::ofstream outfile;
	outfile.open(filePath);
	if (outfile.fail()) {
		LogWriter::logAndPrintError("failed to open file " + filePath + " @write_tecplot_boundary_file\n");
		return;
	}

	std::string title = "tecplot boundary file";
	outfile
		// ���⡢������
		<< "title= " << Quote(title) << "\n"
		<< "variables= " << m_outputScheme_boundary.get_variable_names() << "\n";
	
	const auto& edgeSets = boundary_host_new.edgeSets;
	//m_referenceData.CD_vector = std::vector<myfloat>(edgeSets.size(), 0.0);// ����edgeSet��������ֵ��ʼ��Ϊ0
	//m_referenceData.CL_vector = std::vector<myfloat>(edgeSets.size(), 0.0);
	const myfloat S_ref = GlobalPara::constant::referenceArea;// �ο����
	const myfloat one_on_Sref = 1 / S_ref;
	for (size_t iSet = 0; iSet < edgeSets.size(); iSet++) {
		const auto& edgeSet = edgeSets[iSet];// ��ǰ�߼���
		std::string zone_title = "boundary " + std::to_string(edgeSet.ID) + " " + edgeSet.name;// boundary 1 far
		myint boundary_edge_num = edgeSet.edge_vector.size();// ��ǰedge set��edge����
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
			�����ǰedgeSet�ı߽������
			ͨ����ͬһ��edgeSet����������edge�������ģ����ǰ�ߵ�node1���ں��ߵ�node0
			��Ҳ������������������������
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
			//�������������ɱ߽�ѹ�����ֵõ�
			//��ʽ��
			//dCD = - Cp / S_ref * dS * nx
			//dCL = - Cp / S_ref * dS * ny
			//Ȼ��Ա߽���л�·����
			//*/
			//
			//myint nx = edge_host.normal[0][edgeID];
			//myint ny = edge_host.normal[1][edgeID];
			//myint dS = edge_host.length[edgeID];// �߳���
			//myint Cp = 0.5 * (output_node_field.Cp[node0] + output_node_field.Cp[node1]);
			//myint dC_tmp = (-Cp) * one_on_Sref * dS;// �м������Ϊ�˼��ټ����������ã���������
			//m_referenceData.CD_vector[iSet] += dC_tmp * nx;
			//m_referenceData.CL_vector[iSet] += dC_tmp * ny;
		}
	}

	outfile.close();
}

void FieldWriter::ReferenceData::calculate_edgeSet_force(GPU::EdgeSoA& edge_host, GPU::OutputNodeFieldSoA& output_node_field, GPU::BoundaryV2& boundary_host_new) {
	/*
	�����߽�ߣ�����ѹ��ϵ�����ֵõ�������
	��ʽ��
	dCD = Cp / S_ref * dS * nx
	dCL = Cp / S_ref * dS * ny
	*/
	const auto& edgeSets = boundary_host_new.edgeSets;
	const myfloat S_ref = GlobalPara::constant::referenceArea;// �ο����
	const myfloat one_on_Sref = 1 / S_ref;

	// ��ʼ������size����ʼ��ֵΪ0
	this->CD_vector = std::vector<myfloat>(edgeSets.size(), 0.0);
	this->CL_vector = std::vector<myfloat>(edgeSets.size(), 0.0);

	for (size_t iSet = 0; iSet < edgeSets.size(); iSet++) {
		const auto& edgeSet = edgeSets[iSet];// ��ǰedgeSet
		myint boundary_edge_num = edgeSet.edge_vector.size();// ��ǰedgeSet��edge����
		for (myint iEdge = 0; iEdge < boundary_edge_num; iEdge++) {
			myint edgeID = edgeSet.edge_vector[iEdge];
			myint node0 = edge_host.nodes[0][edgeID];
			myint node1 = edge_host.nodes[1][edgeID];
			myfloat nx = edge_host.normal[0][edgeID];
			myfloat ny = edge_host.normal[1][edgeID];
			myfloat dS = edge_host.length[edgeID];// �߳���
			myfloat Cp = 0.5 * (output_node_field.Cp[node0] + output_node_field.Cp[node1]);
			myfloat dC_tmp = Cp * one_on_Sref * dS;// �м������Ϊ�˼��ټ����������ã���������
			this->CD_vector[iSet] += dC_tmp * nx;
			this->CL_vector[iSet] += dC_tmp * ny;
		}
	}
}

void FieldWriter::write_tecplot_hist_file(std::string filePath, int iteration, myfloat t_physics, myfloat residual[4], GPU::EdgeSoA& edge_host, GPU::OutputNodeFieldSoA& output_node_field, GPU::BoundaryV2& boundary_host_new) {
	/*
	�ж��Ƿ�Ӧ��д�ļ�ͷ
	�ڵ�һ�ε��øú���������£�
	��hist�ļ������ڣ�����Ҫд�ļ�ͷ
	��GlobalPara::basic::_continue==false������Ҫд�ļ�ͷ
	*/
	bool b_should_write_head = false;
	if (!m_called_write_tecplot_hist_file) {
		std::ifstream infile(filePath, std::ios_base::in);// ֻ��
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
			// ���⡢������
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
	// ����ݴ��ļ��������´�����
	// Ŀǰ������ֱ����GPUID

    // ��������
	std::ofstream outfile;
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	BoundaryManager_2D& boundaryManager = pFVM2D->boundaryManager;

	// ���ļ�
	outfile.open(filePath);
	if (!outfile.is_open())std::cout << "Error: Fail to open file " << filePath << std::endl;
	// ���ʱ�䲽��
	outfile << "t, istep, CFL" << std::endl;
	outfile << t_current << " " << i_step << " "<< GlobalPara::time::CFL << std::endl;// CFL���Զ���CFL��������Ҫ��ȡ�ģ�û�иò��Զ�ȡ��Ҳû��
	// ����ڵ�
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
		//һ��̫���ᵼ�¶�ȡʱbuffer���Ȳ������Ӷ���������Ԥ�ϵĽ��(��ѭ������ȡʧ�ܵ�)
		int nEdge = (int)boundaryManager.boundaries[ib].pEdges.size();//�������' '����'\n'
		bool check = 0;//����
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
		if (!check)outfile << std::endl;//���ϴ������' '
	}
	outfile.close();

}

void FieldWriter::allocNodeFieldDataUsingOutputScheme(GPU::OutputNodeFieldSoA& nodeField, myint num_node) {
	/*
	����������������ڴ档��Ҫ�ڶ�ȡ��Config��Mesh����ã���ΪҪ֪����������ͽڵ����
	*/
	// �ж��Ƿ��Ѷ�ȡConfig
	if (!m_initialized) {
		LogWriter::logAndPrintError("outputScheme has not been initialized. @FieldWriter::allocNodeFieldDataUsingOutputScheme\n");
		exit(-1);
	}
	// �ж��Ƿ��ȡ����Mesh����num_node==0��ʾδ��ȡ����
	if (num_node == 0) {
		LogWriter::logAndPrintError("num_node == 0. @FieldWriter::allocNodeFieldDataUsingOutputScheme\n");
		exit(-1);
	}
	//�ж��Ƿ��ظ������ڴ�
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
	����������element_vruvpת��Ϊ�������nodeField��ruvp

	���ܵ�������
	���ǵ�Ԫ���ĸ���洢��һ����nodeID[3]=-1������nodeID[3]==nodeID[2]���ᵼ��nodeID[2]Ȩ��ƫ��
	*/
	auto& node_ruvp = nodeField.ruvp;
	myint num_node = nodeField.num_node;
	myint num_element = element.num_element;
	const int nNodePerElement = 4;
	const int nValuePerNode = 4;
	int* node_neighborElement_num;// ��¼ÿ���ڵ���ھӵ�Ԫ����
	// ������Դ ��ʼ��
	node_neighborElement_num = new int[num_node]();// С�����Զ���ʼ��Ϊ0
	for (int i = 0; i < 4; i++) {
		memset(node_ruvp[i], 0, num_node * sizeof(myfloat));
	}

	// ����Ԫֵ�ӵ����ӽڵ�ֵ��
	for (int iElement = 0; iElement < num_element; iElement++) {
		// ͳ��ÿ���ڵ���ھӵ�Ԫ���������޸Ľڵ�ֵ
		for (int jNode = 0; jNode < nNodePerElement; jNode++) {
			// ��ȡ��Ԫ�ڵ�ID
			int GPUID_of_node = element.nodes[jNode][iElement];// GPUID of node 0
			// nodeID[3]=-1ʱ�������ýڵ�
			if (GPUID_of_node < 0 || GPUID_of_node >= num_node)continue;// ��������ѭ������������ѭ����
			// ��ID��Ӧ���ھӵ�Ԫ����+1
			node_neighborElement_num[GPUID_of_node]++;
			// ��ID��Ӧ�Ľڵ�����ֵ���ϵ�Ԫֵ
			for (int kValue = 0; kValue < nValuePerNode; kValue++) {
				node_ruvp[kValue][GPUID_of_node] += element_vruvp[kValue][iElement];
			}
		}
	}

	// �ڵ�ֵ�����ھӵ�Ԫ�����õ�ƽ��ֵ����Ϊ�ڵ�ruvpֵ
	for (int iNode = 0; iNode < num_node; iNode++) {
		// Ϊ�˱������0����ĸ����0������
		if (node_neighborElement_num[iNode] == 0)continue;
		// node_ruvp�����ھӵ�Ԫ�����õ�ƽ��ֵ
		for (int kValue = 0; kValue < nValuePerNode; kValue++) {
			node_ruvp[kValue][iNode] /= node_neighborElement_num[iNode];
		}
	}

	// �ͷ���Դ
	delete[] node_neighborElement_num;
}

void FieldWriter::update_nodeField_other_variables(GPU::OutputNodeFieldSoA& nodeField) {
	// ����Ƿ������ڴ�
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
			// ����ϵ�� Cp = (p - p_inf)/p_dynamic_inf
			m_referenceData.update_farfield();
			nodeField.Cp[i] = (p - m_referenceData.p_inf) / m_referenceData.p_dynamic_inf;
			// ������ ��write_tecplot_boundary_file
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
	// ��ʼ�����������ֻ��readConfigʱ����
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
	// ����output_type����ȡtecplot�ļ�ͷ��variables������
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

	// float vectorתstring���ÿո�ָ����������з�
	std::stringstream ss;
	ss << std::scientific;
	for (size_t i = 0; i < float_vector.size(); i++) {
		ss << float_vector[i];// �ⲿ�ִ���ռ����18.75%��ʱ��
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
	// iteration��int����������
	float_vector.push_back(histData.physical_time);
	float_vector.push_back(histData.residual[0]);
	float_vector.push_back(histData.residual[1]);
	float_vector.push_back(histData.residual[2]);
	float_vector.push_back(histData.residual[3]);
	float_vector.push_back(histData.CL_wall);
	float_vector.push_back(histData.CD_wall);

	std::stringstream ss;
	ss << std::scientific;
	// iteration��int����������
	ss << histData.iteration << " ";
	// float vectorתstring���ÿո�ָ����������з�
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
	this->residual = _res;// ֻҪû��std::move�����Ǹ��Ƹ�ֵ https://runebook.dev/zh/docs/cpp/container/vector/operator=
	set_CDCL_wall(referenceData, boundary_host_new);
}

void FieldWriter::HistData::set_CDCL_wall(ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new) {
	/*
	1 �������������������������Ϊwall��edgeSet������������ͣ���ΪhistData��������
	2 ������б߽�����������Ҫ�洢�߽�����
	Ŀǰ�Ƚ����1��һ��һ��������
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

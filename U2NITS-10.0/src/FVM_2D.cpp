#include "FVM_2D.h"
#include "output/ConsolePrinter.h"
#include "output/ResidualCalculator.h"
#include "output/LogWriter.h"
#include "global/StringProcessor.h"
#include "global/SystemInfo.h"

FVM_2D* FVM_2D::pFVM2D = nullptr;

FVM_2D::FVM_2D() {
	pFVM2D = this;
}

void FVM_2D::run() {
	//�����񣬳�ʼ����ֵ�����

	//��ȡ�����ļ����������ļ�����ʼ�������ڲ���ֵ
	std::cout << "\nPreProcessing..." << std::endl;
	//std::cout << "\033[36m" << "\nPreProcessing..." << std::endl << "\033[0m";
	LogWriter::writeLog("PreProcessing\n", 1);
	//LogWriter logWritter();
	if (GlobalPara::basic::_continue) {
		//��ȡ�����ļ�
		int ret = readContinueFile();
		//��ȡ�����ļ��������ڲ���ֵ
		if (ret == -1) {
			LogWriter::writeLogAndCout("Warning: Fail to read previous mesh(pause_*.dat)\n");
			int ret2 = readMeshFile(".inp");
			if (ret2 == -1) {
				LogWriter::writeLogAndCout("Program will exit.\n");
				return;
			}
			setInitialCondition();
		}		
	}
	else {
		int ret2 = readMeshFile(".inp");
		if (ret2 == -1) {
			LogWriter::writeLogAndCout("Program will exit.\n");
			return;
		}
		setInitialCondition();
	}
	
	//���߽�edge���ϱ�ǩ����������ڱ߽�������
	{
		int ret = boundaryManager.iniBoundaryEdgeSetID_and_iniBoundaryType(this);
		if (ret != 0)return;
	}

	//��־��¼�߽����
	logBoundaryCondition();

	//�����
	std::cout << "\nSolving..." << std::endl;
	//std::cout << "\033[36m" << "\nSolving..." << std::endl << "\033[0m";
	if (GlobalPara::physicsModel::equation == _EQ_euler) {
		solve_("_Euler.out", "_Euler.info");
	}

}

void FVM_2D::generateMeshScript_gmsh(std::string suffix) {
	//std::string str = "gmsh 2D.geo -2 -order 1 -format inp";
	std::string m = "gmsh";
	std::string op_name = " ./input/" + GlobalStatic::filename + ".geo";//option_name
	std::string op_dimension = " -2";
	std::string op_elementOrder = " -order 1";
	std::string op_outputformat = " -format " + suffix;
	std::string str = m + op_name + op_dimension + op_elementOrder + op_outputformat;
	system(str.c_str());
}

int FVM_2D::readMeshFile(std::string suffix) {
	LogWriter::writeLogAndCout("read mesh file\n");
	std::string dir = GlobalStatic::getExePath_withSlash(1) + "input\\";
	std::ifstream infile(dir + GlobalStatic::filename + suffix);
	if (!infile) {
		LogWriter::writeLogAndCout("Warning: fail to read mesh(*.inp). (FVM_2D::readMeshFile) \n");
		return -1;
	}
	int state = 0;//1-Node, 2-Element
	int substate = 0;//1-Line
	int maxnodeID = 1;
	int maxedgeID = 1;
	int maxelementID = 1;
	const int bufferLength = 300;
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<VirtualEdge_2D> vBoundaryEdges;//��ʱ�������洢�߽�Ԫ�ڵ�ID
	vBoundaryEdges.push_back(VirtualEdge_2D());//���1��nullԪ�أ��Ա�֤�кű�ʾID
	std::vector<int>vInts;//��ʱ����
	std::string setName;//��ʱ����
	while (1) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())break;
		tLine = buffer;
		tWords = GlobalStatic::splitString(tLine);
		if (tWords.size() == 0)continue;//���ڿ��У�ǿ�Ƚ�����һ��ѭ������ֹ��ȡtWords[0]�����ڴ����
		else {
			if (tWords[0] == "*NODE")state = 1;
			if (tWords[0] == "*ELEMENT") {
				state = 2;
				if (tWords[2].substr(0, 7) == "ELSET=L") substate = 1;//ELSET=Line[x]
				else if (tWords[2].substr(0, 7) == "ELSET=S") substate = 2;
			}
			if (tWords[0] == "*ELSET"){
				state = 3;
			}
			if (tWords[0] == "*NSET") state = 4;//Ŀǰû�õ������gmsh����inpʱ���Բ������ڵ�
		}

		//translate//ע������stoi
		if (state == 1 && tWords[0].substr(0,1)!= "*") {//"*NODE" ��ȡ�ڵ�ID������
			Node_2D node;
			node.ID= (int)std::stod(tWords[0]);
			node.x= std::stod(tWords[1]);
			node.y= std::stod(tWords[2]);
			nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 2 && substate == 1 && tWords[0].substr(0, 1) != "*") {//�߽�Ԫ "*ELEMENT""ELSET=Line"
			VirtualEdge_2D ve;//��Ա��nodes[2]
			ve.nodes[0]= std::stoi(tWords[1]);//1, 1, 6; EdgeID,node0,node1
			ve.nodes[1]= std::stoi(tWords[2]);
			vBoundaryEdges.push_back(ve);
		}
		else if (state == 2 && substate == 2 && tWords[0].substr(0, 1) != "*") {//�浥Ԫ "*ELEMENT""ELSET=Surface"
			Element_T3 e;
			e.ID = (int)std::stod(tWords[0]);
			e.nodes[0] = (int)std::stod(tWords[1]);
			e.nodes[1] = (int)std::stod(tWords[2]);
			e.nodes[2] = (int)std::stod(tWords[3]);
			elements.push_back(e);
			maxelementID = (std::max)(maxelementID, e.ID);
		}
		else if (state == 3) {//"*ELSET"
			if (tWords[0].substr(0, 1) == "*") {
				//��ʼ��vBoundarySets
				std::vector<int> start_end = boundaryManager.splitInts(vInts);
				int nSets = (int)start_end.size() / 2;
				for (int is = 0; is < nSets; is++) {//����nSets!=0��˵��setName�Ѿ�����ʼ��
					VirtualBoundarySet_2D vb;
					vb.name = setName;
					vb.ID = (int)boundaryManager.vBoundarySets.size() + 1;
					vb.startID = start_end[is * 2 + 0];//0,2,4,... 0=is*2+0
					vb.endID = start_end[is * 2 + 1];
					boundaryManager.vBoundarySets.push_back(vb);
				}
				//����
				vInts.clear();
				setName = tWords[1].substr(6);//"wall""outlet""inlet"...
			}
			else {
				for (int iword = 0; iword < tWords.size(); iword++) {
					vInts.push_back(std::stoi(tWords[iword]));
				}
			}
		}
	}
	infile.close();

	std::cout << "Generate and initiate elements..." << std::endl;

	iniPNodeTable(maxnodeID);
	iniEdges();
	iniPEdgeTable();
	iniPElementTable(maxelementID);
	iniElement_xy_pEdges();
	iniNode_neighborElements();
	boundaryManager.iniBoundarySetPEdges_in_readMeshFile(this, vBoundaryEdges);//��ʼ��vBoundarySet��pEdges


	return 0;
}

void FVM_2D::setInitialCondition() {
	switch (GlobalPara::initialCondition::type) {
	case 1:
	{
		//1.����
		using namespace GlobalPara::boundaryCondition::_2D;
		for (int ie = 0; ie < elements.size(); ie++) {
			Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, Constant::gamma);
			//if(elements[ie].ID==173)system()
		}
		std::string str;
		str += "InitialCondition:\nUniform flow, ruvp\t" + StringProcessor::DoubleArray_2_String(inf::ruvp, 4) + "\n";
		LogWriter::writeLogAndCout(str);
	}
		break;

	case 2:
	{
		//2.�������͵����еĵ���

		using namespace GlobalPara::boundaryCondition::_2D;
		std::string str;
		str += "InitialCondition:\nUniform flow + isentropicVortex, ruvp of Uniform flow:" + StringProcessor::DoubleArray_2_String(inf::ruvp, 4) + "\n";
		LogWriter::writeLogAndCout(str);
		isentropicVortex_2(5, 5, 5, inf::ruvp);
		//for (int ie = 0; ie < elements.size(); ie++) {
		//	//������+������ Isentropic vortex
		//	double deltau, deltav, deltaT;
		//	double x = elements[ie].x; double y = elements[ie].y;
		//	double xc = 5; double yc = 5;
		//	isentropicVortex(x, y, xc, yc, 1, deltau, deltav, deltaT);
		//	double ruvp[4]{};
		//	ruvp[1] = inf::ruvp[1] + deltau;
		//	ruvp[2] = inf::ruvp[2] + deltav;
		//	double rho0 = inf::ruvp[0];
		//	double p0 = inf::ruvp[3];
		//	double T0 = p0 / rho0 / Constant::gamma;//p=rho*R*T
		//	double deltap = deltaT * rho0 / (1 - Constant::R * rho0 * T0 / pow(p0, Constant::gamma));
		//	double deltarho = deltap * rho0 / pow(p0, Constant::gamma);
		//	ruvp[0] = rho0 + deltarho;
		//	ruvp[3] = p0 + deltap;
		//	Math_2D::ruvp_2_U(ruvp, elements[ie].U, Constant::gamma);
		//}

	}
		break;

	}

}

void FVM_2D::logBoundaryCondition() {
	//��־��¼�߽����
	std::string str;
	str += "BoundaryCondition:\n";
	str += "inlet::ruvp\t" + StringProcessor::DoubleArray_2_String(GlobalPara::boundaryCondition::_2D::inlet::ruvp, 4)
		+ "\noutlet::ruvp\t" + StringProcessor::DoubleArray_2_String(GlobalPara::boundaryCondition::_2D::outlet::ruvp, 4)
		+ "\ninf::ruvp\t" + StringProcessor::DoubleArray_2_String(GlobalPara::boundaryCondition::_2D::inf::ruvp, 4)
		+ "\n";
	LogWriter::writeLog(str);
}

void FVM_2D::solve_(std::string suffix_out, std::string suffix_info) {
	//�����ʽ(tecplot)�����ļ����������ÿһ֡���һ���ļ�
	std::vector<Element_T3> elements_old;

	//���Ʊ���
	double t = t_solve;//��ȡcontinue�ļ�������ʱ�䡣���ڿ����Ƿ���ֹ����
	const double T = GlobalPara::time::T;//������ʱ�䡣���ڿ����Ƿ���ֹ����
	int nFiles = 0;//����ļ�������������ʾ
	bool signal_pause = 0;//��ͣ�źţ����ڿ���
	//��ʱ���������ڼ������ʱ��(����CPU����)
	clock_t start_t;//λ������������У����ڼ���CPUʱ��
	start_t = clock();
	//�������������ļ�ͷ
	if (GlobalPara::initialCondition::type == 2)writeFileHeader_isentropicVortex();
	//���ƽ�����
	ConsolePrinter consolePrinter;
	consolePrinter.restoreCursorPosition();
	std::cout << "Progress:";	ConsolePrinter::drawProgressBar(int(t / T));
	for (int istep = istep_solve; istep < 1e6 && t < T; istep++) {
		elements_old = elements;
		//����ʱ��
		caldt();
		dt = (t + dt > T) ? (T - t) : dt;//if (t + dt > T)dt = T - t;
		t += dt;
		//����ͨ����ʱ���ƽ�
		solver.evolve(dt);

		//������
		if (istep % 10 == 1) {
			consolePrinter.EraseToLastCursorPosition();
			std::cout << "Progress:";	ConsolePrinter::drawProgressBar(int(t / T * 100 + 2));//+x��Ϊ�˷�ֹͣ��99%
			//���ʱ��
			double time_used = (double)(clock() - start_t) / CLOCKS_PER_SEC;//
			//���ʱ��
			std::cout << "\n";
			std::cout << "\ttime used: " <<math.timeFormat((int)time_used)<<" ( " << (int)time_used << " s )\n";
			std::cout << "\tsolution time: " << t << " s / " << T << " s\n";
			std::cout << "\tsolution step: " << istep << "\t";
			std::cout << "\tnumber of output file: " << nFiles << "\n";
			std::cout << "Press ESC to end Computation\n";
		}
		//���Ƿ�ֵ
		if (isNan()) {
			std::cout << "\nWarning: \"NaN\" detected, look into LOG for details. Computation stopped.\n";
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
			writeContinueFile(GlobalStatic::getExePath_withSlash(1) + "output\\nan_" + GlobalStatic::filename + "[" + szBuffer + "].dat", t, istep);
			break;
		}
		//������
		if (istep % GlobalPara::output::step_per_output == 1) {
			//.dat ������ʾ�ļ� tecplot��ʽ
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//��istepС��4λ������0
			writeTecplotFile(GlobalStatic::getExePath_withSlash(1) + "output\\" + GlobalStatic::filename + "[" + szBuffer + "].dat", t);
			nFiles++;
			//.autosave �Զ������ļ� �����ļ�
			updateAutoSaveFile(t, istep);
			//�����е�����ļ� 
			if (GlobalPara::initialCondition::type == 2) {
				//���ʱ��
				//end_t = clock();
				double time_used = (double)(clock() - start_t) / CLOCKS_PER_SEC;
				//�������ļ�
				ResidualCalculator::cal_error_isentropicVortex(0, 0, 10, 10, 5, t, istep, time_used, GlobalPara::boundaryCondition::_2D::inf::ruvp);
			}
		}
		//��esc��ֹ����
		if (_kbhit()) {
			char ch = _getch();
			if (ch == 27) {
				char szBuffer[20];
				sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
				writeContinueFile(GlobalStatic::getExePath_withSlash(1) + "output\\pause_" + GlobalStatic::filename + "[" + szBuffer + "].dat",t,istep);
				std::cout << "[ESC]: Computation Interrupted\n";
				break;
			}
		}
		//��ֹ
		if (t >= T) {
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
			writeContinueFile(GlobalStatic::getExePath_withSlash(1) + "output\\pause_" + GlobalStatic::filename + "[" + szBuffer + "].dat", t, istep);
			std::cout << "Computation finished\n";
			break;
		}
		else if (isStable(elements_old)) {
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
			writeContinueFile(GlobalStatic::getExePath_withSlash(1) + "output\\pause_" + GlobalStatic::filename + "[" + szBuffer + "].dat", t, istep);
			std::cout << "Computation finished as the field is already stable\n";
			break;
		}
	}

}

void FVM_2D::caldt() {
	//����Ϊ�´���
	for (int ie = 0; ie < elements.size(); ie++) {
		double LambdaC_i = elements[ie].calLambdaFlux(this);
		double Omega_i = elements[ie].calArea(this);
		double dt_i = GlobalPara::time::CFL * Omega_i / LambdaC_i;
		dt = (std::min)(dt, dt_i);
	}
}

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

void FVM_2D::writeContinueFile(std::string f_name, double t, int istep) {
	//std::cout << "wirtefile:" << f_name << std::endl;
	int flag_OUT = -1;
	std::ofstream f(f_name);
	if (!f.is_open())std::cout << "Error: Fail to open file " << f_name << std::endl;

	f << "t, istep" << std::endl;
	f << t << " " << istep << std::endl;

	f << "nodes: ID, x, y" << std::endl;
	for (int in = 0; in < nodes.size(); in++) {
		std::vector<double> uu_i = nodes[in].calNodeValue(this);
		f << nodes[in].ID << " ";
		f << nodes[in].x << " ";
		f << nodes[in].y << std::endl;
	}

	f << "elements: ID, nodes, U" << std::endl;
	for (int ie = 0; ie < elements.size(); ie++) {
		f << elements[ie].ID << " ";
		f << elements[ie].nodes[0] << " ";
		f << elements[ie].nodes[1] << " ";
		f << elements[ie].nodes[2] << " ";
		f << elements[ie].U[0] << " ";
		f << elements[ie].U[1] << " ";
		f << elements[ie].U[2] << " ";
		f << elements[ie].U[3] << std::endl;
	}

	f << "boundaries: ID, name; edgeIDs" << std::endl;
	for (int ib = 0; ib < boundaryManager.vBoundarySets.size(); ib++) {
		f << boundaryManager.vBoundarySets[ib].ID << " ";
		f << boundaryManager.vBoundarySets[ib].name << std::endl;
		//һ��̫���ᵼ�¶�ȡʱbuffer���Ȳ������Ӷ���������Ԥ�ϵĽ��(��ѭ������ȡʧ�ܵ�)
		int nEdge = (int)boundaryManager.vBoundarySets[ib].pEdges.size();//�������' '����'\n'
		bool check = 0;//����
		for (int iEdge = 0; iEdge < nEdge;iEdge++) {
			f << boundaryManager.vBoundarySets[ib].pEdges[iEdge]->ID;
			if (iEdge % 10 == 9){
				f << std::endl; 
				check = 1;
			}
			else {
				f << ' ';
				check = 0;
			}
		}
		if(!check)f << std::endl;//���ϴ������' '
	}
	f.close();

}

int FVM_2D::readContinueFile() {
	std::string m_path = GlobalStatic::getExePath_withSlash(1) + "output\\";
	std::vector<std::string> files = GlobalStatic::ls(m_path);//filesΪ�ַ����������洢��·���������ļ���
	int index_maxNstep = -1;
	//��Ѱ�Ƿ���pause_filename[xxx].dat�ļ������У�����xxx����
	int maxNstep = 0;
	int tmp_index = 0;//����ָʾ"[", "]"
	std::string tmp_str;
	for (int i = 0; i < files.size(); i++) {
		//��Ѱ��pause_filename��ͷ���ļ���
		std::string str_match = "pause_" + GlobalStatic::filename;
		if (files[i].substr(0, int(str_match.size())) == str_match) {
			//�޳�"pause_filename["
			int pre_length = int(str_match.size() + 1);//"pause_filename["�ĳ���
			std::string str = files[i].substr(pre_length);
			//�޳�"].dat"
			int post_length = 5;//"].dat"�ĳ���
			for (int i = 0; i < post_length; i++) {
				str.pop_back();
			}
			//��xxx����num
			int num = std::stoi(str);
			if (num > maxNstep) {
				index_maxNstep = i;
				maxNstep = num;
			}
		}
	}
	//���ޣ��򷵻�-1
	if (index_maxNstep == -1)return -1;
	std::ifstream infile(m_path + files[index_maxNstep]);
	if (!infile) {
		return -1;
	}

	std::cout << "Read previous mesh:"<< files[index_maxNstep] << std::endl;

	int state = -1;//1-Node, 2-Element
	int maxnodeID = 1;
	int maxedgeID = 1;
	int maxelementID = 1;
	int iline_set = 0;//��ʱ���������ڼ�¼��ȡsetʱ�����˶�����
	std::vector<std::vector<int>> edges_of_all_sets;//��ʱ�������洢��set��edgeID
	const int bufferLength = 600;//!!!С�ĳ��Ȳ���
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<int> edge_ID_long;
	int nNullRow = 0;
	bool loop = 1;
	while (loop) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())loop=0;
		tLine = buffer;
		tWords = GlobalStatic::splitString(tLine);
		if (nNullRow >= 10)loop = 0;
		if (tWords.size() == 0) {
			nNullRow++;
			continue;
		}//���У�ǿ�ȿ�ʼ��һ��ѭ������ֹtWords[0]�ڴ����
		else {
			if (tWords[0] == "t")state = 0;
			if (tWords[0] == "nodes:")state = 1;
			if (tWords[0] == "edges:")state = 2;
			if (tWords[0] == "elements:")state = 3;
			if (tWords[0] == "boundaries:")state = 4;
		}

		//translate
		if (state == 0 && tWords[0].substr(0, 1) != "t") {//t_solve, istep_solve
			t_solve = std::stod(tWords[0]);
			istep_solve = (int)std::stod(tWords[1]);
		}
		else if (state == 1 && tWords[0].substr(0, 1) != "n") {//nodes: ID, x, y
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 3 && tWords[0].substr(0, 1) != "e") {//elements: ID, nodes, U
			Element_T3 ele;
			ele.ID = (int)std::stod(tWords[0]);
			ele.nodes[0] = (int)std::stod(tWords[1]);
			ele.nodes[1] = (int)std::stod(tWords[2]);
			ele.nodes[2] = (int)std::stod(tWords[3]);
			ele.U[0] = std::stod(tWords[4]);
			ele.U[1] = std::stod(tWords[5]);
			ele.U[2] = std::stod(tWords[6]);
			ele.U[3] = std::stod(tWords[7]);
			elements.push_back(ele);
			maxelementID = (std::max)(maxelementID, ele.ID);
		}
		else if (state == 4 && tWords[0].substr(0, 1) != "b") {//boundaries: ID, name, edgeIDs
			iline_set++;
			if (!isdigit(tWords[1][0])) {//���������� 1.vBoundarySets������Ŀ������Ŀ�Ĳ��ֳ�ʼ�� 2.edges_of_all_sets������Ŀ
				//set��ID,name��ʼ��
				VirtualBoundarySet_2D vb;
				vb.ID = (int)std::stod(tWords[0]);
				vb.name = tWords[1];
				boundaryManager.vBoundarySets.push_back(vb);
				//edges_of_all_sets������Ŀ
				edges_of_all_sets.push_back(std::vector<int>());
			}
			else {//�������� 1.edges_of_all_sets�ĳ�ʼ��
				//�������һ��BoundarySet
				std::vector<int>& edges_of_current_set = edges_of_all_sets[edges_of_all_sets.size() - 1];
				std::vector<int> intVector = StringProcessor::Words2Ints(tWords);
				boundaryManager.attachToVector(edges_of_current_set, intVector);
			}
		}
	}

	infile.close();

	std::cout << "Generate and initiate elements..." << std::endl;

	iniPNodeTable(maxnodeID);
	iniEdges();
	iniPEdgeTable();
	iniPElementTable(maxelementID);
	iniElement_xy_pEdges();
	iniNode_neighborElements();
	boundaryManager.iniBoundarySetPEdges_in_readContinueFile(this, edges_of_all_sets);

	return 0;
}

void FVM_2D::updateAutoSaveFile(double t, int istep) {
	//��ౣ��global::autosaveFileNum��autosave�ļ����ٶ���д֮ǰ��
	writeContinueFile(GlobalStatic::getExePath_withSlash(1) + "output\\autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalStatic::filename + ".dat", t, istep);
	iCurrentAutosaveFile++;
	if (iCurrentAutosaveFile >= GlobalStatic::autosaveFileNum)iCurrentAutosaveFile = 0;
}

void FVM_2D::writeBinaryFile(std::string f_name) {

}

void FVM_2D::calculateNodeValue() {

	if (GlobalPara::output::output_var_ruvp[0]) 		rho_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[1]) 		u_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[2]) 		v_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[3]) 		p_nodes.resize(nodes.size());

	for (int i = 0; i < nodes.size(); i++) {
		std::vector<double> U_node = nodes[i].calNodeValue(this);
		double U[4]{ U_node[0],U_node[1],U_node[2],U_node[3] };
		double ruvp[4];
		math.U_2_ruvp(U, ruvp, Constant::gamma);
		if (GlobalPara::output::output_var_ruvp[0]) 		rho_nodes[i] = ruvp[0];
		if (GlobalPara::output::output_var_ruvp[1]) 		u_nodes[i] = ruvp[1];
		if (GlobalPara::output::output_var_ruvp[2]) 		v_nodes[i] = ruvp[2];
		if (GlobalPara::output::output_var_ruvp[3]) 		p_nodes[i] = ruvp[3];
	}

}

void FVM_2D::iniPNodeTable(int maxnodeID) {
	pNodeTable.resize(maxnodeID + 1);
	for (int in = 0; in < nodes.size(); in++) {
		pNodeTable[nodes[in].ID] = &(nodes[in]);
	}
}

void FVM_2D::iniEdges() {
	for (int ie = 0; ie < elements.size(); ie++) {
		int n0, n1, n2;
		n0 = elements[ie].nodes[0];
		n1 = elements[ie].nodes[1];
		n2 = elements[ie].nodes[2];
		iniEdges_registerSingle(n0, n1, &(elements[ie]));
		iniEdges_registerSingle(n1, n2, &(elements[ie]));
		iniEdges_registerSingle(n2, n0, &(elements[ie]));
	}
}

void FVM_2D::iniPEdgeTable() {
	pEdgeTable.resize(edges.size() + 1);
	for (int i = 0; i < edges.size(); i++) {
		pEdgeTable[edges[i].ID] = &(edges[i]);
	}
}

void FVM_2D::iniNode_neighborElements() {
	for (int ie = 0; ie < elements.size(); ie++) {
		elements[ie].nodes;
		for (int jn = 0; jn < 3; jn++) {
			Node_2D* pN = getNodeByID(elements[ie].nodes[jn]);
			pN->neighborElements.push_back(&(elements[ie]));
		}
	}
}

Edge_2D* FVM_2D::isEdgeExisted(int n0, int n1) {
	//n0, n1��ʾ�ڵ�ID
	if(edges.size()==0)return nullptr;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].nodes[0] == n0 && edges[i].nodes[1] == n1)return &(edges[i]);
		if (edges[i].nodes[0] == n1 && edges[i].nodes[1] == n0)return &(edges[i]);
	}
	return nullptr;
}

void FVM_2D::iniEdges_registerSingle(int n0, int n1, Element_T3* pE) {
	Edge_2D* pEdge = isEdgeExisted(n0, n1);
	if (pEdge == nullptr) {
		Edge_2D tmp_edge;
		tmp_edge.ID = (int)edges.size() + 1;
		tmp_edge.nodes[0] = n0;
		tmp_edge.nodes[1] = n1;
		tmp_edge.pElement_L = pE;
		edges.push_back(tmp_edge);
	}
	else {
		pEdge->pElement_R = pE;
	}
}

void FVM_2D::iniPElementTable(int maxelementID) {
	pElement_T3Table.resize(maxelementID + 1);
	for (int ie = 0; ie < elements.size(); ie++) {
		pElement_T3Table[elements[ie].ID] = &(elements[ie]);
	}
}

void FVM_2D::iniElement_xy_pEdges() {
	for (int ie = 0; ie < elements.size(); ie++) {
		elements[ie].inixy(this); //�Ѿ�calxy���Ժ󲻱ص��ġ����Ǳ���Ҫ�ȶ�ȡnode�ٶ�ȡelement
		elements[ie].iniPEdges(this);
	}
}

Node_2D* FVM_2D::getNodeByID(int ID) {
	//if (pNodeTable.size() == 0) {
	//	std::cout << "Error: uninitialized pNodeTable. Initializing..." << std::endl;
	//	iniPNodeTable();
	//	//return nullptr;
	//}
	return pNodeTable[ID];
}

void FVM_2D::isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT) {
	double xbar = x - xc;
	double ybar = y - yc;
	double r2 = xbar * xbar + ybar * ybar;
	const double gamma = Constant::gamma;
	const double PI = Constant::PI;
	deltau = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * (-ybar);
	deltav = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * xbar;
	deltaT = -(gamma - 1) / chi / chi * 8 * gamma * PI * PI * exp(1 - r2);
}

void FVM_2D::isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0) {
	//ruvp0������������
	for (int ie = 0; ie < elements.size(); ie++) {
		Element_T3& e = elements[ie];
		double rho, u, v, p;
		double xbar, ybar, r2, du, dv, dT;
		const double PI = Constant::PI;
		const double gamma = Constant::gamma;
		xbar = e.x - xc;
		ybar = e.y - yc;
		r2 = xbar * xbar + ybar * ybar;
		du = chi / 2. / PI * exp(0.5 * (1. - r2)) * (-ybar);
		dv = chi / 2. / PI * exp(0.5 * (1. - r2)) * xbar;
		u = ruvp0[1] + du;
		v = ruvp0[2] + dv;
		dT = -(gamma - 1.) * chi * chi / (8. * gamma * PI * PI) * exp(1. - r2);
		rho = pow(ruvp0[3] + dT, 1. / (gamma - 1.));
		p = rho * (ruvp0[3] + dT);
		double ruvp[4]{ rho,u,v,p };
		Math_2D::ruvp_2_U(ruvp, e.U, gamma);
	}
	//lambda���ʽ [ capture ] ( params ) opt -> ret { body; }; http://c.biancheng.net/view/3741.html
	//���� auto f = [](int a) -> int { return a + 1; };
	//auto fun_WA = [](const double* xy) {
}


void FVM_2D::writeFileHeader_isentropicVortex() {
	std::ofstream outfile(GlobalStatic::exePath_withSlash + "output\\" + "error_isentropicVortex_" + GlobalStatic::filename + ".txt", std::ios::app);//׷��ģʽ
	outfile << SystemInfo::getCurrentDateTime() << "\n";
	outfile
		<< "istep" << "\t"
		<< "cpu_time[s]" << "\t"
		<< "error_L1" << "\t"
		<< "error_L2" << "\t"
		<< "error_max"<< "\n";
	outfile.close();
}

bool FVM_2D::isStable(std::vector<Element_T3> old) {
	double sum = 0;
	double res;
	for (int ie = 0; ie < elements.size(); ie++) {
		//sum += abs(elements[ie].calculateu() - old[ie].calculateu());
		res = elements[ie].U[1] - old[ie].U[1];
		sum += abs(res);
	}
	//std::cout << sum << std::endl;
	if (sum <= 1e-12)return 1;
	else return 0;
}

bool FVM_2D::isNan() {
	bool is_nan = 0;
	for (int ie = 0; ie < elements.size(); ie++) {
		//ֻ��Ҫ���1��� rho
		if (isnan(elements[ie].U[0])) {
			is_nan = 1;
			const double& x = elements[ie].x;
			const double& y = elements[ie].y;
			std::string str = 
				"Warning: \"NaN\" detected in element (x=" + std::to_string(x) + ", y=" + std::to_string(y) 
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
			LogWriter::writeLog(str, 1);
			//break;
		}
	}
	return is_nan;
}

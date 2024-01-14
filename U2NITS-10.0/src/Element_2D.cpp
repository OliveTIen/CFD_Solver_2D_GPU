#include "Element_2D.h"
#include "FVM_2D.h"
#include "Edge_2D.h"

//�ú����Ծ�̬��Ա��ʼ������Ž��ֲ�namespace
namespace {
    Eigen::MatrixXi iniSiti_Q4() {
        Eigen::MatrixXi si_ti(2, 4);//(ά��, ������)
        si_ti <<
            -1, +1, +1, -1,//si
            -1, -1, +1, +1;//ti
        return si_ti;
    }

    Eigen::MatrixXd iniGaussPoint_Q4() {
        //��˹���ֵ� 1/sqrt(3)=0.577350269189625(7,8), doubleȡ0.577350269189626
        Eigen::MatrixXd GaussPointMatrix(8, 3);//(������, ά��)
        double ops3 = 0.577350269189626;//onePerSqrt3
        GaussPointMatrix <<
            -ops3, -ops3,
            ops3, -ops3,
            ops3, ops3,
            -ops3, ops3;
        return GaussPointMatrix;
    }
}

//Eigen::MatrixXi Element_Q4::si_ti = iniSiti_Q4();
//Eigen::MatrixXd Element_Q4::GaussPointMatrix = iniGaussPoint_Q4();

double Element_T3::calArea(FVM_2D* f) {
    //�������������
    double xn[3]{}, yn[3]{};
    for (int i = 0; i < 3; i++) {
        xn[i] = f->getNodeByID(nodes[i])->x;
        yn[i] = f->getNodeByID(nodes[i])->y;
    }
    return 0.5 * abs(xn[0] * (yn[1] - yn[2]) + xn[1] * (yn[2] - yn[0]) + xn[2] * (yn[0] - yn[1]));
}

void Element_T3::inixy(FVM_2D* f) {
    //��ʼ����Ԫ��������
    double tmp_x = 0;
    double tmp_y = 0;
    for (int i = 0; i < 3; i++) {
        tmp_x += f->getNodeByID(nodes[i])->x;
        tmp_y += f->getNodeByID(nodes[i])->y;
    }
    x = tmp_x / 3.0;
    y = tmp_y / 3.0;
    //hasCalxyed = 1;
}

std::vector<Element_T3*> Element_T3::findNeighbor() {
    //ʹ��ǰ�豣֤pEdges�����µģ���pEdges��pElement_L/R�����µ�
    std::vector<Element_T3*> pNeighborElements(3);
    for (int iedge = 0; iedge < 3; iedge++) {
        if (pEdges[iedge] == nullptr) {
            std::cout << "Error: Element_T3::pEdges[" << iedge << "] is uninitialized(Element ID: " << this->ID << ")" << std::endl;
            pNeighborElements[iedge] = nullptr;
            //continue;//��������ѭ��
        }
        else {
            if (pEdges[iedge]->pElement_L == this) {
                if (pEdges[iedge]->pElement_R == nullptr) {
                    pNeighborElements[iedge] = nullptr;
                }
                else {
                    pNeighborElements[iedge] = pEdges[iedge]->pElement_R;
                }
            }
            else if(pEdges[iedge]->pElement_R == this) {
                if (pEdges[iedge]->pElement_L == nullptr) {
                    pNeighborElements[iedge] = nullptr;
                }
                else {
                    pNeighborElements[iedge] = pEdges[iedge]->pElement_L;
                }
            }
        }
    }
    return pNeighborElements;
}

double Element_T3::calDistanceFromNearestNeighbor(FVM_2D* f) {
    std::vector<Element_T3*> n = findNeighbor();
    double dis = 1e10;
    double tmpx, tmpy, tmpdis;
    //this->calxy(f);��ȡ�ļ�ʱ�Ѿ�calxy��
    for (int in = 0; in < n.size(); in++) {
        if (n[in] != nullptr) {
            //n[in]->calxy(f);��ȡ�ļ�ʱ�Ѿ�calxy��
            tmpx = n[in]->x;
            tmpy = n[in]->y;
            tmpdis = sqrt((tmpx - x) * (tmpx - x) + (tmpy - y) * (tmpy - y));
            dis = (std::min)(dis, tmpdis);
        }
    }
    return dis;
}

void Element_T3::updateSlope_Barth(FVM_2D* f) {
    //����UxUy����������

    //3���ھ� ������
    //   |x0-x y0-y|       |U0i-Ui|
    //   |x1-x y1-y| |Uxi|=|U1i-Ui|
    //   |x2-x y2-y| |Uyi| |U2i-Ui|
    //   
    
    //2���ھ�
    //   
    //   |x0-x y0-y| |Uxi|=|U0i-Ui|
    //   |x1-x y1-y| |Uyi| |U1i-Ui|
    //   
    
    //1���ھ� ���� ��������
    //   
    //   |x0-x y0-y| |Uxi|=|U0i-Ui|
    //               |Uyi|
    //   

    std::vector<Element_T3*>neighbors_origin = findNeighbor();
    std::vector<Element_T3*>neighbors;
    //���nullptr
    for (int i = 0; i < 3; i++) {
        if (neighbors_origin[i] != nullptr) {
            //neighbors_origin[i]->calxy(f);//��ʼ������ calxy!!!
            neighbors.push_back(neighbors_origin[i]);
        }
    }
    //this->calxy(f);// ��ȡ�ļ�ʱ�Ѿ�calxy
    const int nNeighbor = (int)neighbors.size();

    for (int i_ = 0; i_ < 4; i_++) {
        Eigen::MatrixXd xy(nNeighbor, 2);
        Eigen::VectorXd Uxy(2);//������
        Eigen::VectorXd UU(nNeighbor);
        //initialize xy, UU
        for (int j = 0; j < nNeighbor; j++) {
            xy(j, 0) = neighbors[j]->x - this->x;//xj-x �Ѿ�calxy
            xy(j, 1) = neighbors[j]->y - this->y;//yj-y
            UU(j) = neighbors[j]->U[i_] - this->U[i_];//Uji-Ui
        }
        //solve
        if (nNeighbor == 3) {
            Uxy = xy.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(UU);//����ֵ�ֽ�
            this->Ux[i_] = Uxy[0];
            this->Uy[i_] = Uxy[1];
        }
        else if (nNeighbor == 2) {
            Uxy = xy.householderQr().solve(UU);//�����Է�����
            this->Ux[i_] = Uxy[0];
            this->Uy[i_] = Uxy[1];
        }
        else {
            std::cout << "δ��ɣ�Element_T3::calSlope_Barth(FVM_2D* f) nNeighbor==1" << std::endl;
            //ָ���ݶȷ���
        }

    }
    //������
    restructor_in_updateSlope_Barth(f);
}

void Element_T3::restructor_in_updateSlope_Barth(FVM_2D* f) {
    //����������������Ux,Uy�����������������������㴦��Ux,Uy�����н��ԣ�������ȫ�����н���
    //����ƫ�����½�
    std::vector<Element_T3*> neighbors = findNeighbor();
    double UU[3][4]{};//�ھӺ���ֵ��������ֵ�Ĳ�
    double UUup = 0;
    double UUdown = 0;
    for (int in = 0; in < neighbors.size(); in++) {
        if (neighbors[in] != nullptr) {
            for (int j = 0; j < 4; j++) {
                UU[in][j] = neighbors[in]->U[j] - this->U[j];
                UUup = (std::max)(UUup, UU[in][j]);
                UUdown = (std::min)(UUdown, UU[in][j]);
            }
        }
    }
    //���㶥��ƫ��
    double U_node[3][4]{};
    double UU_node[3][4]{};//���㺯��ֵ��������ֵ�Ĳ�
    double UUup_node = 0;
    double UUdown_node = 0;
    for (int i_node = 0; i_node < 3; i_node++) {
        Node_2D* pNode = f->pNodeTable[nodes[i_node]];
        get_U(pNode->x, pNode->y, U_node[i_node]);
        for (int j = 0; j < 4; j++) {
            UU_node[i_node][j] = U_node[i_node][j] - this->U[j];
            UUup_node = (std::max)(UUup_node, UU_node[i_node][j]);
            UUdown_node = (std::min)(UUdown_node, UU_node[i_node][j]);
        }
    }
    //����Ux, Uy
    double ratio = 1;
    if (UUup_node > UUup) {
        ratio = (std::max)(ratio, UUup_node / UUup);
    }
    if (UUdown_node < UUdown) {
        ratio = (std::max)(ratio, UUdown_node / UUdown);
    }

    for (int j = 0; j < 4; j++) {
        Ux[j] /= ratio;
        Uy[j] /= ratio;
    }
}

void Element_T3::get_U(double xpoint, double ypoint, double* _U) {
    //�����ع�
    if (GlobalPara::space::flag_reconstruct == _REC_constant) {
        for (int j = 0; j < 4; j++) {
            _U[j] = U[j];//�˴�UΪ��Ա����
        }
    }
    //�����ع�
    else if (GlobalPara::space::flag_reconstruct == _REC_linear) {
        for (int j = 0; j < 4; j++) {
            _U[j] = U[j] + Ux[j] * (xpoint - x) + Uy[j] * (ypoint - y);
        }
    }
}

Eigen::Vector4d Element_T3::get_U(double xpoint, double ypoint) {
    //��ʼ��_U����
    double _U[4];
    get_U(xpoint, ypoint, _U);//��һ�����غ���
    //����ת����
    Eigen::Vector4d _U_vector;
    for (int j = 0; j < 4; j++) {
        _U_vector[j] = _U[j];
    }
    return _U_vector;
}

std::vector<double> Element_T3::U2uv(const Eigen::Vector4d& Uc) {
    //�غ�����,��u,��v,��E
    std::vector<double> uv;//���غ���u v
    uv[0] = Uc[1]/Uc[0];
    uv[1] = Uc[2]/Uc[0];
    return uv;
}

void Element_T3::generateElementEdge(FVM_2D* f) {
    generateElementEdge_registerSingle(f, nodes[0], nodes[1], 0);
    generateElementEdge_registerSingle(f, nodes[1], nodes[2], 1);
    generateElementEdge_registerSingle(f, nodes[2], nodes[0], 2);
}

void Element_T3::generateElementEdge_registerSingle(FVM_2D* f, int ID_0, int ID_1, int iEdge) {
    //���룺��ע��edge������node���
    //������
    //1.���edge������edge��Element_R��Ϣ�����޸�edges��pEdgeTable
    //2.��Element.pEdgesע��
    
    //����f_edges����ĳԪ�ص�nodes�뵱ǰ��ȷ��nodes��ͬ(�����Ⱥ�)��˵���Ѿ�ע����
    std::vector<Edge_2D>& f_edges = f->edges;//���õı�����ָ�볣��
    std::vector<Edge_2D*>& f_pEdgeTable = f->pEdgeTable;
    for (int i_edge = 0; i_edge < f_edges.size(); i_edge++) {
        Edge_2D& edge = f_edges[i_edge];//���ص��ľֲ������ķ��������ͷŵĿ�������Ϊ�������ܴ�����
        //if ((ID_0 == edge.nodes[0] && ID_1 == edge.nodes[1]) || (ID_0 == edge.nodes[1] && ID_1 == edge.nodes[0])) {
        //�Ѿ�ע��������ֻ�貹��element_R��Ϣ
        //����˵�����Ѿ�ע������Ҷ�����ʱ����ת������ע���edge����ע���edge��Ȼ�Ƿ����෴��
        if (ID_0 == edge.nodes[1] && ID_1 == edge.nodes[0]) {
            edge.pElement_R = this;
            //����pEdges
            pEdges[iEdge] = &(f_edges[i_edge]);
            //ֱ����ֹ
            return;
        }
    }
    //δ��Ѱ����˵������ע��� L R H ID nodes
    Edge_2D new_edge;
    new_edge.pElement_L = this;
    new_edge.ID = (int)f_edges.size() + 1;//���ϸ���inp�ļ��е�ID��
    int edgeID = new_edge.ID;
    new_edge.nodes[0] = ID_0;
    new_edge.nodes[1] = ID_1;
    f_edges.push_back(new_edge);
    Edge_2D* pNewEdge = &(f_edges[f_edges.size() - 1]);
    //����table
    if (f_pEdgeTable.size() < edgeID + 1) {
        f_pEdgeTable.resize(edgeID + 1);//��ID�ţ���ҪID+1�Ŀռ�
        f_pEdgeTable[edgeID] = pNewEdge;
    }
    //����pEdges
    pEdges[iEdge] = pNewEdge;//ȡ��&(f_edges[f_edges.size() - 1])����˲�����ֲ�����������
}

void Element_T3::iniPEdges(FVM_2D* f) {
    pEdges[0] = f->isEdgeExisted(nodes[0], nodes[1]);
    pEdges[1] = f->isEdgeExisted(nodes[1], nodes[2]);
    pEdges[2] = f->isEdgeExisted(nodes[2], nodes[0]);
}


double Element_T3::calLambda(const double gamma) {
    double rho = U[0];
    double rho_u = U[1];//rho*u
    double rho_v = U[2];
    double rho_E = U[3];
    double E = rho_E / rho;
    double u = rho_u / rho;
    double v = rho_v / rho;
    double V2 = u * u + v * v;//|U|^2
    double p = 0.5 * rho * (gamma - 1) * (2 * E - V2);
    return sqrt(V2) + sqrt(gamma * p / rho);
}

double Element_T3::calLambdaFlux(FVM_2D* f) {
    double LambdaC = 0;
    for (int ie = 0; ie < 3; ie++) {
        //eU
        double ex, ey;
        double eU[4]{};//���е㴦��,��u,��v,��E
        pEdges[ie]->getxy(f, ex, ey);
        get_U(ex, ey, eU);
        //en
        std::vector<double> en = pEdges[ie]->getDirectionN(f);
        //eabs
        double u = eU[1] / eU[0];
        double v = eU[2] / eU[0];
        double E = eU[3] / eU[0];
        double eabs = abs(u * en[0] + v * en[1]);//|euv��en|
        //ec
        double V2 = u * u + v * v;
        double& rho = eU[0];
        double p = 0.5 * rho * (Constant::gamma - 1) * (2 * E - V2);
        double ec = sqrt(Constant::gamma * p / rho);
        double dl = pEdges[ie]->getLength(f);
        LambdaC += (eabs + ec) * dl;
    }
    return LambdaC;
}

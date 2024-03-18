#ifndef _CONSTEXPRS_H_
#define _CONSTEXPRS_H_

constexpr auto _NULL = 0;
//������� 1-999
constexpr auto _OUT_u = 1;
constexpr auto _OUT_v = 2;
constexpr auto _OUT_rho = 4;
constexpr auto _OUT_p = 8;
constexpr auto _OUT_T = 16;
constexpr auto _OUT_Unon_all = 101;
//�������� 1001-1999
constexpr auto _EQ_convective = 1001;//������ɢ����
constexpr auto _EQ_euler = 1002;//ŷ������
constexpr auto _EQ_NS = 1003;
//��������� 2001-2999
constexpr auto _SOL_LocalLaxFriedrichs = 2001;//��������� Local Lax-Friedrichs;
constexpr auto _SOL_Roe = 2002;
//�ع���ʽ 3001-3099
constexpr auto _REC_constant = 3001;//�����ع�
constexpr auto _REC_linear = 3002;  //�����ع�
//ʱ���ƽ���ʽ 3101-3199 evolution
constexpr auto _EVO_explicit = 3101;//��ʽʱ���ƽ�
constexpr auto _EVO_rk3 = 3103;//RK3
constexpr auto _EVO_rk5 = 3105;
//������ 4001-4999 limiter
constexpr auto _LIM_minmod = 4001;
//��Ԫ���� 5001-5999
constexpr auto _ELE_D2 = 5001;//2�ڵ��߶�Ԫ
constexpr auto _ELE_T3 = 5101;//3�ڵ�����Ԫ
constexpr auto _ELE_Q4 = 5102;//4�ڵ��ı�Ԫ
//�߽����� 6001-6999
constexpr auto _BC_wall_nonViscous = 6001;//�̱ڣ���ճ(�л���)
constexpr auto _BC_wall_adiabat = 6002;//�̱ڣ��޻��ƣ�����
constexpr auto _BC_wall_isothermal = 6003;//�̱ڣ��޻��ƣ�����
constexpr auto _BC_inlet = 6011;//���
constexpr auto _BC_outlet = 6012;//����
constexpr auto _BC_inf = 6013;//Զ���߽�
constexpr auto _BC_symmetry = 6021;//�ԳƱ߽�
constexpr auto _BC_periodic_0 = 6100;//���ڱ߽�
constexpr auto _BC_periodic_1 = 6101;
constexpr auto _BC_periodic_2 = 6102;
constexpr auto _BC_periodic_3 = 6103;
constexpr auto _BC_periodic_4 = 6104;
constexpr auto _BC_periodic_5 = 6105;
constexpr auto _BC_periodic_6 = 6106;
constexpr auto _BC_periodic_7 = 6107;
constexpr auto _BC_periodic_8 = 6108;
constexpr auto _BC_periodic_9 = 6109;


#endif

#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>



namespace GlobalPara {
	namespace constant {
		extern double const R;// ���峣��
		extern double const PI;
		extern double T0;// ��ƽ���¶Ȳο�ֵ
		extern double p0;// ��ƽ��ѹ���ο�ֵ
		extern double c0;// ��ƽ�����ٲο�ֵ
		extern double gamma;
		extern double epsilon;
		extern double Re;
		extern double Pr;
		extern double mu;// ����ճ��ϵ��
	}

	//�ļ���
	namespace basic {
		extern int dimension;
		extern bool _continue;
		extern std::string filename;
		extern std::string meshFileType;
		extern bool useGPU;// ������extern�����򱨴��ض��塱
	}
	namespace space {
		namespace _1D {
			extern int nElement;
			extern double x1;
			extern double x2;
		}
		extern int flag_reconstruct;// �ع�������[todo]��ͨ�����췽���׻�����fun3d��δ�ҵ���ѡ��
	}
	namespace time {
		extern double CFL;
		extern double T;
		extern double residual;// �в�����
		extern int time_advance;// ʱ���ƽ���ʽ
		// �������ݴ�ContinueFile�ж�ȡ������Toml
		extern double t_previous;// ����ʱ����ʼʱ�� ��readContineFile()�ж�ȡ
		extern int istep_previous;// ����ʱ����ʼ�� ��readContineFile()�ж�ȡ
	}
	namespace physicsModel {
		extern int equation;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				extern bool use_ruvp;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
			namespace outlet {
				extern bool use_ruvp;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
			namespace inf {
				extern bool use_ruvp;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
		}
	}
	namespace initialCondition {
		extern int type;
	}
	namespace output {
		extern int step_per_print;
		extern int step_per_output_field;
		extern int step_per_output_hist;
		extern bool output_var_ruvp[4];
		extern int autosaveFileNum;
		extern int maxIteration;
	}
	namespace inviscid_flux_method {
		extern int flux_conservation_scheme;// ��ճͨ�� �غ��ʽ ��������������� LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// ͨ��������
	}
}


#endif
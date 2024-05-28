#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H
#include <string>
#include <vector>
/*
Ŀǰ�������ֻ࣬�þ�̬����
��ֻ����һ�ֲ��ϡ���������vector
[ע��]�����޸�ΪGPU����ʱ��GPU����ʱ��SDeviceParaҪ�����µ�mu
*/

class MaterialManager {
private:
	struct FluidMaterial {
	public:
		std::string name = "air";
		double mu0 = 17.9e-6; // ����ճ��ϵ���ο�ֵ
		FluidMaterial(std::string _name, double _mu0_dynamic_viscosity) {
			name = _name;
			mu0 = _mu0_dynamic_viscosity;
		}
		// �����ο��ܶ�rho���ο��ٶ�u���ο�����L����ŵ��Re������mu0���洢
		FluidMaterial(std::string _name, double rho, double u, double L, double Re) {
			name = _name;
			mu0 = rho * u * L / Re;
		}
		FluidMaterial() {};
		//FluidMaterial(FluidMaterial&) = delete;
		//FluidMaterial(FluidMaterial&&) = delete;
	};

	static MaterialManager* p_instance;
	std::vector<FluidMaterial> m_materials;

	MaterialManager() {};

public:
	static MaterialManager* getInstance();
	void addMaterial(FluidMaterial material);
	// ���س����ã�������Ϊ��ֵ�����ɱ��޸�
	const FluidMaterial& getMaterial(size_t index = 0);
	void initialize_using_config(double rho_ref, double U_ref, double L_ref);
};

#endif // !MATERIAL_MANAGER_H

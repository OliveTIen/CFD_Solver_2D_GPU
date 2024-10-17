#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "ShaderParticle.h"
#include "../io/ConsoleProxy.h"
#include "../object/Object.h"

constexpr double PI = 3.1415926535897932;

class ParticleSystem :public  Object {
public:
	// initialize
	std::vector<glm::vec3> vertices;
	glm::uvec3 initialNum = glm::uvec3(5,5,5); // number of particles in each dimension
	glm::vec3 initialCube[2]{ glm::vec3(0.3f,5.5f,0.3f),glm::vec3(0.7f,9.9f,0.7f) };// cube range
	// background grid
	glm::vec3 bgGridCube[2]{ glm::vec3(0.0f,0.0f,0.0f),glm::vec3(1.0f,10.0f,1.0f) };
	glm::uvec3 bgGridElemNum = glm::uvec3(8, 8, 8);// element num
	glm::vec3 bgGridElemSize;
	std::vector<glm::ivec3> particleGridIds;
	std::vector<std::vector<int>> gridParticleIds;
	std::vector<std::vector<std::pair<int, float>>> particleNeighbors;
	// physical data
	std::vector<glm::vec3> velocities;
	std::vector<glm::vec3> accelerations;
	std::vector<float> densities;
	std::vector<float> pressures;
	std::vector<float> dt_queue = std::vector<float>(5, 0);

	bool enable_rendering = false;

	struct {
		bool useBgGrid = false;
		bool printOnConsole = false;
		bool applyRigidCollision = false;
		bool applyPressureForce = true;
		bool applyViscousForce = true;
		int bounceType = 1;// 0-acceleration 1-speedChange

		float gravity = 9.8f;
		glm::vec3 gravity_direction{0, -1, 0};

		float particle_mass = 1.0f;


		float smooth_length = 0.05f;
		float sphere_radius = 0.01f;
		float rigidCollision_radius = 0.08f;
		float rigidCollision_spring = 50.0f;
		float air_drag = 1.0;
		float bounce0_spring = 1.0e4;// apply penalty force to implement bouncing
		float bounce0_norm_drag = 1.0;// drag in normal direction
		float bounce0_wall_friction = 0.1;
		float bounce1_damping_coefficient = 0.9f;
		float bounce1_shearSpeedAddition = 0.0f;// convert normal speed to shear speed
		float time_step_length = 0.01;
		
		float m_standard_distance = 0.19f;// m_standard_distance = pow(particle_mass / m_restDensity, 1.f/3.f);
		float m_rest_density = particle_mass / pow(m_standard_distance, 3);// pD^3 = m / restDensity
		float m_gasConstantK = 1.0f;
		float m_epsilon = 1e-3 * m_standard_distance;
		float m_viscosity = 1.0f;
		bool m_clip_smooth_function = true;
		
	} physicalParameters;

	// rendering data
	ShaderParticle* shader;
	unsigned int vao;
	unsigned int vbo;
	glm::vec4 userColor = glm::vec4(0.5f, 0.5f, 1.0f, 1.0f);
	glm::mat4 matrix_model = glm::mat4(1.0f);
	float pointRadius_render = 40;// point size in world space
	float refDistance = 1.0f;// the distance where gl_PointSize = pointRadius
	COORD current_cursor_position;

	// functions
	ParticleSystem() { init(); }
	virtual ~ParticleSystem() override { cleanup(); };
	virtual void init() override;
	virtual void reset() override { resetParticles(); }
	virtual void update(float dt) override;
	virtual void render(const glm::mat4& projection, const glm::mat4& view) const override;
	virtual void cleanup() override;

	static float getDensity(float particleMass, float smoothingLength);
	static float getKernelPoly6(float h);
	static float getWPoly6(float h, float r);
	//static float getRho(float m, float h, float r);
	
	static float getPressureForce(float particleMass, float smoothingLength);
	static float getViscousForce(float particleMass, float dynamicViscosity, float smoothingLength);

private:
	bool buffer_allocated = false;

	void resetParticles();
	int find_biggest_pow_of_2(int upper_limit);
	void initializeBgGrid();
	void updateParticleGridIds();
	std::vector<std::pair<int,float>> searchNeighborsInRadius(int particleId, float radius);
	bool in_bgGridCube(const glm::vec3& point, float ball_radius);
	// flag can be -1,0 or 1, and 0 means is in cube
	void in_bgGridCube(const glm::vec3& point, float ball_radius, int& flag_x, int& flag_y, int& flag_z);
	void setBoundaryCondition_Wall(int particle, float ball_radius, float dt);
	void setParameters();
};
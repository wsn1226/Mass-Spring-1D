#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::MatrixXd &H, const Eigen::VectorXd &q, double stiffness)
{
    H(0, 0) = stiffness;
}
#include <dV_spring_particle_particle_dq.h>
#include <iostream>

void dV_spring_particle_particle_dq(Eigen::VectorXd &dV, const Eigen::VectorXd &q, double stiffness)
{
    dV(0) = stiffness * q(0);
}
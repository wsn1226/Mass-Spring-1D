#include <Eigen/Dense>

// Input:
//   q - generalized coordiantes for the mass-spring system
//   qdot - generalized velocity for the mass spring system
//   dt - the time step in seconds
//   mass - the mass
//   force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//   stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.
// Output:
//   q - set q to the updated generalized coordinate using Backward Euler time integration
//   qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template <typename FORCE, typename STIFFNESS>
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force, STIFFNESS &stiffness)
{
    Eigen::VectorXd f(1);
    Eigen::MatrixXd H(1, 1);
    force(f, q, qdot);
    stiffness(H, q, qdot);
    double coe = mass / (mass - H(0, 0) * dt * dt);

    double q_new = (q(0) + qdot(0) * dt) * coe;
    double qdot_new = (qdot(0) + (f(0) / mass) * dt) * coe;

    q(0) = q_new;
    qdot(0) = qdot_new;
}
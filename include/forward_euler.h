#include <Eigen/Dense>

// Input:
//   q - generalized coordiantes for the mass-spring system
//   qdot - generalized velocity for the mass spring system
//   dt - the time step in seconds
//   mass - the mass
//   force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
// Output:
//   q - set q to the updated generalized coordinate using Forward Euler time integration
//   qdot - set qdot to the updated generalized velocity using Forward Euler time integration

template <typename FORCE>
inline void forward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force)
{

    Eigen::VectorXd f(1);
    force(f, q, qdot);

    double q_new = q(0) + qdot(0) * dt;
    double qdot_new = qdot(0) + (f(0) / mass) * dt;

    q(0) = q_new;
    qdot(0) = qdot_new;

    std::cout << q(0) << std::endl;
}
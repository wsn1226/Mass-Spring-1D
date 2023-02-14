// Input:
//   q - generalized coordiantes for the mass-spring system
//   qdot - generalized velocity for the mass spring system
//   dt - the time step in seconds
//   mass - the mass
//   force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
// Output:
//   q - set q to the updated generalized coordinate using Runge-Kutta time integration
//   qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template <typename FORCE>
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force)
{
    Eigen::VectorXd f(1);
    Eigen::VectorXd q_new(1);
    Eigen::VectorXd qdot_new(1);
    double k_v1, k_x1, k_v2, k_x2, k_v3, k_x3, k_v4, k_x4;

    force(f, q, qdot);
    k_v1 = (f(0) / mass);
    k_x1 = qdot(0);
    q_new(0) = q(0) + 0.5 * dt * k_x1;
    qdot_new(0) = qdot(0) + 0.5 * dt * k_v1;

    force(f, q_new, qdot_new);
    k_v2 = (f(0) / mass);
    k_x2 = qdot_new(0);
    q_new(0) = q(0) + 0.5 * dt * k_x2;
    qdot_new(0) = qdot(0) + 0.5 * dt * k_v2;

    force(f, q_new, qdot_new);
    k_v3 = (f(0) / mass);
    k_x3 = qdot_new(0);
    q_new(0) = q(0) + dt * k_x3;
    qdot_new(0) = qdot(0) + dt * k_v3;

    force(f, q_new, qdot_new);
    k_v4 = (f(0) / mass);
    k_x4 = qdot_new(0);

    q(0) = q(0) + (1.0 / 6.0) * dt * (k_x1 + 2 * k_x2 + 2 * k_x3 + k_x4);
    qdot(0) = qdot(0) + (1.0 / 6.0) * dt * (k_v1 + 2 * k_v2 + 2 * k_v3 + k_v4);
}
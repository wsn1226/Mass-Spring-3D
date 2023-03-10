#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <EigenTypes.h>

// Input:
//   q - generalized coordiantes for the mass-spring system
//   qdot - generalized velocity for the mass spring system
//   dt - the time step in seconds
//   mass - the mass matrix
//   force(q, qdot) - a function that computes the force acting on the mass-spring system. This takes q and qdot as parameters.
//   stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.
//   tmp_force - scratch space to collect forces
//   tmp_stiffness - scratch space to collect stiffness matrix
// Output:
//   q - set q to the updated generalized coordinate using linearly implicit time integration
//   qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template <typename FORCE, typename STIFFNESS>
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                                    const Eigen::SparseMatrixd &mass, FORCE &force, STIFFNESS &stiffness,
                                    Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness)
{

    stiffness(tmp_stiffness, q, qdot);
    force(tmp_force, q, qdot);
    Eigen::SparseMatrix<double> A(mass.rows(), mass.cols());
    A = mass - dt * dt * tmp_stiffness; // M-dt2*k
    // std::cout << tmp_stiffness.coeff(0, 0) << std::endl;
    Eigen::VectorXd x, b;
    b = mass * qdot + dt * tmp_force;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    x = solver.solve(b);
    q = q + dt * x;
    qdot = x;

    // Eigen::VectorXd q_new = q + qdot * dt;
    // Eigen::VectorXd qdot_new = qdot + tmp_force * dt;

    /*     q = q_new;
        qdot = qdot_new; */

    // std::cout << qdot(5) << std::endl;
}
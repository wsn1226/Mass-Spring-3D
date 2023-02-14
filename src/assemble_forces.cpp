#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0,
                     double mass, double k)
{
    Eigen::Vector6d indi_force;
    Eigen::Vector3d q0;
    Eigen::Vector3d q1;
    Eigen::Vector3d indi_force0;
    Eigen::Vector3d indi_force1;
    Eigen::Vector3d g_force;
    f.resize(V.rows() * 3);
    f.setZero();
    for (int i = 0; i < E.rows(); i++)
    {
        int ind0 = E(i, 0);
        int ind1 = E(i, 1);
        q0 << q(3 * ind0), q(3 * ind0 + 1), q(3 * ind0 + 2);
        q1 << q(3 * ind1), q(3 * ind1 + 1), q(3 * ind1 + 2);
        // q0 = Eigen::Vector3d(V(ind0, 0), V(ind0, 1), V(ind0, 2));
        // q1 = Eigen::Vector3d(V(ind1, 0), V(ind1, 1), V(ind1, 2));
        dV_spring_particle_particle_dq(indi_force, q0, q1, l0[i], k);
        indi_force0 << indi_force[0], indi_force[1], indi_force[2];
        indi_force1 << indi_force[3], indi_force[4], indi_force[5];
        dV_gravity_particle_dq(g_force, mass, Eigen::Vector3d(0., -9.8, 0.));
        indi_force0 -= g_force;
        indi_force1 -= g_force;
        f[3 * ind0] += indi_force0[0];
        f[3 * ind0 + 1] += indi_force0[1];
        f[3 * ind0 + 2] += indi_force0[2];

        f[3 * ind1] += indi_force1[0];
        f[3 * ind1 + 1] += indi_force1[1];
        f[3 * ind1 + 2] += indi_force1[2];
    }
    f *= -1.0;
    //    std::cout << "heyhey0" << std::endl;
};
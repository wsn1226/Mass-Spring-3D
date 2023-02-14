#include <d2V_spring_particle_particle_dq2.h>
#include <iostream>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness)
{
    Eigen::Matrix66d BTB, A, B, C, a;
    BTB << Eigen::MatrixXd::Identity(3, 3), -1.0 * Eigen::MatrixXd::Identity(3, 3), -1.0 * Eigen::MatrixXd::Identity(3, 3), Eigen::MatrixXd::Identity(3, 3);
    Eigen::Vector3d rel_pos = q1 - q0;
    Eigen::Vector6d q;
    q << q0, q1;

    double current_len = sqrt(rel_pos.transpose() * rel_pos); // sqrt(qtBTBq)
    double rel_len = current_len - l0;

    A = BTB;
    A *= stiffness; // A = kBTB

    B.setIdentity();
    B /= current_len;

    C = (-1.0 / std::pow(current_len, 3)) * q * q.transpose() * BTB;

    a = B + C;

    H = A - l0 * A * a;
}
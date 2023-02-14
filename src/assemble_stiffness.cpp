#include <assemble_stiffness.h>
#include <iostream>

typedef Eigen::Triplet<double> T;
typedef std::tuple<uint16_t, uint16_t> tpl;

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0,
                        double k)
{
    Eigen::Vector3d q0, q1;
    std::vector<T> K_entries;
    std::map<tpl, double> Hmap;
    K.resize(3 * V.rows(), 3 * V.rows());
    K.setZero();
    for (int i = 0; i < E.rows(); i++)
    {
        int ind0 = E(i, 0);
        int ind1 = E(i, 1);
        q0 << q(3 * ind0), q(3 * ind0 + 1), q(3 * ind0 + 2);
        q1 << q(3 * ind1), q(3 * ind1 + 1), q(3 * ind1 + 2);
        Eigen::Matrix66d PerSpr_H;
        d2V_spring_particle_particle_dq2(PerSpr_H, q0, q1, l0[i], k);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                tpl index = tpl(3 * ind0 + i, 3 * ind0 + j);
                if ((Hmap.find(index) == Hmap.end()))
                {
                    Hmap[index] = PerSpr_H(i, j);
                }
                else
                {
                    Hmap[index] += PerSpr_H(i, j);
                }
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                tpl index = tpl(3 * ind0 + i, 3 * ind1 + j);
                if ((Hmap.find(index) == Hmap.end()))
                {
                    Hmap[index] = PerSpr_H(i, j + 3);
                }
                else
                {
                    Hmap[index] += PerSpr_H(i, j + 3);
                }
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                tpl index = tpl(3 * ind1 + i, 3 * ind0 + j);
                if ((Hmap.find(index) == Hmap.end()))
                {
                    Hmap[index] = PerSpr_H(i + 3, j);
                }
                else
                {
                    Hmap[index] += PerSpr_H(i + 3, j);
                }
            }
        }
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                tpl index = tpl(3 * ind1 + i, 3 * ind1 + j);
                if ((Hmap.find(index) == Hmap.end()))
                {
                    Hmap[index] = PerSpr_H(i + 3, j + 3);
                }
                else
                {
                    Hmap[index] += PerSpr_H(i + 3, j + 3);
                }
            }
        }
    }

    for (auto const &ele : Hmap)
    {
        K_entries.push_back(T(std::get<0>(ele.first), std::get<1>(ele.first), ele.second));
    }

    // std::cout << "heyhey" << std::endl;
    K.setFromTriplets(K_entries.begin(), K_entries.end());
    K *= -1.0;
    // std::cout << "heyheyhey" << std::endl;
};
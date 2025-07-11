#include "Random.h"

namespace Meso {
    namespace Random {
        std::mt19937 sequential_gen(5489U);
        std::uniform_real_distribution<real> uniform_unit_real(0.0, 1.0);

        int RandInt(int a, int b)
        {
            std::uniform_int_distribution<int> uid(a, b);
            int ret;
#pragma omp critical
            {
                ret = uid(sequential_gen);
            }
            return ret;
        }

        real Random(void)
        {
            real ret;
#pragma omp critical
            {
                ret = uniform_unit_real(sequential_gen);
            }
            return ret;
        }

        real Uniform(real a, real b) {
            return a + (b - a) * Random();
        }

        VectorXd Random_VectorXd(int n, real a, real b) {
            VectorXd x(n);
            for (int i = 0; i < n; i++) x[i] = Uniform(a, b);
            return x;
        }
        VectorXi Random_VectorXi(int n, int a, int b) {
            VectorXi x(n);
            for (int i = 0; i < n; i++) x[i] = RandInt(a, b);
            return x;
        }
        real Random_Sign(void) {
            return Random() < (real)0.5 ? (real) - 1 : (real) 1;
        }
    }
}

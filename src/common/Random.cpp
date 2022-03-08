#include "Random.h"

namespace Meso {
    namespace Random {
        std::mt19937 sequential_gen(5489U);
        std::uniform_real_distribution<real> uniform_unit_real(0.0, 1.0);

        int RandInt(int a, int b)
        {
            std::uniform_int_distribution<int> uid(a, b);
#pragma omp critical
            return uid(sequential_gen);
        }

        real Random(void)
        {
#pragma omp critical
            return uniform_unit_real(sequential_gen);
        }

        real Uniform(real a, real b) {
            return a + (b - a) * Random();
        }
    }
}

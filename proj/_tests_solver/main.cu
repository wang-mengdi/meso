#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "PoissonMapping.h"
#include "PoissonTests.h"
#include "MultigridTests.h"
#include "ExteriorDerivativeTests.h"

using namespace Meso;

int main(){

    //A visual test

    Test_MGPCG_Dirichlet<double>(Vector2i(256, 256), false);
    Test_MGPCG_Dirichlet<double>(Vector2i(512, 512), false);
    Test_MGPCG_Dirichlet<double>(Vector3i(128, 128, 128), false);
    Test_MGPCG_Dirichlet<double>(Vector3i(256, 256, 256), false);

   Test_MGPCG_Dirichlet_Neumann<double>(Vector2i(256, 256), false);
    Test_MGPCG_Dirichlet_Neumann<double>(Vector2i(512, 512), false);
    Test_MGPCG_Dirichlet_Neumann<double>(Vector3i(128, 128, 128), false);
    Test_MGPCG_Dirichlet_Neumann<double>(Vector3i(256, 256, 256), false);

    Test_MGPCG_Neumann<double>(Vector2i(256, 256), false);
    Test_MGPCG_Neumann<double>(Vector2i(512, 512), false);
    Test_MGPCG_Neumann<double>(Vector3i(256, 256, 256), false);

    return 0;
}
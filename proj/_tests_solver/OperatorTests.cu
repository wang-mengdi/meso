#include "OperatorTests.h"
#include "AuxFunc.h"
//#include "Random.h"

namespace Meso {

	void Test_Coarsener2(const Vector2i counts)
	{
		//Info("counts: {}", counts);
		Grid<2> grd(counts);
		//Field<bool, 2> finer_data(Grid<2>(counts));
		Field<bool, 2> finer_data(grd);

		int fnx = finer_data.grid.counts[0], fny = finer_data.grid.counts[1];
		for (int i = 0; i < fnx; i++) {
			real frac = (i + 0.0) / fnx;
			for (int j = 0; j < fny; j++) {
				//bool value = (Random::Random() <= frac);
				//finer_data(Vector2i(i, j)) = value;
			}
		}

		//Field<bool, 2> coarser_data(Grid<2>(counts / 2));
		//int cnx = coarser_data.grid.counts[0], cny = coarser_data.grid.counts[1];
		//for (int i = 0; i < cnx; i++) {
		//	for (int j = 0; j < cny; j++) {
		//		bool value = false;
		//		for (int t0 = 0; t0 < 2; t0++) {
		//			for (int t1 = 0; t1 < 2; t1++) {
		//				Vector2i finer_sub(i * 2 + t0, j * 2 + t1);
		//				//if (finer_data.grid.Valid(finer_sub)) value |= finer_data(finer_sub);
		//			}
		//		}
		//		coarser_data(Vector2i(i, j)) = value;
		//	}
		//}
	}

}
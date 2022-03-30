#include "OperatorTests.h"
#include "AuxFunc.h"
#include "Random.h"

namespace Meso {

	void Test_Coarsener2(const Vector2i counts)
	{
		Field<bool, 2> finer_data{ Grid<2>(counts) };

		int fnx = finer_data.grid.counts[0], fny = finer_data.grid.counts[1];
		for (int i = 0; i < fnx; i++) {
			real frac = (i + 0.0) / fnx;
			for (int j = 0; j < fny; j++) {
				bool value = (Random::Random() <= frac);
				finer_data(Vector2i(i, j)) = value;
			}
		}

		Field<bool, 2> coarser_data(Grid<2>(counts / 2));
		int cnx = coarser_data.grid.counts[0], cny = coarser_data.grid.counts[1];
		for (int i = 0; i < cnx; i++) {
			for (int j = 0; j < cny; j++) {
				bool fixed = true;
				for (int t0 = 0; t0 < 2; t0++) {
					for (int t1 = 0; t1 < 2; t1++) {
						Vector2i finer_sub(i * 2 + t0, j * 2 + t1);
						if (finer_data.grid.Valid(finer_sub) && !finer_data(finer_sub)) fixed = false;
					}
				}
				coarser_data(Vector2i(i, j)) = fixed;
			}
		}

		FieldDv<bool, 2> finer_device; finer_device = finer_data;
		FieldDv<bool, 2> coarser_device{ coarser_data.grid };
		Coarsener<2>::Apply(coarser_device, finer_device);

		Field<bool, 2> coarsen_result; coarsen_result = coarser_device;

		if (ArrayFunc::Equals<bool>(coarsen_result.data, coarser_data.data)) {
			Pass("Test_Coarsener2 passed {}", counts);
		}
		else {
			Error("Test_Coarsener2 failed {}", counts);
			Field<bool, 2> finer_temp; finer_temp.Copy(finer_device);
			Info("finer_data: \n{}", finer_data);
			Info("coarser_data:\n{}", coarser_data);
			Info("finer on device:\n{}", finer_temp);
			Info("coarsen_result:\n{}", coarsen_result);
			
		}
	}

}
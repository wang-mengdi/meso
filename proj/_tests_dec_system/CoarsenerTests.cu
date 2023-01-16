#include "CoarsenerTests.h"
#include "AuxFunc.h"
#include "Random.h"

namespace Meso {
	/*
	void Test_Coarsener2(const Vector2i counts)
	{
		Field<bool, 2> finer_data{ Grid<2>(counts) };
		Vector2i finer_counts = finer_data.grid.Counts();

		int fnx = finer_counts[0], fny = finer_counts[1];
		for (int i = 0; i < fnx; i++) {
			real frac = (i + 0.0) / fnx;
			for (int j = 0; j < fny; j++) {
				bool value = (Random::Random() <= frac);
				finer_data(Vector2i(i, j)) = value;
			}
		}

		Vector2i coarse_counts = MathFunc::Round_Up_To_Align<2>(counts / 2, Grid<2>::Block_Size());
		Field<bool, 2> coarser_data{ Grid<2>(coarse_counts) };
		int cnx = coarse_counts[0], cny = coarse_counts[1];
		for (int i = 0; i < cnx; i++) {
			for (int j = 0; j < cny; j++) {
				//bool fixed = true;
				bool fixed = false;
				for (int t0 = 0; t0 < 2; t0++) {
					for (int t1 = 0; t1 < 2; t1++) {
						Vector2i finer_sub(i * 2 + t0, j * 2 + t1);
						//if (finer_data.grid.Valid(finer_sub) && !finer_data(finer_sub)) fixed = false;
						if (!finer_data.grid.Valid(finer_sub) || finer_data(finer_sub)) fixed = true;
					}
				}
				coarser_data(Vector2i(i, j)) = fixed;
			}
		}

		MaskedPoissonMapping<real, 2> finer, coarser;
		finer.Init(finer_data.grid);
		finer.fixed = finer_data;
		coarser.Init(coarser_data.grid);

		Coarsener<2>::Apply(coarser, finer);

		Field<bool, 2> coarsen_result; coarsen_result = coarser.fixed;

		if (ArrayFunc::Equals<bool>(coarsen_result.Data(), coarser_data.Data())) {
			Pass("Test_Coarsener2 passed {}", counts);
		}
		else {
			Error("Test_Coarsener2 failed {}", counts);
			Field<bool, 2> finer_temp = finer.fixed;
			Info("finer_data: \n{}", finer_data);
			Info("coarser_data:\n{}", coarser_data);
			Info("finer on device:\n{}", finer_temp);
			Info("coarsen_result:\n{}", coarsen_result);
			
		}
	}

	void Test_Coarsener3(const Vector3i counts)
	{
		Field<bool, 3> finer_data{ Grid<3>(counts) };
		Vector3i finer_counts = finer_data.grid.Counts();

		int fnx = finer_counts[0], fny = finer_counts[1], fnz = finer_counts[2];
		for (int i = 0; i < fnx; i++) {
			real frac = (i + 0.0) / fnx;
			for (int j = 0; j < fny; j++) {
				for (int k = 0; k < fnz; k++) {
					bool value = (Random::Random() <= frac);
					finer_data(Vector3i(i, j, k)) = value;
				}
			}
		}

		Vector3i coarser_counts = MathFunc::Round_Up_To_Align<3>(counts / 2, Grid<3>::Block_Size());
		Field<bool, 3> coarser_data{ Grid<3>(coarser_counts) };
		
		int cnx = coarser_counts[0], cny = coarser_counts[1], cnz = coarser_counts[2];
		for (int i = 0; i < cnx; i++) {
			for (int j = 0; j < cny; j++) {
				for (int k = 0; k < cnz; k++) {
					//bool fixed = true;
					bool fixed = false;
					for (int t0 = 0; t0 < 2; t0++) {
						for (int t1 = 0; t1 < 2; t1++) {
							for (int t2 = 0; t2 < 2; t2++) {
								Vector3i finer_sub(i * 2 + t0, j * 2 + t1, k * 2 + t2);
								//if (finer_data.grid.Valid(finer_sub) && !finer_data(finer_sub)) fixed = false;
								if (!finer_data.grid.Valid(finer_sub) || finer_data(finer_sub)) fixed = true;
							}
						}
					}
					coarser_data(Vector3i(i, j, k)) = fixed;
				}
			}
		}

		MaskedPoissonMapping<real, 3> finer, coarser;
		finer.Init(finer_data.grid);
		finer.fixed = finer_data;
		coarser.Init(coarser_data.grid);

		Coarsener<3>::Apply(coarser, finer);

		Field<bool, 3> coarsen_result; coarsen_result = coarser.fixed;

		if (ArrayFunc::Equals<bool>(coarsen_result.Data(), coarser_data.Data())) {
			Pass("Test_Coarsener3 passed {}", counts);
		}
		else {
			Error("Test_Coarsener3 failed {}", counts);
		}
	}
	*/
}
//////////////////////////////////////////////////////////////////////////
// Test restrictor
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Restrictor.h"
#include "Field.h"
#include "Random.h"
using namespace Meso;

template<class T>
T Weight_By_Offset(const int i) {
	if (i == -1 || i == 2) return 1.0 / 8;
	else if (i == 0 || i == 1) return 3.0 / 8;
}

template<class T, int d>
void Test_Restrictor(const Vector<int, d> counts) {
	Typedef_VectorD(d);
	Grid<d> finer_grid(counts);
	Grid<d> coarser_grid(counts / 2);
	Field<T, d> finer_host(finer_grid);
	Field<T, d> coarser_host(coarser_grid);
	Random::Fill_Random_Array<T>(finer_host.data);
	coarser_host.Calc_Cells(
		[&](const VectorDi coarser_coord) {
			T result = 0;
			if constexpr (d == 2) {
				for (int i = -1; i < 3; i++) {
					for (int j = -1; j < 3; j++) {
						VectorDi finer_coord = coarser_coord * 2 + Vector2i(i, j);
						real w = Weight_By_Offset<T>(i) * Weight_By_Offset<T>(j);
						if (finer_grid.Valid(finer_coord)) result += w * finer_host(finer_coord);
						else result += 0;
					}
				}
			}
			else if constexpr (d == 3) {
				for (int i = -1; i < 3; i++) {
					for (int j = -1; j < 3; j++) {
						for (int k = -1; k < 3; k++) {
							VectorDi finer_coord = coarser_coord * 2 + Vector3i(i, j, k);
							real w = Weight_By_Offset<T>(i) * Weight_By_Offset<T>(j) * Weight_By_Offset<T>(k);
							if (finer_grid.Valid(finer_coord)) result += w * finer_host(finer_coord);
							else result += 0;
						}
					}
				}
			}
			return result;
		}
	);
	
	Restrictor<T, d> restrictor;
	restrictor.Init(coarser_grid, finer_grid);
	FieldDv<T, d> finer_dev = finer_host;
	FieldDv<T, d> coarser_dev = coarser_host;
	restrictor.Apply(coarser_dev.data, finer_dev.data);
	Field<T, d> restrictor_result = coarser_dev;

	Info("coarser_host: \n{}\n", coarser_host);
	Info("restrictor result: \n{}\n", restrictor_result);

	if (!ArrayFunc::IsApprox<T>(restrictor_result.data, coarser_host.data)) {
		Error("Test_Restrictor failed for counts={}", counts);
	}
	else Pass("Test_Restrictor passed for counts={}", counts);
}
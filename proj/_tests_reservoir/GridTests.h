//////////////////////////////////////////////////////////////////////////
// Test grid structure
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
using namespace Meso;

template<class T, int d>
void Test_Grid_Index(const Vector<int, d> counts) {
	Grid<d> grid(counts);
	Field<int, d> visited(grid, -1);
	//test ind-cell and cell-ind
	int n = grid.DoF();
	for (int i = 0; i < n; i++) {
		Vector<int, d> cell = grid.Coord(i);
		int i1 = grid.Index(cell);
		Assert(i == i1, "Test_Grid_Index2 failed: cell {} and index {} doesn't match", cell, i);
		visited(cell) = 1;
	}
	if constexpr (d == 2) {
		for (int i = 0; i < counts[0]; i++) {
			for (int j = 0; j < counts[1]; j++) {
				Vector2i cell(i, j);
				Assert(visited(cell) == 1, "Test_Grid_Index2 failed: cell {} not visited", cell);
			}
		}
	}
	else if constexpr (d == 3) {
		for (int i = 0; i < counts[0]; i++) {
			for (int j = 0; j < counts[1]; j++) {
				for (int k = 0; k < counts[2]; k++) {
					Vector3i cell(i, j, k);
					Assert(visited(cell) == 1, "Test_Grid_Index2 failed: cell {} not visited", cell);
				}
			}
		}
	}

	FaceField<int, d> fvisited(grid, -1);
	for (int axis = 0; axis < d; axis++) {
		int fn = grid.Face_DoF(axis);
		for (int i = 0; i < fn; i++) {
			Vector<int, d> face = grid.Face_Coord(axis, i);
			int i1 = grid.Face_Index(axis, face);
			Assert(i == i1, "Test_Grid_Index2 failed: axis={}, face {} and index {} doesn't match", axis, face, i);
			fvisited(axis, face) = 1;
		}
		Vector<int, d> fcnt = grid.Face_Counts(axis);
		if constexpr (d == 2) {
			for (int i = 0; i < fcnt[0]; i++) {
				for (int j = 0; j < fcnt[1]; j++) {
					Vector2i face(i, j);
					Assert(fvisited(axis, face) == 1, "Test_Grid_Index2 failed: axis {} face {} not visited", axis, face);
				}
			}
		}
		else if constexpr (d == 3) {
			for (int i = 0; i < fcnt[0]; i++) {
				for (int j = 0; j < fcnt[1]; j++) {
					for (int k = 0; k < fcnt[2]; k++) {
						Vector3i face(i, j, k);
						Assert(fvisited(axis, face) == 1, "Test_Grid_Index2 failed: axis {} face {} not visited", axis, face);
					}
				}
			}
		}
	}


	Pass("Test_Grid_Index2 passed for counts={}", counts);
}
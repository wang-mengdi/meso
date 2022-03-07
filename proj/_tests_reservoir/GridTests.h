//////////////////////////////////////////////////////////////////////////
// Test grid structure
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
using namespace Meso;

template<class T>
void Test_Grid_Index2(const Vector2i counts) {
	Grid<2> grid(counts);
	Field<int, 2> visited(grid, -1);
	//test ind-cell and cell-ind
	int n = grid.DoF();
	for (int i = 0; i < n; i++) {
		Vector2i cell = grid.Coord(i);
		int i1 = grid.Index(cell);
		Assert(i == i1, "Test_Grid_Index2 failed: cell {} and index {} doesn't match", cell, i);
		visited(cell) = 1;
	}
	for (int i = 0; i < counts[0]; i++) {
		for (int j = 0; j < counts[1]; j++) {
			Vector2i cell(i, j);
			Assert(visited(cell) == 1, "Test_Grid_Index2 failed: cell {} not visited", cell);
		}
	}

	FaceField<int, 2> fvisited(grid, -1);
	for (int axis = 0; axis < 2; axis++) {
		int fn = grid.Face_DoF(axis);
		for (int i = 0; i < fn; i++) {
			Vector2i face = grid.Face_Coord(axis, i);
			int i1 = grid.Face_Index(axis, face);
			Assert(i == i1, "Test_Grid_Index2 failed: axis={}, face {} and index {} doesn't match", axis, face, i);
			fvisited(axis, face) = 1;
		}
		Vector<int, 2> fcnt = grid.counts; fcnt[axis]++;
		for (int i = 0; i < fcnt[0]; i++) {
			for (int j = 0; j < fcnt[1]; j++) {
				Vector2i face(i, j);
				Assert(fvisited(axis, face) == 1, "Test_Grid_Index2 failed: axis {} face {} not visited", axis, face);
			}
		}
	}


	Pass("Test_Grid_Index2 passed for counts={}", counts);
}
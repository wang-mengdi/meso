//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Interpolation.h"
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

namespace Meso {

	namespace VTKFunc {
		template<class T, int d, DataHolder side>
		void Output_VTS(const FaceField<T, d, side>& F, std::string file_name) {
			Typedef_VectorD(d);
			const auto grid = F.grid;
			int nx, ny, nz;
			if constexpr (d == 2) {
				nx = grid.counts[0];
				ny = grid.counts[1];
				nz = 1;
			}
			else {
				nx = grid.counts[0];
				ny = grid.counts[1];
				nz = grid.counts[2];
			}

			// setup VTK
			vtkXMLStructuredGridWriter* writer = vtkXMLStructuredGridWriter::New();
			vtkStructuredGrid* structured_grid = vtkStructuredGrid::New();
			structured_grid->SetDimensions(nx, ny, nz);
			vtkPoints* nodes = vtkPoints::New();
			nodes->Allocate(nx * ny * nz);
			//vtkDoubleArray* prsArray = vtkDoubleArray::New();
			//prsArray->SetNumberOfComponents(1);
			//prsArray->SetName("Pressure");
			vtkDoubleArray* velArray = vtkDoubleArray::New();
			velArray->SetName("Velocity");
			velArray->SetNumberOfComponents(3);

			FaceField<T, d> F_host = F;
			Field<VectorD, d> vf_host(grid);
			vf_host.Calc_Cells(
				[&](const VectorDi& cell) {
					VectorD pos = grid.Position(cell);
					return IntpLinear::Face_Vector(F_host, pos);
				}
			);

			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						VectorDi cell = VectorFunc::Vi<d>(i, j, k);
						VectorD pos = grid.Position(cell);
						Vector3 pos3 = VectorFunc::V<3>(pos);
						nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);

						//VectorD vec = Interpolation<PointIntpLinear>::Face_Vector<T, d, HOST>(field_host, pos);
						VectorD vec = vf_host(cell);
						Vector3 vec3 = VectorFunc::V<3>(vec);
						velArray->InsertNextTuple3(vec3[0], vec3[1], vec3[2]);
						//velArray->InsertNextTuple3(pos3[0], pos3[1], pos3[2]);
					}
				}
			}

			structured_grid->SetPoints(nodes);
			structured_grid->GetPointData()->AddArray(velArray);
			structured_grid->GetPointData()->SetActiveVectors("velocity");
			//structured_grid->GetPointData()->AddArray(prsArray);

#if (VTK_MAJOR_VERSION >=6)
			writer->SetInputData(structured_grid);
#else
			writer->SetInput(structured_grid);
#endif

			writer->SetFileName(file_name.c_str());
			writer->SetDataModeToAscii();
			writer->Write();

			structured_grid->Delete();
			writer->Delete();
			nodes->Delete();
			velArray->Delete();
		}
	}

}
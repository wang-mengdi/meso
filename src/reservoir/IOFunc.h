//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Mengdi Wang, Yitong Deng, Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Interpolation.h"
#include "Mesh.h"
#include <vtkNew.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
namespace Meso {

	namespace VTKFunc {
		template<class T, int d, DataHolder side>
		void Write_VTS(const FaceField<T, d, side>& F, std::string file_name) {
			Assert(!F.Empty(), "VTKFunc::Output_VTS error: empty FaceField");
			Typedef_VectorD(d);
			const auto grid = F.grid;
			int nx, ny, nz;
			VectorDi counts = grid.Counts();
			if constexpr (d == 2) {
				nx = counts[0];
				ny = counts[1];
				nz = 1;
			}
			else {
				nx = counts[0];
				ny = counts[1];
				nz = counts[2];
			}

			// setup VTK
			vtkNew<vtkXMLStructuredGridWriter> writer;
			vtkNew<vtkStructuredGrid> structured_grid;
			structured_grid->SetDimensions(nx, ny, nz);
			vtkNew<vtkPoints> nodes;
			nodes->Allocate(nx * ny * nz);
			//vtkDoubleArray* prsArray = vtkDoubleArray::New();
			//prsArray->SetNumberOfComponents(1);
			//prsArray->SetName("Pressure");
			vtkNew<vtkDoubleArray> velArray;
			velArray->SetName("Velocity");
			velArray->SetNumberOfComponents(3);

			FaceField<T, d> F_host = F;
			Field<VectorD, d> vf_host(grid);
			vf_host.Calc_Nodes(
				[&](const VectorDi& cell) {
					VectorD pos = grid.Position(cell);
					return IntpLinear::Face_Vector(F_host, pos);
				}
			);

			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						VectorDi cell = MathFunc::Vi<d>(i, j, k);
						VectorD pos = grid.Position(cell);
						Vector3 pos3 = MathFunc::V<3>(pos);
						nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);

						//VectorD vec = Interpolation<PointIntpLinear>::Face_Vector<T, d, HOST>(field_host, pos);
						VectorD vec = vf_host(cell);
						Vector3 vec3 = MathFunc::V<3>(vec);
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
			writer->SetDataModeToBinary();
			writer->Write();
		}

		template<int d>
		void Write_VTU_Particles(const Array<Vector<real, d>>& pos, const Array<Vector<real, d>>& vel, std::string file_name) {
			//output particles
			Typedef_VectorD(d);
			// setup VTK
			vtkNew<vtkXMLUnstructuredGridWriter> writer;
			vtkNew<vtkUnstructuredGrid> unstructured_grid;

			vtkNew<vtkPoints> nodes;
			nodes->Allocate(pos.size());
			vtkNew<vtkDoubleArray> posArray;
			posArray->SetName("Position");
			posArray->SetNumberOfComponents(d);
			vtkNew<vtkDoubleArray> velArray;
			velArray->SetName("Velocity");
			velArray->SetNumberOfComponents(d);

			for (int i = 0; i < pos.size(); i++) {
				Vector3 pos3 = MathFunc::V<3>(pos[i]);
				nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);

				Vector3 vel3 = MathFunc::V<3>(vel[i]);
				velArray->InsertNextTuple3(vel3[0], vel3[1], vel3[2]);
			}

			unstructured_grid->SetPoints(nodes);

			unstructured_grid->GetPointData()->AddArray(velArray);
			unstructured_grid->GetPointData()->SetActiveVectors("velocity");

			writer->SetFileName(file_name.c_str());
			writer->SetInputData(unstructured_grid);
			writer->Write();
		}
	}


	namespace OBJFunc {

		// reference: https://github.com/tinyobjloader/tinyobjloader#example-code-new-object-oriented-api
		// `mesh->indice->normal` represents the normal on vetices, which is used for rendering. 
		// Therefore, we only load vertices and faces from .obj file, then compute normals on each faces.


		//template<class T, int d, int ed>
		//bool Read_Obj(const std::string& filename, const VertexMatrix<T, d>& vertex_matrix, const ElementMatrix<ed>& element_matrix) {
		//	igl::readOBJ(filename, vertex_matrix, element_matrix);
		//}


		template<class T, int d, int ed>
		bool Write_OBJ(const std::string& filename, const VertexMatrix<T, d>& vertex_matrix, const ElementMatrix<ed>& element_matrix) {
			if constexpr (d == 2) {
				Assert(ed == 2, "Vertex have dim={}, but element has dim={}", d, ed);
				FILE* fp = fopen(filename.c_str(), "w");
				Assert(fp != nullptr, "Failed to open file [ {} ] for write.\n", filename);

				for (size_t i = 0; i < vertex_matrix.rows(); i++)
					fprintf(fp, "v %f %f\n", vertex_matrix.row(i)[0], vertex_matrix.row(i)[1]);

				fprintf(fp, "\n");
				for (size_t i = 0; i < element_matrix.rows(); i++)
					fprintf(fp, "l %d %d\n", element_matrix.row(i)[0] + 1, element_matrix.row(i)[1] + 1);

				fclose(fp);
				return true;
			}
			else if constexpr (d == 3) {
				Assert(ed > 2, "Vertex have dim={}, but element has dim={}", d, ed);
				return igl::writeOBJ<VertexMatrix<T, d>, ElementMatrix<ed>>(filename, vertex_matrix, element_matrix);
			}
		}
	}


}
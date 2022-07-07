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
		void VTS_Add_Field(vtkStructuredGrid& vtk_grid, const Field<T, d, side>& F, const std::string name) {
			Typedef_VectorD(d);
			const auto grid = F.grid;
			int nx, ny, nz;
			VectorDi counts = grid.Counts();
			if constexpr (d == 2) nx = counts[0], ny = counts[1], nz = 1;
			else nx = counts[0], ny = counts[1], nz = counts[2];

			constexpr bool is_vec3 = std::is_same<T, Vector3>::value;

			vtkNew<vtkDoubleArray> f_array;
			f_array->SetName(name.c_str());
			if (is_vec3) f_array->SetNumberOfComponents(3);
			else f_array->SetNumberOfComponents(1);

			Field<T, d> F_host;
			if constexpr (side == HOST) F_host.Shallow_Copy(F);
			else F_host.Deep_Copy(F);

			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						VectorDi cell = MathFunc::Vi<d>(i, j, k);
						if constexpr (std::is_same<T, Vector3>::value) {
							Vector3 vec3 = F_host(cell);
							f_array->InsertNextTuple3(vec3[0], vec3[1], vec3[2]);
						}
						else {
							f_array->InsertNextValue(F_host(cell));
							//f_array->InsertComponent(F_host(cell));
						}
					}
				}
			}

			vtk_grid.GetPointData()->AddArray(f_array);
			if (is_vec3) vtk_grid.GetPointData()->SetActiveVectors(name.c_str());
			else vtk_grid.GetPointData()->SetActiveScalars(name.c_str());
		}

		template<class T, int d, DataHolder side>
		void VTS_Add_Vector(vtkStructuredGrid& vtk_grid, const FaceField<T, d, side>& F, const std::string name) {
			//vtk_grid must be already allocated and initialized with points
			Typedef_VectorD(d);
			const auto grid = F.grid;

			FaceField<T, d> F_host;
			if constexpr (side == HOST) F_host.Shallow_Copy(F);
			else F_host.Deep_Copy(F);

			Field<Vector3, d> vf_host(grid);
			vf_host.Calc_Nodes(
				[&](const VectorDi& cell) {
					VectorD pos = grid.Position(cell);
					Vector<real, d> vec = IntpLinear::Face_Vector(F_host, pos).template cast<real>();
					return MathFunc::V<3>(vec);
				}
			);

			VTS_Add_Field(vtk_grid, vf_host, name);
		}

		template<int d>
		void VTS_Init_Grid(vtkStructuredGrid& vtk_grid, const Grid<d> grid) {
			Typedef_VectorD(d);
			int nx, ny, nz;
			VectorDi counts = grid.Counts();
			if constexpr (d == 2) nx = counts[0], ny = counts[1], nz = 1;
			else nx = counts[0], ny = counts[1], nz = counts[2];

			vtk_grid.SetDimensions(nx, ny, nz);
			vtkNew<vtkPoints> nodes;
			nodes->Allocate(nx * ny * nz);

			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						VectorDi cell = MathFunc::Vi<d>(i, j, k);
						VectorD pos = grid.Position(cell);
						Vector3 pos3 = MathFunc::V<3>(pos);
						nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
					}
				}
			}

			vtk_grid.SetPoints(nodes);
		}

		template<class T, int d, DataHolder side>
		void Write_VTS(const FaceField<T, d, side>& F, std::string file_name) {
			Assert(!F.Empty(), "VTKFunc::Output_VTS error: empty FaceField");
			Typedef_VectorD(d);

			// setup VTK
			vtkNew<vtkXMLStructuredGridWriter> writer;
			vtkNew<vtkStructuredGrid> structured_grid;
			VTS_Init_Grid(*structured_grid, F.grid);
			VTS_Add_Vector(*structured_grid, F, "velocity");

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
			velArray->SetNumberOfComponents(3);

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
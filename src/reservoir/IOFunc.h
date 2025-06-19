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
#include <vtkXMLPolyDataWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkLine.h>

#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include "tiny_obj_loader.h"

namespace Meso {

	namespace VTKFunc {
		template<class T, int d, DataHolder side>
		void VTS_Add_Field(vtkStructuredGrid& vtk_grid, const Field<T, d, side>& F, const std::string name, bool cell_data=true/*point_data otherwise*/) {
			Typedef_VectorD(d);
			const auto grid = F.grid;
			int nx, ny, nz;
			VectorDi counts = grid.Counts();
			if constexpr (d == 2) nx = counts[0], ny = counts[1], nz = 1;
			else nx = counts[0], ny = counts[1], nz = counts[2];

			constexpr bool is_vec3 = std::is_same<T, Vector3>::value;
			constexpr bool is_vec2 = std::is_same<T, Vector2>::value;

			vtkNew<vtkDoubleArray> f_array;
			f_array->SetName(name.c_str());
			if (is_vec3||is_vec2) f_array->SetNumberOfComponents(3);
			else f_array->SetNumberOfComponents(1);

			Field<T, d> F_host;
			if constexpr (side == HOST) F_host.Shallow_Copy(F);
			else F_host.Deep_Copy(F);

			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						VectorDi cell = MathFunc::Vi<d>(i, j, k);
						if constexpr (is_vec2) {
							Vector2 vec2 = F_host(cell);
							f_array->InsertNextTuple3(vec2[0], vec2[1],0);
						}else if constexpr (is_vec3) {
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
			
			if(!cell_data){
				vtkPointData* point_data = vtk_grid.GetPointData();
				point_data->AddArray(f_array);
			}
			else {
				vtkCellData* cell_data = vtk_grid.GetCellData();
				cell_data->AddArray(f_array);
			}
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
					Vector<real, d> vec = IntpLinearClamp::Face_Vector(F_host, pos).template cast<real>();
					return MathFunc::V<3>(vec);
				}
			);

			VTS_Add_Field(vtk_grid, vf_host, name, false);
		}

		template<class T, int d, DataHolder side>
		void VTS_Add_Colocated_Vector(vtkStructuredGrid& vtk_grid, const ArrayF<Field<T, d, side>, d>& F, const std::string name) {
			//vtk_grid must be already allocated and initialized with points
			Typedef_VectorD(d);
			const auto grid = F[0].grid;

			ArrayF<Field<T, d>, d> F_host;
			if constexpr (side == HOST)
				for (int axis = 0; axis < d; axis++)
					F_host[axis].Shallow_Copy(F[axis]);
			else
				for (int axis = 0; axis < d; axis++)
					F_host[axis].Deep_Copy(F[axis]);

			Field<Vector3, d> vf_host(grid);
			vf_host.Calc_Nodes(
				[&](const VectorDi& cell) {
					int index = grid.Index(cell);
					VectorD vec;
					for (int axis = 0; axis < d; axis++)
						vec[axis] = F_host[axis](cell);
					return MathFunc::V<3>(vec);
				}
			);

			VTS_Add_Field(vtk_grid, vf_host, name, false);
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
		template<class T, int d, DataHolder side> void Write_VTS(const FaceField<T, d, side>& F, const fs::path& path) { Write_VTS<T, d, side>(F, path.string()); }

		template<class T, int d, DataHolder side>
		void Write_VTS(const Field<T, d, side>& F, std::string file_name) {
			Assert(!F.Empty(), "VTKFunc::Output_VTS error: empty Field");
			Typedef_VectorD(d);

			// setup VTK
			vtkNew<vtkXMLStructuredGridWriter> writer;
			vtkNew<vtkStructuredGrid> structured_grid;
			VTS_Init_Grid(*structured_grid, F.grid);
			VTS_Add_Field(*structured_grid, F, "scalar", false);

#if (VTK_MAJOR_VERSION >=6)
			writer->SetInputData(structured_grid);
#else
			writer->SetInput(structured_grid);
#endif

			writer->SetFileName(file_name.c_str());
			writer->SetDataModeToBinary();
			writer->Write();
		}
		template<class T, int d, DataHolder side> void Write_VTS(const Field<T, d, side>& F, const fs::path& path) { Write_VTS<T, d, side>(F, path.string()); }

		template<class T, int d, DataHolder side>
		void Write_VTS(const ArrayF<Field<T, d, side>, d>& F, std::string file_name) {
			Assert(!F[0].Empty(), "VTKFunc::Output_VTS error: empty Colocated Vector Field");
			Typedef_VectorD(d);

			// setup VTK
			vtkNew<vtkXMLStructuredGridWriter> writer;
			vtkNew<vtkStructuredGrid> structured_grid;
			VTS_Init_Grid(*structured_grid, F[0].grid);
			VTS_Add_Colocated_Vector(*structured_grid, F, "velocity");

#if (VTK_MAJOR_VERSION >=6)
			writer->SetInputData(structured_grid);
#else
			writer->SetInput(structured_grid);
#endif

			writer->SetFileName(file_name.c_str());
			writer->SetDataModeToBinary();
			writer->Write();
		}

		void Write_VTS(vtkStructuredGrid& vtk_grid, std::string file_name);

		template<class T, int d,int ed>
		void VTP_Init_Mesh(vtkPolyData& vtp, const VertexMatrix<T, d>& vertices, const ElementMatrix<ed>& elements) {
			if (elements.rows() == 0) { Warn("VTKFunc::VTP_Init_Mesh: no vertices to output"); return; }
			Typedef_VectorD(d);
			vtkNew<vtkPoints> points;
			for (int i = 0; i < vertices.rows();i++)
			{
				Vector3 vec=MathFunc::V<3>(VectorD(vertices.row(i)));
				points->InsertNextPoint(vec[0], vec[1], vec[2]);
			}
			vtp.SetPoints(points);

			if constexpr(d==2){
				vtkNew<vtkCellArray> lines;
				for (int i = 0; i < elements.rows(); i++) {
					vtkNew<vtkLine> line;
					line->GetPointIds()->SetId(0, elements.row(i)[0]);
					line->GetPointIds()->SetId(1, elements.row(i)[1]);
					lines->InsertNextCell(line);
				}
				vtp.SetLines(lines);
			}
		}

		void Write_VTP(vtkPolyData& vtp, std::string file_name);

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

		template<int d>
		void Write_VTU_Particles(const Array<Vector<real, d>>& pos, const Array<Vector<real, d>>& vel, const fs::path& path) { Write_VTU_Particles<d>(pos, vel, path.string()); }

		template<int d>
		void Write_VTU_Particles(const Array<Vector<real, d>>& pos, const fs::path& path) {
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

			for (int i = 0; i < pos.size(); i++) {
				Vector3 pos3 = MathFunc::V<3>(pos[i]);
				nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
			}

			unstructured_grid->SetPoints(nodes);

			writer->SetFileName(path.string().c_str());
			writer->SetInputData(unstructured_grid);
			writer->Write();
		}

		template<class T, int d>
		void Output_Mesh_As_VTU(const VertexMatrix<T, d>& verts, const ElementMatrix<d>& elements, const std::string file_name) {
			Typedef_VectorD(d);

			// setup VTK
			vtkNew<vtkXMLUnstructuredGridWriter> writer;
			vtkNew<vtkUnstructuredGrid> unstructured_grid;

			vtkNew<vtkPoints> nodes;
			nodes->Allocate(verts.rows());
			vtkNew<vtkCellArray> cellArray;

			for (int i = 0; i < verts.rows(); i++) {
				Vector<real, d> pos = verts.row(i).template cast<real>();
				Vector3 pos3 = MathFunc::V<3>(pos);
				nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
			}
			unstructured_grid->SetPoints(nodes);

			if constexpr (d == 2) {
				for (int i = 0; i < elements.rows(); i++) {
					vtkNew<vtkLine> line;
					line->GetPointIds()->SetId(0, elements(i, 0));
					line->GetPointIds()->SetId(1, elements(i, 1));
					cellArray->InsertNextCell(line);
				}
				unstructured_grid->SetCells(VTK_LINE, cellArray);
			}
			else if constexpr (d == 3) {
				for (int i = 0; i < elements.rows(); i++) {
					vtkNew<vtkTriangle> triangle;
					triangle->GetPointIds()->SetId(0, elements(i, 0));
					triangle->GetPointIds()->SetId(1, elements(i, 1));
					triangle->GetPointIds()->SetId(2, elements(i, 2));
					cellArray->InsertNextCell(triangle);
				}
				unstructured_grid->SetCells(VTK_TRIANGLE, cellArray);
			}

			writer->SetFileName(file_name.c_str());
			writer->SetInputData(unstructured_grid);
			writer->Write();
		}

		template<class T, int d> void Output_Mesh_As_VTU(const VertexMatrix<T, d>& verts, const ElementMatrix<d>& elements, const fs::path& path) { Output_Mesh_As_VTU<T, d>(verts, elements, path.string()); }
	}


	namespace OBJFunc {

		// reference: https://github.com/tinyobjloader/tinyobjloader#example-code-new-object-oriented-api
		// `mesh->indice->normal` represents the normal on vetices, which is used for rendering. 
		// Therefore, we only load vertices and faces from .obj file, then compute normals on each faces.

		template<class T, int d, int ed>
		bool Read_Obj(const std::string& filename, VertexMatrix<T, d>& vertex_matrix, ElementMatrix<ed>& element_matrix) {
			MatrixXd vertex_matrix_tmp;
			bool result;
			if constexpr (d == 2 && ed==2) {
				tinyobj::ObjReaderConfig reader_config;

				tinyobj::ObjReader reader;
				if (!reader.ParseFromFile(filename , reader_config)) {
					if (!reader.Error().empty()) {
						Error("TinyObjReader: {}", reader.Error());
					}
				}
				auto& attrib = reader.GetAttrib();
				auto& shapes = reader.GetShapes(); //assume that there is only one shape
				Assert(shapes.size() == 1, "Read_Obj: Does support not multiple shapes {}",shapes.size());
				int line_num = shapes[0].lines.num_line_vertices.size();
				vertex_matrix.resize(attrib.GetVertices().size()/3, 2);
				element_matrix.resize(line_num, 2);
				size_t index_offset = 0;
				for (size_t l = 0; l < line_num; l++) {
					size_t lv = size_t(shapes[0].lines.num_line_vertices[l]);
					Assert(lv == 2); //we know that each segment only has two vertices
					// Loop over vertices in the face.
					for (size_t v = 0; v < lv; v++) {
						// access to vertex
						tinyobj::index_t idx = shapes[0].lines.indices[index_offset+v];
						element_matrix(l , v) = idx.vertex_index;
					}
					index_offset += lv;
				}
				
				//there are three dimensions each point
				for (size_t v = 0; v < attrib.vertices.size()/3; v++) {
					vertex_matrix(v, 0) = attrib.vertices[v*3];
					vertex_matrix(v, 1) = attrib.vertices[v*3 + 1];
				}
			}
			else if constexpr (d==3 && ed==3) {
				result = igl::readOBJ(filename, vertex_matrix, element_matrix);
			}
			else {
				Error("Demension with d={} and ed={} is not supported yet.", d, ed);
			}
			return result;
		}

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
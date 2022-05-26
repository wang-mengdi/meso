//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Mengdi Wang, Yitong Deng
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

#define TINYOBJLOADER_IMPLEMENTATION   
#include <tiny_obj_loader.h>

namespace Meso {

	namespace VTKFunc {
		template<class T, int d, DataHolder side>
		void Write_VTS(const FaceField<T, d, side>& F, std::string file_name) {
			Assert(!F.Empty(), "VTKFunc::Output_VTS error: empty FaceField");
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

		// Using tinyobj 1.0.7 
		// ONLY support traiangle mesh now

		template<class T>
		void Read_Meshes(const std::string& file_name, Array<std::shared_ptr<T>>& meshes) {

			tinyobj::attrib_t attrib;
			std::vector<tinyobj::shape_t> shapes; std::vector<tinyobj::material_t> materials;
			std::string warn;std::string err;
			bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, file_name.c_str());
			if (!warn.empty()) std::cout << warn << std::endl;
			if (!err.empty()) std::cerr << err << std::endl;
			if (!ret) exit(1);

			meshes.resize(shapes.size());
			size_t v_begin = 0; size_t v_end = 0;

			for (size_t s = 0; s < shapes.size(); s++) {
				meshes[s] = std::make_shared<T>();

				const tinyobj::mesh_t& mesh = shapes[s].mesh;
				// three vertices form a triagnle face 
				for (size_t i = 0; i < mesh.indices.size() / 3; i++) {
					meshes[s]->faces.push_back(Vector3i(
						mesh.indices[i * 3 + 0].vertex_index - v_begin,
						mesh.indices[i * 3 + 1].vertex_index - v_begin,
						mesh.indices[i * 3 + 2].vertex_index - v_begin));

					v_end = std::max({ v_end, v_begin + (size_t)meshes[s]->faces.back().maxCoeff() });
				}
				for (size_t i = v_begin; i <= v_end; i++) {
					(*meshes[s]->vertices).push_back(Vector3(attrib.vertices[i * 3 + 0], attrib.vertices[i * 3 + 1], attrib.vertices[i * 3 + 2]));
				}
				v_begin = v_end + 1;
			}
		}

		template<class T>
		void Read_Mesh(const std::string& file_name, std::shared_ptr<T>& mesh) {
			Array<std::shared_ptr<T>> meshes;
			Read_Meshes<T>(file_name, meshes);
			Assert(meshes.size() == 1, "Wrong mesh number {} in {}, which should be 1", meshes.size(), file_name);
			mesh = meshes[0];
		}

		bool Write_Obj(const std::string& filename, const tinyobj::attrib_t& attributes, const std::vector<tinyobj::shape_t>& shapes, const std::vector<tinyobj::material_t>& materials) {
			FILE* fp = fopen(filename.c_str(), "w");
			if (!fp) {
				fprintf(stderr, "Failed to open file [ %s ] for write.\n", filename.c_str());
				return false;
			}

			for (size_t k = 0; k < attributes.vertices.size(); k += 3)
				fprintf(fp, "v %f %f %f\n", attributes.vertices[k + 0], attributes.vertices[k + 1], attributes.vertices[k + 2]);

			fprintf(fp, "\n");

			for (size_t i = 0; i < shapes.size(); i++) {
				fprintf(fp, "\n");
				fprintf(fp, "g Unknown\n");

				int face_index = 0;
				for (size_t k = 0; k < shapes[i].mesh.indices.size(); k += shapes[i].mesh.num_face_vertices[face_index++]) {
					unsigned char v_per_f = shapes[i].mesh.num_face_vertices[face_index];
					fprintf(fp, "f");
					for (int l = 0; l < v_per_f; l++) {
						const tinyobj::index_t& ref = shapes[i].mesh.indices[k + l];
						fprintf(fp, " %d", ref.vertex_index + 1);
					}
					fprintf(fp, "\n");
				}
			}

			fclose(fp);
			return true;
		}

		template<class T>
		bool Write_Meshes(const std::string& filename, const Array<std::shared_ptr<T>>& meshes) {
			tinyobj::attrib_t attribute;
			std::vector<tinyobj::shape_t> shapes(meshes.size());

			for (size_t m = 0; m < meshes.size(); m++) {
				auto& mesh = meshes[m];
				int offset = attribute.vertices.size();

				for (auto& vertex : mesh->Vertices()) {
					attribute.vertices.push_back(vertex[0]);
					attribute.vertices.push_back(vertex[1]);
					attribute.vertices.push_back(vertex[2]);
				}

				shapes[m] = tinyobj::shape_t();
				for (size_t f = 0; f < mesh->Elements().size(); f++) {
					shapes[m].mesh.num_face_vertices.push_back(3);
					shapes[m].mesh.indices.push_back({ mesh->Elements()[f][0] + offset / 3, -1, -1 });
					shapes[m].mesh.indices.push_back({ mesh->Elements()[f][1] + offset / 3, -1, -1 });
					shapes[m].mesh.indices.push_back({ mesh->Elements()[f][2] + offset / 3, -1, -1 });
				}
			}

			return Write_Obj(filename, attribute, shapes, std::vector<tinyobj::material_t>());
		}

		template<class T>
		bool Write_Mesh(const std::string& filename, const std::shared_ptr<T> mesh) {
			Array<std::shared_ptr<T>> meshes; meshes.push_back(mesh);
			return Write_Meshes<T>(filename, meshes);
		}

	}


}
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

#define TINYOBJLOADER_IMPLEMENTATION   
#include "tiny_obj_loader.h"

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


	namespace OBJFunc {

		// reference: https://github.com/tinyobjloader/tinyobjloader#example-code-new-object-oriented-api
		// `mesh->indice->normal` represents the normal on vetices, which is used for rendering. 
		// Therefore, we only load vertices and faces from .obj file, then compute normals on each faces.
		
		// Idea: we can save face normal in .obj file correponding to each face, where len(vn) == len(f)

		Array<std::shared_ptr<TriangleMesh<3>>> Read_Mesh_From_Obj_File(const std::string& file_name) {
			Array<std::shared_ptr<TriangleMesh<3>>> meshes;

			tinyobj::ObjReaderConfig reader_config;
			tinyobj::ObjReader reader;

			if (!reader.ParseFromFile(file_name, reader_config)) {
				if (!reader.Error().empty()) std::cerr << "TinyObjReader: " << reader.Error();
				exit(1);
			}
			if (!reader.Warning().empty()) std::cout << "TinyObjReader: " << reader.Warning();


			const tinyobj::attrib_t& attrib = reader.GetAttrib();
			const std::vector<tinyobj::shape_t>& shapes = reader.GetShapes();
			meshes.resize(shapes.size());

			size_t v_begin = 0; size_t v_end = 0;

			for (size_t s = 0; s < shapes.size(); s++) {
				meshes[s] = std::make_shared<TriangleMesh<3>>();

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
			return meshes;
		}

		

		static std::string GetFileBasename(const std::string& FileName)
		{
			if (FileName.find_last_of(".") != std::string::npos)
				return FileName.substr(0, FileName.find_last_of("."));
			return "";
		}
		bool WriteObj(const std::string& filename, const tinyobj::attrib_t& attributes, const std::vector<tinyobj::shape_t>& shapes, const std::vector<tinyobj::material_t>& materials) {
			FILE* fp = fopen(filename.c_str(), "w");
			if (!fp) {
				fprintf(stderr, "Failed to open file [ %s ] for write.\n", filename.c_str());
				return false;
			}

			std::string basename = GetFileBasename(filename);
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
		bool Write_Mesh_To_Obj_File(const std::string& filename, Array<std::shared_ptr<TriangleMesh<3>>>& meshes) {
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
				for (size_t f = 0; f < mesh->Faces().size(); f++) {
					shapes[m].mesh.num_face_vertices.push_back(3);
					shapes[m].mesh.indices.push_back({ mesh->Faces()[f][0] + offset / 3, -1, -1 });
					shapes[m].mesh.indices.push_back({ mesh->Faces()[f][1] + offset / 3, -1, -1 });
					shapes[m].mesh.indices.push_back({ mesh->Faces()[f][2] + offset / 3, -1, -1 });
				}
			}

			return OBJFunc::WriteObj(filename, attribute, shapes, std::vector<tinyobj::material_t>());
		}}
}
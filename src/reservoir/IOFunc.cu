//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#define TINYOBJLOADER_IMPLEMENTATION   
#include "IOFunc.h"

namespace Meso {
	namespace OBJFunc {
		//bool Write_Obj(const std::string& filename, const tinyobj::attrib_t& attributes, const std::vector<tinyobj::shape_t>& shapes, const std::vector<tinyobj::material_t>& materials) {
			/*FILE* fp = fopen(filename.c_str(), "w");
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
		}*/
	}
}
#include "WaveFunc.h"
namespace Meso {
	namespace WaveFunc {
		__host__ __device__ C Vector2C_Dot(const Vector2C& a, const Vector2C& b)
		{
			return a[0] * b[0] + a[1] * b[1];
		}

		__host__ __device__ Vector4 Vector2C_To_Vector4(const Vector<C, 2>& v)
		{
			return Vector4(v[0].real(), v[0].imag(), v[1].real(), v[1].imag());
		}

		__host__ __device__ Vector<C, 2> Vector4_To_Vector2C(const Vector4& v)
		{
			Vector<C, 2> c(C(v[0], v[1]), C(v[2], v[3])); return c;
		}
		
		__global__ void W2V_Mapping_Kernel2_Padding0(const Grid<2> grid, real* face_x, real* face_y, const Vector2C* cell, const real h_bar_over_dx)
		{
			Typedef_VectorD(2);
			VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);

			if (!grid.Valid(coord)) return;
			const Vector2C& cell_data = cell[grid.Index(coord)];
			const Vector2C cell_data_conj(thrust::conj(cell_data[0]), thrust::conj(cell_data[1]));


			int face_ind_x = Neighbor_Face_Idx(grid, coord, 0, 1);
			int face_ind_y = Neighbor_Face_Idx(grid, coord, 1, 1);

			if (face_ind_x != -1 && grid.Valid(coord + VectorDi::Unit(0))) {
				const Vector2C& nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(0))];
				const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
				face_x[face_ind_x] = h_bar_over_dx * thrust::arg(re_product);
			}
			if (face_ind_y != -1 && grid.Valid(coord + VectorDi::Unit(1))) {
				const Vector2C& nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(1))];
				const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
				face_y[face_ind_y] = h_bar_over_dx * thrust::arg(re_product);
			}
		}

		// for blockDim = (4, 4, 4)
		__global__ void W2V_Mapping_Kernel3_Padding0(const GridIndexer<3> grid, real* face_x, real* face_y, real* face_z, const Vector2C* cell, const real h_bar_over_dx)
		{
			Typedef_VectorD(3);
			VectorDi coord = GPUFunc::Thread_Coord<3>(blockIdx, threadIdx);

			const int idx = threadIdx.x;
			const int idy = threadIdx.y;
			const int idz = threadIdx.z;

			if (!grid.Valid(coord)) return;
			const Vector2C cell_data = cell[grid.Index(coord)];
			const Vector2C cell_data_conj(thrust::conj(cell_data[0]), thrust::conj(cell_data[1]));

			int face_ind_x = Neighbor_Face_Idx(grid, coord, 0, 1);
			int face_ind_y = Neighbor_Face_Idx(grid, coord, 1, 1);
			int face_ind_z = Neighbor_Face_Idx(grid, coord, 2, 1);

			if (face_ind_x != -1 && grid.Valid(coord + VectorDi::Unit(0))) {
				const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(0))];
				const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
				face_x[face_ind_x] = h_bar_over_dx * thrust::arg(re_product);

			}
			if (face_ind_y != -1 && grid.Valid(coord + VectorDi::Unit(1))) {
				const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(1))];
				const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
				face_y[face_ind_y] = h_bar_over_dx * thrust::arg(re_product);
			}
			if (face_ind_z != -1 && grid.Valid(coord + VectorDi::Unit(2))) {
				const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(2))];
				const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
				face_z[face_ind_z] = h_bar_over_dx * thrust::arg(re_product);
			}
		}
	}
}
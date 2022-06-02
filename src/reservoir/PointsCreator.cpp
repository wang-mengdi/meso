#include "PointsCreator.h"

namespace Meso {
	using namespace CommonConstants;

	void Initialize_Lattice_Points(const Vector3& origin, const int nx, const int ny, const Vector3 dx, const Vector3 dy, Points& pts, Array<Vector3>& pos) {
		pts.Resize((nx + 1) * (ny + 1));
		int idx = 0;
		for (int i = 0; i <= nx; i++) {
			for (int j = 0; j <= ny; j++) {
				Vector3 curr_pos = origin + i * dx + j * dy;
				pos[idx] = curr_pos;
				idx++;
			}
		}
	}

	real Initialize_Box_Points_2D(const Vector2& origin, const Vector2i& size, const Vector2& _dx, Points& pts, Array<Vector2>& pos, bool keep_existing) {
		int nx = size[0]; int ny = size[1]; real dx = _dx[0]; real dy = _dx[1];
		pts.Resize(nx * ny);
		int idx = 0;
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				Vector2 curr_pos = origin + (real)i * dx * Vector2::Unit(0) + (real)j * dy * Vector2::Unit(1);
				pos[idx] = curr_pos;
				idx++;
			}
		}
		return (dx + dy) / 2.;
	}

	real Initialize_Box_Rim_Points_2D(const Vector2& origin, const Vector2i& pad_size, const Vector2i& int_size, const Vector2& _dx, Points& pts, Array<Vector2>& pos, bool keep_existing) {
		int nx_in = int_size[0]; int ny_in = int_size[1];
		int nx_out = 2 * pad_size[0] + nx_in; int ny_out = 2 * pad_size[1] + ny_in; 
		real dx = _dx[0]; real dy = _dx[1];
		int begin = 0;
		if (keep_existing) begin = pts.Size();
		pts.Resize(begin + nx_out * ny_out - nx_in * ny_in);
		int idx = begin;
		for (int i = -pad_size[0]; i < nx_in + pad_size[0]; i++) {
			for (int j = -pad_size[1]; j < ny_in + pad_size[1]; j++) {
				if ((i < 0) || (i >= nx_in) || (j < 0) || (j >= ny_in)) {
					Vector2 curr_pos = origin + (real)i * dx * Vector2::Unit(0) + (real)j * dy * Vector2::Unit(1);
					pos[idx] = curr_pos;
					idx++;
				}
			}
		}
		return (dx + dy) / 2.;
	}

	real Initialize_Box_Points_3D(const Vector3& origin, const Vector3i& size, const Vector3& _dx, Points& pts, Array<Vector3>& pos, bool keep_existing) {
		int nx = size[0]; int ny = size[1]; int nz = size[2];
		real dx = _dx[0]; real dy = _dx[1]; real dz = _dx[2];
		pts.Resize(nx * ny * nz);
		int idx = 0;
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nz; k++) {
					Vector3 curr_pos = origin + (real)i * dx * Vector3::Unit(0) + (real)j * dy * Vector3::Unit(1) + (real)k * dz * Vector3::Unit(2);
					pos[idx] = curr_pos;
					idx++;
				}
			}
		}
		return (dx + dy + dz) / 3.;
	}

	real Initialize_Box_Rim_Points_3D(const Vector3& origin, const Vector3i& pad_size, const Vector3i& int_size, const Vector3& _dx, Points& pts, Array<Vector3>& pos, bool keep_existing) {
		int nx_in = int_size[0]; int ny_in = int_size[1]; int nz_in = int_size[2];
		int nx_out = 2 * pad_size[0] + nx_in; int ny_out = 2 * pad_size[1] + ny_in; int nz_out = 2 * pad_size[2] + nz_in;
		real dx = _dx[0]; real dy = _dx[1]; real dz = _dx[2];
		int begin = 0;
		if (keep_existing) begin = pts.Size();
		int new_size = begin + nx_out * ny_out * nz_out - nx_in * ny_in * nz_in;
		pts.Resize(begin + nx_out * ny_out * nz_out - nx_in * ny_in * nz_in);
		int idx = begin;
		for (int i = -pad_size[0]; i < nx_in + pad_size[0]; i++) {
			for (int j = -pad_size[1]; j < ny_in + pad_size[1]; j++) {
				for (int k = -pad_size[2]; k < nz_in + pad_size[2]; k++) {
					if ((i < 0) || (i >= nx_in) || (j < 0) || (j >= ny_in) || (k < 0) || (k >= nz_in)) {
						Vector3 curr_pos = origin + (real)i * dx * Vector3::Unit(0) + (real)j * dy * Vector3::Unit(1) + (real)k * dz * Vector3::Unit(2);
						pos[idx] = curr_pos;
						idx++;
					}
				}
			}
		}
		return (dx + dy + dz) / 3.;
	}


	// algorithm: https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
	// note that the num_pts you put in is only an average, it is not necessarily the num of pts you get
	real Initialize_Sphere_Points_Regular(const Vector3& center, const real R, const int num_pts, Points& pts, Array<Vector3>& pos) {
		Info("Sphere init 0");
		Array<real> dxs;
		Array<Vector3> tmp_pos;
		real a = 4 * pi / num_pts;
		real d = sqrt(a);
		int M_theta = (int)round(pi / d);
		real d_theta = pi / M_theta;
		real d_phi = a / d_theta;
		for (int m = 0; m < M_theta; m++) {
			real theta = pi * ((real)m + 0.5) / M_theta;
			int M_phi = (int)round(2 * pi * sin(theta) / d_phi);
			for (int n = 0; n < M_phi; n++) {
				real phi = 2 * pi * n / M_phi;
				real x = sin(theta) * cos(phi);
				real y = sin(theta) * sin(phi);
				real z = cos(theta);
				Vector3 p(R * x, R * y, R * z);
				p += center;
				tmp_pos.push_back(p);
			}
			dxs.push_back(0.5 * (d_phi + d_theta));
		}
		real avg_dx = ArrayFunc::Mean<real>(dxs);

		int actual_num_pts = tmp_pos.size();
		pts.Resize(actual_num_pts);

#pragma omp parallel for
		for (int i = 0; i < actual_num_pts; i++) {
			pos[i] = tmp_pos[i];
		}

		return R * avg_dx;
	}
}
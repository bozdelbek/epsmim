#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

#define I(a, b) ( (a) * Ny + (b) )

// only for processed_phase use
#define HD(i, j) (I((i) * 2    , (j)    )) // horizontal down (i)
#define HU(i, j) (I((i) * 2 - 2, (j)    )) // horizontal down (i - 1)
#define VR(i, j) (I((i) * 2 - 1, (j)    )) // vertical right (j)
#define VL(i, j) (I((i) * 2 - 1, (j) - 1)) // vertical left (j - 1)

typedef struct {
	double *prev;
	double *curr;
	double *next;
	double *phase;
	double *processed_phase;
	int Nx;
	int Ny;
	int Sx;
	int Sy;
} modeling_plane;

int write_to_file(char *filename, double *arr, int size) {
	int flags = O_WRONLY | O_CREAT | O_TRUNC;
	mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH;
	int fd = open(filename, flags, mode);
	if (fd == -1) {
		perror("open");
		return -1;
	}
	if (write(fd, arr, size * sizeof(double)) == -1) {
		perror("write");
		close(fd);
		return -2;
	}
	close(fd);
	return 0;
}

void process_phase(modeling_plane *plane) {
	int Nx = plane->Nx;
	int Ny = plane->Ny;
	double *phase = plane->phase;
	double *processed_phase = plane->processed_phase;


	double hd, hu, vl, vr;
	for (int i = 1; i < Nx - 1; i++) {
		for (int j = 1; j < Ny - 1; j++) {
			hd = phase[I(i    , j    )] + phase[I(i    , j - 1)];
			hu = phase[I(i - 1, j    )] + phase[I(i - 1, j - 1)];
			vr = phase[I(i    , j    )] + phase[I(i - 1, j    )];
			vl = phase[I(i    , j - 1)] + phase[I(i - 1, j - 1)];
			processed_phase[HD(i, j)] = hd;
			processed_phase[HU(i, j)] = hu;
			processed_phase[VR(i, j)] = vr;
			processed_phase[VL(i, j)] = vl;
		}
	}
}

int init_modeling_plane(modeling_plane *plane, int Nx, int Ny, int Sx, int Sy) {
	plane->Nx = Nx;
	plane->Ny = Ny;
	plane->Sx = Sx;
	plane->Sy = Sy;
	plane->processed_phase = (double*)malloc(Ny * (2 * Nx - 1) * sizeof(double));
	plane->prev = (double*)malloc(Nx * Ny * sizeof(double));
	plane->curr = (double*)malloc(Nx * Ny * sizeof(double));
	// plane->next = (double*)malloc(Nx * Ny * sizeof(double));
	plane->next = plane->prev;
	plane->phase = (double*)malloc(Nx * Ny * sizeof(double));
	if (plane->prev == NULL || plane->curr == NULL || plane->next == NULL || plane->phase == NULL) {
		perror("malloc");
		return -1;
	}
	// maybe i can use memset?
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			plane->prev[I(i,j)] = 0.0;
			plane->curr[I(i,j)] = 0.0;
			plane->next[I(i,j)] = 0.0;
			plane->phase[I(i,j)] = 0.01;
		}
	}
	for (int i = Nx / 2; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			plane->phase[I(i,j)] = 0.02;
		}
	}

	process_phase(plane);

	return 0;
}

void _process_n(double *phase_lwr, double *phase_mdl, double *res) {
	for (int i = 0; i < 4; i++) {
		res[i] = phase_lwr[i] + phase_mdl[i];
	}
	for (int i = 0; i < 3; i++) {
		res[4 + i] = phase_lwr[i] + phase_lwr[i + 1];
		res[4 + i + 3] = phase_mdl[i] + phase_mdl[i + 1];
	}
}




double f(int n, double tou) {
	double rev_gammasq = 1 / 16.0;
	double tmp = (2 * M_PI * (n * tou - 1.5));
	double exp_arg = - ( tmp * tmp * rev_gammasq);
	return exp(exp_arg) * sin(tmp) * 0.5;
}

void print_arr(double *arr, int Nx, int Ny) {
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			printf("%f ", arr[i * Ny + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void calc_step(modeling_plane *plane, double tou) {
	double *prev = plane->prev;
	double *curr = plane->curr;
	double *next = plane->next;
	double *phase = plane->phase;
	double *proc_phase = plane->processed_phase;
	int Nx = plane->Nx;
	int Ny = plane->Ny;
	int Sx = plane->Sx;
	int Sy = plane->Sy;

	static int n = 1;

	double tousq = tou * tou;
	double hy = 4.0 / (double)(Nx - 1);
	double hx = 4.0 / (double)(Ny - 1);
	
	double phixt = tou / hx;
	double phix = phixt * phixt * 0.5;
	double phiyt = tou / hy;
	double phiy = phiyt * phiyt * 0.5;

	int i = 1;
	int j = 1;

	double elem_x = 0.0;
	double elem_y = 0.0;

	double *curr_lwr = NULL;
	double *curr_mdl = NULL;
	double *curr_upr = NULL;

	double *phase_lwr = NULL;
	double *phase_mdl = NULL;

	double *prev_mdl = NULL;
	double *next_mdl = NULL;
	for (i = 1; i < Nx - 1; i++) {

		curr_lwr = curr + (i - 1) * Ny;
		curr_mdl = curr + (i) * Ny;
		curr_upr = curr + (i + 1) * Ny;

		phase_lwr = phase + (i - 1) * Ny;
		phase_mdl = phase + (i) * Ny;

		next_mdl = next + (i) * Ny;
		prev_mdl = prev + (i) * Ny;

		for (j = 1; j < Ny - 1; j++) {


			/*double hu, hd, vl, vr;
			hu = phase_lwr[j - 1] + phase_lwr[j];
			hd = phase_mdl[j - 1] + phase_mdl[j];
			vl = phase_lwr[j - 1] + phase_mdl[j - 1];
			vr = phase_lwr[j] + phase_mdl[j];*/

			double hun, hdn, vln, vrn;
			hdn = proc_phase[HD(i, j)];
			hun = proc_phase[HU(i, j)];
			vrn = proc_phase[VR(i, j)];
			vln = proc_phase[VL(i, j)];
			/*if (hd != hdn) {
				printf("i,j (%d, %d), hd: %f, hdn: %f\n", i, j, hd, hdn);
				print_arr(phase, Nx, Ny);
				print_arr(proc_phase, 2 * Nx - 1, Ny);
			}
			assert(hd == hdn);
			assert(hu == hun);
			assert(vl == vln);
			assert(vr == vrn);*/

			/*elem_x = (curr_mdl[j+1] - curr_mdl[j]) * (phase_lwr[j] + phase_mdl[j]) +
						(curr_mdl[j - 1] - curr_mdl[j]) * (phase_lwr[j - 1] + phase_mdl[j - 1]);
			elem_y = (curr_upr[j] - curr_mdl[j]) * (phase_mdl[j - 1] + phase_mdl[j]) +
						(curr_lwr[j] - curr_mdl[j]) * (phase_lwr[j - 1] + phase_lwr[j]);*/

			elem_x = (curr_mdl[j+1] - curr_mdl[j]) * (vrn) +
						(curr_mdl[j - 1] - curr_mdl[j]) * (vln);
			elem_y = (curr_upr[j] - curr_mdl[j]) * (hdn) +
						(curr_lwr[j] - curr_mdl[j]) * (hun);

			
			elem_x *= phix;
			elem_y *= phiy;

			next_mdl[j] = 2.0 * curr_mdl[j] - prev_mdl[j] + elem_x + elem_y;

		}
	}
	next[I(Sx, Sy)] += tousq * f(n, tou);

	plane->prev = curr;
	plane->curr = next;
	plane->next = prev;
	plane->next = plane->prev;
	n++;
}

void print_m(double *arr, int Nx, int Ny) {
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (arr[I(i,j)] > 1000 || arr[I(i,j)] < -1000) {
				printf("(%d, %d) %f\n", i, j, arr[I(i,j)]);
			}
			
		}
	}
}

void print_csv(double *arr, int Nx, int Ny) {
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			printf("%f;", arr[I(i,j)]);
		}
		printf("\b\n");
	}
}

int main(int argc, char *argv[]) {
	// double tou = 0.01;
	int Nx = 0, Ny = 0, Nt = 0;
	int opt = 0;
	while ( (opt = getopt(argc, argv, "x:y:t:")) != -1 ) {
		switch (opt) {
			case 'x':
				Nx = atoi(optarg);
				break;
			case 'y':
				Ny = atoi(optarg);
				break;
			case 't':
				Nt = atoi(optarg);
				break;
			case '?':
				printf("error: no such arg\n");
				break;
		}
	}
	printf("Nx: %d, Ny: %d, Nt: %d\n", Nx, Ny, Nt);

	int Sx = Nx / 2;
	int Sy = Ny / 2;
	
	modeling_plane plane;
	if (init_modeling_plane(&plane, Nx, Ny, Sx, Sy) == -1) {
		fprintf(stderr, "error in init_modeling_plane\n");
		exit(-1);
	}

	double tou = 0.01;
	int iters = (double)Nt / tou;
	printf("Nt: %d, tou: %f, iters: %d\n", Nt, tou, iters);

	char fname[100] = { 0 };

	for (int i = 1; i < iters; i++) {
		calc_step(&plane, tou);
		// sprintf(fname, "prev%d", i);
		// write_to_file(fname, plane.prev, Nx * Ny);
		sprintf(fname, "data/ncurr%d", i);
		if (i % 50 == 0) {
			write_to_file(fname, plane.curr, Nx * Ny);
		}
		/*if (i == Nt) {
			print_csv(plane.curr, Nx, Ny);
		}*/
		// sprintf(fname, "next%d", i);
		// write_to_file(fname, plane.next, Nx * Ny);
		/*if (i % 10 == 0) {
			write_to_file(fname, plane.curr, Nx * Ny);
		}*/
	}

	return 0;
}
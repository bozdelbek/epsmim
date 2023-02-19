#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

#define I(a, b) ( (a) * Nx + (b) )

typedef struct {
	double *prev;
	double *curr;
	double *next;
	double *phase;
	int Nx;
	int Ny;
	int Sx;
	int Sy;
} modeling_plane;

int write_to_file(char *filename, double *arr, int size) {
	int flags = O_WRONLY | O_CREAT | O_TRUNC | (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
	int fd = open(filename, flags);
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

int init_modeling_plane(modeling_plane *plane, int Nx, int Ny, int Sx, int Sy) {
	plane->Nx = Nx;
	plane->Ny = Ny;
	plane->Sx = Sx;
	plane->Sy = Sy;
	plane->prev = (double*)malloc(Nx * Ny * sizeof(double));
	plane->curr = (double*)malloc(Nx * Ny * sizeof(double));
	plane->next = (double*)malloc(Nx * Ny * sizeof(double));
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
	return 0;
}

double f(int n, double tou) {
	double rev_gammasq = 1 / 16.0;
	double tmp = (2 * M_PI * (n * tou - 1.5));
	double exp_arg = - ( tmp * tmp * rev_gammasq);
	return exp(exp_arg) * sin(tmp) * 0.5;
}

void calc_step(modeling_plane *plane, double tou) {
	double *prev = plane->prev;
	double *curr = plane->curr;
	double *next = plane->next;
	double *phase = plane->phase;
	int Nx = plane->Nx;
	int Ny = plane->Ny;
	int Sx = plane->Sx;
	int Sy = plane->Sy;
	static int n = 1;
	double tousq = tou * tou;
	double rev_hxsq = (double)((Nx - 1)*(Nx - 1)) / 16.0;
	double rev_hysq = (double)((Ny - 1)*(Ny - 1)) / 16.0;
	double hy = 4.0 / (double)(Nx - 1);
	double hx = 4.0 / (double)(Ny - 1);
	double hxsq = hx * hx;
	double hysq = hy * hy;
	for (int i = 1; i < Nx - 1; i++) {
		for (int j = 1; j < Ny - 1; j++) {

			//don't forget to mul by hx\hy
			// double elem_x = (curr[I(i, j+1)] - curr[I(i, j)]) * (phase[I(i-1, j)] + phase[I(i, j)]) +
			// 				(curr[I(i, j-1)] - curr[I(i, j)]) * (phase[I(i-1, j-1)] + phase[I(i, j-1)]);
			// double elem_y = (curr[I(i+1, j)] - curr[I(i, j)]) * (phase[I(i, j-1)] + phase[I(i, j)]) +
			// 				(curr[I(i-1, j)] - curr[I(i, j)]) * (phase[I(i-1, j-1)] + phase[I(i-1, j)]);
			double elem_x = (curr[I(i, j+1)] - curr[I(i, j)]) + (curr[I(i, j-1)] - curr[I(i, j)]);
			double elem_y = (curr[I(i+1, j)] - curr[I(i, j)]) +	(curr[I(i-1, j)] - curr[I(i, j)]);
			elem_x *= 0.02;
			elem_y *= 0.02;
			// next[I(i, j)] = 2.0 * curr[I(i, j)] - prev[I(i, j)] + tousq * (elem_x * rev_hxsq + elem_y * rev_hysq) * 0.5;
			next[I(i, j)] = 2.0 * curr[I(i, j)] - prev[I(i, j)] + tousq * (elem_x / (2 * hxsq) + elem_y / (2 * hysq));
			if (next[I(i,j)] < -10 || 10 < next[I(i,j)]) {
				printf("(%d, %d, [%d]) %lf %lf %lf\n", i, j, n, elem_x, elem_y, next[I(i,j)]);
				usleep(100000);
			}
			// assert(next[I(i,j)] < 1000.0);
			/*if (-1000.0 > next[I(i,j)] || next[I(i,j)] > 1000.0) {
				printf("i: %d, j: %d (U = %f)\n", i, j, next[I(i,j)]);
				sleep(1);
				// assert(next[I(i,j)] < 1000);
			}*/
		}
	}
	next[I(Sx, Sy)] += tousq * f(n, tou);
	plane->prev = curr;
	plane->curr = next;
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

	double tou = 0.1;

	char fname[100] = { 0 };

	for (int i = 1; i < 300; i++) {
		calc_step(&plane, tou);
		// sprintf(fname, "prev%d", i);
		// write_to_file(fname, plane.prev, Nx * Ny);
		sprintf(fname, "curr%d", i);
		write_to_file(fname, plane.curr, Nx * Ny);
		if (i == 1) {
			print_m(plane.curr, Nx, Ny);
		}
		// sprintf(fname, "next%d", i);
		// write_to_file(fname, plane.next, Nx * Ny);
		/*if (i % 10 == 0) {
			write_to_file(fname, plane.curr, Nx * Ny);
		}*/
	}

	return 0;
}
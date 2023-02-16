#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#define I(a, b) ( (a) * Nx + (b) )

double f(int n, double tou) {
	double rev_gammasq = 1 / 16.0;
	double tmp = (2 * M_PI * (n * tou - 1.5));
	double exp_arg = - ( tmp * tmp * rev_gammasq);
	return exp(exp_arg) * sin(tmp);
}

void calc_step(double *prev, double *curr, double *next, double *phase, double tou, int Nx, int Ny, int Sx, int Sy) {
	static int n = 1;
	double tousq = tou * tou;
	double rev_hxsq = (double)(Nx - 1)*(Nx - 1) / 16.0;
	double rev_hysq = (double)(Ny - 1)*(Ny - 1) / 16.0;
	for (int i = 1; i < Nx - 1; i++) {
		for (int j = 1; j < Ny - 1; j++) {
			//don't forget to mul by hx\hy
			double elem_x = (curr[I(i, j+1)] - curr[I(i, j)]) * (phase[I(i-1, j)] + phase[I(i, j)]) +
							(curr[I(i, j-1)] - curr[I(i, j)]) * (phase[I(i-1, j-1)] + phase[I(i, j-1)]);
			double elem_y = (curr[I(i+1, j)] - curr[I(i, j)]) * (phase[I(i, j-1)] + phase[I(i, j)]) +
							(curr[I(i-1, j)] - curr[I(i, j)]) * (phase[I(i-1, j-1)] + phase[I(i-1, j)]);
			next[I(i, j)] = 2.0 * curr[I(i, j)] - prev[I(i, j)] + tousq * (elem_x * rev_hxsq + elem_y * rev_hysq);
		}
	}
	next[Sy * Nx + Sx] += tousq * f(n, tou);
	n++;
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
	
	double *prev = (double*)malloc(Nx * Ny * sizeof(double));
	double *curr = (double*)malloc(Nx * Ny * sizeof(double));
	double *next = (double*)malloc(Nx * Ny * sizeof(double));
	double *phase = (double*)malloc(Nx * Ny * sizeof(double));
	// ---------------
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			prev[i * Ny + j] = 0.0;
			curr[i * Ny + j] = 0.0;
			phase[i * Ny + j] = 0.0;
		}
	}

	return 0;
}
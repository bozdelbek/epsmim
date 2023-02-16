#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

double f(int n, double tou) {
	double rev_gammasq = 1 / 16.0;
	double tmp = (2 * M_PI * (n * tou - 1.5));
	double exp_arg = - ( tmp * tmp * rev_gammasq);
	return exp(exp_arg) * sin(tmp);
}

double mysin(double x) {
	return sin(x);
}

double myexp(double x) {
	return exp(-x);
}

int main2() {
	for (double x = 0.0; x < 500.0; x += 0.1) {
		printf("%f %f\n", x, f(x, 0.01));
		// printf("%f %f\n", x, myexp(x) * mysin(10.0*x));
	}
	return 0;
}

int main() {
	printf("%d\n", sizeof(double));
	double arr[100] = { 0 };
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			arr[i * 10 + j] = (double)i * j;
		}
	}
	int fd = open("file.dat", O_RDWR | O_CREAT);
	write(fd, arr, 100 * sizeof(double));
	close(fd);
	return 0;
}
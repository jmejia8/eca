#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double randm(){
	return (double)rand() / (double)RAND_MAX ;
}

void printArray(double *array, int N, int D){
	int i, l = N*D;
	for (i = 0; i < l; ++i) {
		if (i> 0 && i % D ==0)
			printf("\n");

		printf("%.3e, ", array[i]);
	}

	printf("\n");
}

void arrayRand(double *array, double a, double b, int N){
	int i;

	double l = b - a;
	for (i = 0; i < N; ++i)
		array[i] = a + l*randm();

}

void evalPopulation(double (*f)(double*, int), double* array, double* fitness, int N, int D){
	int i;
	for (i = 0; i < N; ++i) {
		fitness[i] = (*f)(&array[i*D], D);
	}
}

void randperm(int *r, int n){
	int i;
	for (i = n-1; i >= 0; --i){
		//generate a random number [0, n-1]
		int j = rand() % (i+1);

		//swap the last element with element at random index
		int temp = r[i];
		r[i] = r[j];
		r[j] = temp;
	}
}

void genU(int *U, int N){
	int i;
	for (i = 0; i < N; ++i)
		U[i] = i;

	randperm(U, N);
}

double sum(double *x, int n){
	int i;
	double s = 0;
	for (i = 0; i < n; ++i)
		s += x[i];

	return s;
}

void center(double *c, double *population, double *m, int *U, int* best, int* worst, int K, int D){
	int i, j;
	double M = 0;
	double mini = m[U[0]];

	best[0]  = 0;
	worst[0] = 0;
	
	double m_best = mini + m[best[0]];
	double m_worst= mini + m[worst[0]];

	for (i = 0; i < D; ++i){
		c[i] = 0;

		if (mini > m[U[i]])	{
			mini = m[U[i]];
		}
	}

	if (mini >= 0) mini = 0;
	else mini = 2.0*abs(mini);
	
	for (i = 0; i < K; ++i) {
		double *x = &population[ U[i] ];
		double mass = mini + m[U[i]];

		for (j = 0; j < D; ++j)
			c[j] += x[j] * mass;

		M += mass;

		if (mass > m_best){
			m_best = mass;
			best[0] = U[i];
		}

		if (mass < m_worst){
			m_worst = mass;
			worst[0] = U[i];
		}

	}

	for (i = 0; i < D; ++i) c[i] /= M;

}


void genCenters(double *c, double *population, double *m, int *U, int *best, int *worst, int K, int N, int D){
	int i,ii = 0, k;
	for (i = 0; i < N; ++i) {
		k = ii*K;
		
		if (k >= N-K){
			genU(U, N);
			ii = 0;
			k = 0;
		}

		center(&c[i*D], population, m, &U[k], &best[i], &worst[i], K, D);

		ii++;
	}
}

void genEtas(double* etas, double eta_min, double eta_max, int N){
	arrayRand(etas, eta_min, eta_max, N);
}

void variationOperator(double* h, double *x, double *c, double eta, double* u, int D){
	int i;

	for (i = 0; i < D; ++i)
		h[i] = x[i] + eta*(c[i] - u[i]);

}

void genh(double* h, double *x, double *c, double *etas, int* worst, int N, int D){
	int i;

	for (i = 0; i < N; ++i) {
		int ii = i*D;
		variationOperator(&h[ii], &x[ii], &c[ii], etas[i], &x[worst[i]], D);
	}
}

void replace(double *x, double *h, int D){
	int i;
	for (i = 0; i < D; ++i)
		x[i] = h[i];
}


int maxind(double *f, int N){
	int i;
	int max = 0;
	for (i = 1; i < N; ++i) {
		if (f[i] > f[max]) {
			max = i;
		}
	}

	return max;
}

int minind(double *f, int N){
	int i;
	int min = 0;
	for (i = 1; i < N; ++i) {
		if (f[i] < f[min]) {
			min = i;
		}
	}

	return min;
}

void mutation(double *h, double *x, double *c, double* u_best, double* u_worst, double eta, double a, double b, double P_bin, int D){
	int i;
	double l = b - a;
	for (i = 0; i < D; ++i) {
		if (randm() < P_bin){
			h[i] = u_best[i];
			continue;
		}

		h[i] = x[i] + eta*(c[i] - u_worst[i]);
		
		if ( h[i] < a || b < h[i])
			h[i] = a + l*randm();

	}
}

void eca(double (*f)(double*, int), int D, int N){
	int K = 7;
	double eta_max= 2.0;
	double P_bin  = 0.03;
	int max_evals = 10000*D;
	double a = -10, b = 10;

	int stop = 0;

	double* population = (double *) malloc(sizeof(double)*D*N);
	double* fitness    = (double *) malloc(sizeof(double)*N);


	int *U = (int *) malloc(sizeof(int) * N);
	int best, worst;
	
	double *c  = (double *) malloc(sizeof(double) * D);
	double* h  = (double *) malloc(sizeof(double) * D);
	
	// init population
	arrayRand(population, a, b, N * D);
	
	// evaluate population
	evalPopulation(f,population, fitness, N, D);

	int nevals = N;
	int t,i;
	do
	{
		int ii = 0;
		genU(U, N);
	
		for (i = 0; !stop && i < N; ++i) {
			int k = ii*K;
		
			if (k >= N-K){
				genU(U, N);
				ii = 0;
				k = 0;
			}

			// current solution
			double *x = &population[i];

			// generate a center of mass
			center(c, population, fitness, &U[k], &best, &worst, K, D);

			// stepsize
			double eta = eta_max*randm();
			double *u_worst = &population[worst];
			double *u_best = &population[worst];

			// variation operator
			mutation(h, x, c, u_best, u_worst, eta, a, b, P_bin, D);

			double fh = (*f)(h, D);
			nevals += 1;

			// compare solutions
			if (fitness[i] < fh) {
				int w = minind(fitness, N);
				replace(&population[w], h, D);
				fitness[w] = fh;
			}

			stop = nevals >= max_evals;
			++ii;
		}

		++t;
	} while (!stop);

	int best_i = maxind(fitness, N);

	printArray(&population[best_i], 1, D);
	printf("%e\n", fitness[best_i]);
	printf("%d\n", nevals);
}

double sphere(double* array, int D){
	int i; double s = 0;
	
	for (i = 0; i < D; ++i)
		s += pow(array[i]-1, 2);

	return -s;
}

double gauss(double* array, int D){
	return exp(-sphere(array, D));
}


int main(int argc, char const *argv[])
{
	srand(time(NULL));
	int D = 10;
	int N = 7*D;
	eca(sphere, D, N);

	return 0;
}
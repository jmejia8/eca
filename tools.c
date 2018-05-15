#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef struct _result
{
    double f;
    double *x;
} Result;

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

void evalPopulation(double (*f)(double*, int), double* array, double* fitness,
                    int N, int D){
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

void center(double *c, double *population, double *m, int *U, int* best,
            int* worst, int K, int D){
    int i, j;
    double M = 0;
    double mini = m[U[0]];

    best[0]  = 0;
    worst[0] = 0;
    
    double m_best = mini + m[best[0]];
    double m_worst= mini + m[worst[0]];

    for (i = 0; i < D; ++i){
        c[i] = 0;

        if (mini > m[U[i]]) {
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


void genCenters(double *c, double *population, double *m, int *U, int *best,
                int *worst, int K, int N, int D){
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

void mutation(double *h, double *x, double *c, double* u_best, double* u_worst,
              double eta, double a, double b, double P_bin, int D){
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
#include "tools.c"

/*************************************************************************/
/*  eca.c                                                                */
/*************************************************************************/
/* Evolutionary Centers Algorithm                                        */
/*************************************************************************/
/*                                                                       */ 
/* Parameters (suggested):                                               */
/*      Dimension: D                                                     */
/*      K-value:                                                         */
/*            K = 7                                                      */
/*      Population size:                                                 */
/*            N = K*D                                                    */
/*      stepsize:                                                        */
/*            eta_max = 2.0                                              */
/*      binomial probability:                                            */
/*            P_bin = 0.03                                               */
/*      Max. number of evaluations:                                      */
/*            max_evals = 10000*D                                        */
/*                                                                       */
/* Bounds:                                                               */
/*      Lower: low_bound                                                 */
/*      Upper: up_bound                                                  */
/*                                                                       */
/* Search Type:                                                          */
/*     Maximize:                                                         */
/*         searchType = 1                                                */
/*     minimize:                                                         */
/*         searchType = 0                                                */
/*                                                                       */
/*************************************************************************/

Result eca(double (*f)(double*, int),
            int D,
            int N,
            int K,
            double eta_max,
            double P_bin,
            int    max_evals,
            double low_bound,
            double up_bound,
            int    searchType
            )
{
    double a = low_bound, b = up_bound;

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
    evalPopulation(f, population, fitness, N, D, searchType);

    // current number of evaluations
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

            double fh;
            if (searchType) // maximize
                fh = (*f)(h, D);
            else // maximize
                fh = -(*f)(h, D);

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

    Result result;
    result.x = (double *) malloc(sizeof(double) * D);

    if (searchType) // maximize
        result.f = fitness[best_i];
    else // minimize
        result.f = -fitness[best_i];

    best_i *= D;
    for (i = 0; i < D; ++i) {
        result.x[i] = population[best_i + i];
    }

    printf("===========[ ECA results ]=============\n");
    printf("| nevals = %d\n", nevals);
    printf("| f      = %e\n", result.f);
    printf("| x      = "); printArray(result.x, 1, D);
    printf("=======================================\n");

    return result;
}
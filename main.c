#include "eca.c"
#include "test_functions.c"


int main(int argc, char const *argv[])
{
    srand(time(NULL));

    // ECA parameters
    int D = 10;
    int K = 7;
    int N = K*D;
    double eta_max = 2.0;
    double P_bin = 0.03;
    int    max_evals = 10000*D;
    double low_bound = -100;
    double up_bound = 100;
    int    searchType = 0; // minimize

    // optimize
    Result result = eca(sphere, D, N, K,
                        eta_max,
                        P_bin,
                        max_evals,
                        low_bound,
                        up_bound,
                        searchType);

    return 0;
}
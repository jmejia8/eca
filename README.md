# Evolutionary Centers Algorithm

ECA is a physics-inspired algorithm based on the center of mass concept on 
a D-dimensional space for real-parameter single-objective optimization. The 
general idea is to promote the creation of an irregular body using K mass points
in the current population, then the center of mass is calculated to get a new direction 
for the next population... [read more.](https://www.dropbox.com/s/kqc22ki2edjtt0y/ECA-optimization.pdf)

## Parameters
- Parameters (suggested):
    - Dimension: `D`
    - K-value:
           `K = 7`
    - Population size:
           `N = K*D`
    - stepsize:
           `eta_max = 2.0`
    - binomial probability:
           `P_bin = 0.03`
     - Max. number of evaluations:
           `max_evals = 10000*D`

- Bounds:
     - Lower: `low_bound`
     - Upper: `up_bound`

- Search Type:
    - Maximize:
        - `searchType = 1`
    - minimize:
        - `searchType = 0`


## Example

You can write C code to use ECA in your project:

```C
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
    double low_bound = -10;
    double up_bound = 10;
    int    searchType = 1; // maximize

    // optimize
    Result result = eca(gauss, D, N, K,
                        eta_max,
                        P_bin,
                        max_evals,
                        low_bound,
                        up_bound,
                        searchType);

    return 0;
}
```

Also, you can build the ECA library by running `Make`
#include "eca.c"


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

    Result result = eca(sphere, D, N);

    return 0;
}
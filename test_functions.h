// For minimizing
double sphere(double* array, int D){
    int i; double s = 0;
    
    for (i = 0; i < D; ++i)
        s += pow(array[i]-1, 2);

    return s;
}

// For maximizing
double gauss(double* array, int D){
    return exp(-sphere(array, D));
}

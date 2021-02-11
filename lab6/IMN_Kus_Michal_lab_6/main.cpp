#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "mgmres.h"

using namespace std;

double delta = 0.1;
double eps1 = 1.0;
double eps2 = 1.0;

double localEpsilon(const int nx, const int l) {
    int j = floor(l / (nx+1.0));
    int i = l - j*(nx+1.0);
    return (i <= nx/2.)? eps1 : eps2;
}

int fillSparseMatrix(string matrixFileName, string vectorFileName,
                     const int N, const vector<double> boundCond,
                     vector<double> &a, vector<int> &ia, vector<int> &ja,
                     vector<double> &b, const int nx, const int ny,
                     vector<double> &rho) {
    fstream fs;
    fs.open(vectorFileName, fstream::out | fstream::trunc);
    int k = -1;
    int i;
    int j;

    for(int l = 0; l < N; l++)
    {
        int boundary = 0;
        double Vb = 0.0;
        
        j = floor(l / (nx+1.0));
        i = l - j*(nx+1.0);

        if(i == 0) 
        {
            boundary = 1;
            Vb = boundCond[0];
        }

        if(j == ny)
        {
            boundary = 1;
            Vb = boundCond[1];
        }

        if(i == nx) 
        {
            boundary = 1;
            Vb = boundCond[2];
        }

        if(j == 0)
        {
            boundary = 1;
            Vb = boundCond[3];
        }   

        b[l] = -rho[l];
        if(boundary == 1)
            b[l] = Vb;

        ia[l] = -1;

        if(l - nx - 1 >= 0 && boundary == 0)
        {
            k++;
            if(ia[l] < 0)
                ia[l] = k;

            a[k] = localEpsilon(nx, l) / (delta * delta);
            ja[k] = l - nx - 1;
        }

        if(l - 1 >= 0 && boundary == 0)
        {
            k++;
            if(ia[l] < 0)
                ia[l] = k;

            a[k] = localEpsilon(nx, l) / (delta * delta);
            ja[k] = l - 1;
        }


        k++;
        if(ia[l] < 0)
            ia[l] = k;

        if(boundary == 0)
            a[k] = -(2*localEpsilon(nx, l) + localEpsilon(nx, l+1) + localEpsilon(nx, l+nx+1))/(delta * delta);
        else
            a[k] = 1;
        
        ja[k] = l;

        if(l < N && boundary == 0)
        {
            k++;
            a[k] = localEpsilon(nx, l+1)/(delta * delta);
            ja[k] = l+1;
        }

        if(l < N - nx - 1 && boundary == 0)
        {
            k++;
            a[k] = localEpsilon(nx, l+nx+1)/(delta * delta);
            ja[k] = l + nx + 1;
        }
        fs << l << ", " << i << ", " << j << ", " << b[l] << "\n";
    }
    fs.close();
    
    int nonZeros = k+1;
    ia[N] = nonZeros;
    fs.open(matrixFileName, fstream::out | fstream::trunc);
    for(int k = 0; k < 5*N; k++)
    {
        fs << k << ", " << a[k] << "\n";
    }
    fs.close();

    return nonZeros;
}

void solveSet(string matrixFileName, string vectorFileName, string outputFileName,
              const int N, const vector<double> boundCond,
              vector<double> &aV, vector<int> &iaV, vector<int> &jaV,
              vector<double> &bV, vector<double> &VV, const int nx, const int ny, vector<double> &rho) {
    fstream fs;
    int itr_max = 500;
    int mr = 500;
    double tol_abs = 1e-8;
    double tol_rel = 1e-8;
    vector<double> bVector;
    vector<double> VVector;

    bVector.resize(N, 0.0);
    VVector.resize(N, 0.0);

    fillSparseMatrix(matrixFileName, vectorFileName, N, boundCond, aV, iaV, jaV, bVector, nx, ny, rho);

    int* ia = &iaV[0];
    int* ja = &jaV[0];
    double* a = &aV[0];

    double* b = &bVector[0];
    double* V = &VVector[0];

    pmgmres_ilu_cr(N, ia[N], ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

    fs.open(outputFileName, fstream::out);
    for(int l = 0; l < N; l++)
    {
        int j = floor(l / (nx + 1.));
        int i = l - j*(nx+1);
        fs << delta*j << ", " << delta*i << ", " << V[l] << "\n"; 
    }
    fs.close();
}

int main() {
    int nx = 4;
    int ny = 4;
    double V1 =  10.0;
    double V3 = 10.0;
    double V2 = -10.0;
    double V4 = -10.0;
    int N = (nx + 1)*(ny + 1);

    vector<double> a, b, V;
    vector<double> Vb{V1, V2, V3, V4};
    vector<int> ia, ja;
    vector<double> rhoMatrix;

    string matrixFileName = "matrix.txt";
    string vectorFileName = "vector.txt";

    rhoMatrix.resize(N, 0.0);
    a.resize(5*N, 0.0);
    ia.resize(N+1, -1);
    ja.resize(5*N, 0);
    b.resize(N, 0.0);
    V.resize(N, 0.0);

    solveSet(matrixFileName, vectorFileName, "map_4.txt", N, Vb, a, ia, ja, b, V, nx, ny, rhoMatrix);

    nx = 50;
    ny = 50;
    N = (nx + 1)*(ny + 1);
    a.resize(5*N, 0.0);
    ia.resize(N+1, -1);
    ja.resize(5*N, 0);
    b.resize(N, 0.0);
    V.resize(N, 0.0);
    rhoMatrix.resize(N, 0.0);

    solveSet(matrixFileName, vectorFileName, "map_50.txt", N, Vb, a, ia, ja, b, V, nx, ny, rhoMatrix);

    nx = 100;
    ny = 100;
    N = (nx + 1)*(ny + 1);
    a.resize(5*N, 0.0);
    ia.resize(N+1, -1);
    ja.resize(5*N, 0);
    b.resize(N, 0.0);
    V.resize(N, 0.0);
    rhoMatrix.resize(N, 0.0);

    solveSet(matrixFileName, vectorFileName, "map_100.txt", N, Vb, a, ia, ja, b, V, nx, ny, rhoMatrix);

    nx = 200;
    ny = 200;
    N = (nx + 1)*(ny + 1);
    a.resize(5*N, 0.0);
    ia.resize(N+1, -1);
    ja.resize(5*N, 0);
    b.resize(N, 0.0);
    V.resize(N, 0.0);
    rhoMatrix.resize(N, 0.0);

    solveSet(matrixFileName, vectorFileName, "map_200.txt", N, Vb, a, ia, ja, b, V, nx, ny, rhoMatrix);

    nx = 100;
    ny = 100;
    N = (nx + 1)*(ny + 1);
    a.resize(5*N, 0.0);
    ia.resize(N+1, -1);
    ja.resize(5*N, 0);
    b.resize(N, 0.0);
    V.resize(N, 0.0);
    rhoMatrix.resize(N, 0.0);
    Vb = {0, 0, 0, 0};
    double xmax = delta * nx;
    double ymax = delta * ny;
    double sigma = xmax/10.;
    double rho1, rho2;

    // a
    eps1 = 1.;
    eps2 = 1.;

    for(int l = 0; l < N; l++)
    {
        int j = floor(l / (nx + 1.));
        int i = l - j*(nx+1);
        rho1 = exp(-pow(j - 0.25*ymax, 2)/(sigma * sigma) - pow(i - 0.5*xmax, 2)/(sigma * sigma));
        rho2 = -exp(-pow(j - 0.75*ymax, 2)/(sigma * sigma) - pow(i - 0.5*xmax, 2)/(sigma * sigma));
        rhoMatrix[l] = rho1 + rho2;
    }

    solveSet(matrixFileName, vectorFileName, "map_100_ro1.txt", N, Vb, a, ia, ja, b, V, nx, ny, rhoMatrix);
    // b
    eps1 = 1.;
    eps2 = 2.;
    N = (nx + 1)*(ny + 1);
    a.resize(5*N, 0.0);
    ia.resize(N+1, -1);
    ja.resize(5*N, 0);
    b.resize(N, 0.0);
    V.resize(N, 0.0);

    solveSet(matrixFileName, vectorFileName, "map_100_ro2.txt", N, Vb, a, ia, ja, b, V, nx, ny, rhoMatrix);
    // c
    eps1 = 1.;
    eps2 = 10.;
    N = (nx + 1)*(ny + 1);
    a.resize(5*N, 0.0);
    ia.resize(N+1, -1);
    ja.resize(5*N, 0);
    b.resize(N, 0.0);
    V.resize(N, 0.0);

    solveSet(matrixFileName, vectorFileName, "map_100_ro10.txt", N, Vb, a, ia, ja, b, V, nx, ny, rhoMatrix);

    return 0;
}
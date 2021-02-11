#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

#define OUT

using namespace std;

struct Params {

    Params(double d, double r, double u, int x, int y, int i, int j, int max) {
        delta = d;
        rho = r;
        mi = u;
        nx = x;
        ny = y;
        i1 = i;
        j1 = j;
        IT_MAX = max;
    }

    double delta;
    double rho;
    double mi;
    int nx;
    int ny;
    int i1;
    int j1;
    int IT_MAX;
};

void printMatrix(const string& filename, const Params& params, vector<vector<double>>& matrix) {
    fstream fs;
    fs.open(filename, fstream::out);
    for(int i = 0; i < params.nx + 1; i++){
        for(int j = 0; j < params.ny + 1; j++){
            fs << params.delta*i << ", " << params.delta*j << ", " << matrix[i][j] << endl;
        }
    }
    fs.close();
}

void setBoundaryConditionsPsi(const Params& params, vector<vector<double>>& psi, double Qwe) {
    double x, y;
    double y_ny = params.delta * params.ny; 
    double y_j1 = params.delta * params.j1;
    double Qwy = Qwe * (y_ny*y_ny*y_ny - y_j1*y_j1*y_j1 - 3*y_j1*y_ny*y_ny + 3*y_j1*y_j1*y_ny)/(y_ny*y_ny*y_ny);

// A
    for(int j = params.j1; j < params.ny + 1; j++)
    {
        y = j * params.delta;
        psi[0][j] = (Qwe/(2.*params.mi))*(((y*y*y)/3.) - (y_j1 + y_ny) * (y*y)/2. 
                    + y*y_j1*y_ny);
    } 

// C
    for(int j = 0; j < params.ny + 1; j++)
    {
        y = j * params.delta;
        psi[params.nx][j] = (Qwy/(2.*params.mi))*(((y*y*y)/3.) - y_ny * (y*y)/2.) 
                    + (Qwe*y_j1*y_j1*(-y_j1 + 3*y_ny))/(12.*params.mi);
    } 

// B
    for(int i = 1; i < params.nx; i++)
    {
        psi[i][params.ny] = psi[0][params.ny];
    }

// D
    for(int i = params.i1 + 1; i < params.nx; i++)
    {
        psi[i][0] = psi[0][params.j1];
    }

// E
    for(int j = 1; j < params.j1 + 1; j++)
    {
        psi[params.i1][j] = psi[0][params.j1];
    }

// F
    for(int i = 1; i < params.i1 + 1; i++)
    {
        psi[i][params.j1] = psi[0][params.j1];
    }
}

void setBoundaryConditionsZeta(const Params& params, vector<vector<double>>& psi, vector<vector<double>>& zeta, double Qwe) {
    double x, y;
    double y_ny = params.delta * params.ny; 
    double y_j1 = params.delta * params.j1;
    double Qwy = Qwe * (y_ny*y_ny*y_ny - y_j1*y_j1*y_j1 - 3*y_j1*y_ny*y_ny + 3*y_j1*y_j1*y_ny)/(y_ny*y_ny*y_ny);

// A
    for(int j = params.j1; j < params.ny + 1; j++)
    {
        y = j * params.delta;
        zeta[0][j] = (Qwe/(2.*params.mi))*(2.*y - y_j1 - y_ny);
    } 

// C
    for(int j = 0; j < params.ny + 1; j++)
    {
        y = j * params.delta;
        zeta[params.nx][j] = (Qwy/(2.*params.mi))*(2.*y - y_ny);
    } 

// B
    for(int i = 1; i < params.nx; i++)
    {
        zeta[i][params.ny] = (2./(params.delta*params.delta))*(psi[i][params.ny-1] - psi[i][params.ny]);
    }

// D
    for(int i = params.i1 + 1; i < params.nx; i++)
    {
        zeta[i][0] = (2./(params.delta*params.delta))*(psi[i][1] - psi[i][0]);
    }

// E
    for(int j = 1; j < params.j1; j++)
    {
        zeta[params.i1][j] = (2./(params.delta*params.delta))*(psi[params.i1+1][j] - psi[params.i1][j]);
    }

// F
    for(int i = 1; i < params.i1+1; i++)
    {
        zeta[i][params.j1] = (2./(params.delta*params.delta))*(psi[i][params.j1+1] - psi[i][params.j1]);
    }

    zeta[params.i1][params.j1] = 0.5*(zeta[params.i1-1][params.j1] + zeta[params.i1][params.j1-1]);
}

void relaxation(const Params& params, vector<vector<double>>& psi, vector<vector<double>>& zeta, double Qwe) {
    double omega, x, y, gamma;
    int j2 = params.j1 + 2;
    setBoundaryConditionsPsi(params, psi, Qwe);

    fstream fs;
    string filename = to_string(Qwe) + "err.txt";
    fs.open(filename, fstream::out);

    for(int it = 1; it <= params.IT_MAX; it++)
    {
        if(it < 2000)
            omega = 0.;
        else
            omega = 1.;

        for(int i = 1; i < params.nx; i++)
        {
            x = params.delta * i;
            for(int j = 1; j < params.ny; j++)
            {
                y = params.delta * j;
                if((i <= params.i1 && j > params.j1) || (i > params.i1))
                {
                    psi[i][j] = 0.25*(psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - params.delta*params.delta*zeta[i][j]);
                    zeta[i][j] = 0.25*(zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1]);
                    zeta[i][j] -= omega*((params.rho)/(16.*params.mi)) *
                                 ((psi[i][j+1] - psi[i][j-1])*(zeta[i+1][j] - zeta[i-1][j]) - (psi[i+1][j] - psi[i-1][j])*(zeta[i][j+1] - zeta[i][j-1]));
                }
            }
        }

        setBoundaryConditionsZeta(params, psi, zeta, Qwe);
        gamma = 0;
        for(int i = 1; i < params.nx; i++)
        {
            gamma += psi[i+1][j2] + psi[i-1][j2] + psi[i][j2+1] + psi[i][j2-1] - 4.*psi[i][j2] - params.delta*params.delta*zeta[i][j2];
        }

        fs << it << ", " << gamma << endl;
    }
    fs.close();
}

void clearTables(const Params& params, vector<vector<double>>& psi, vector<vector<double>>& zeta, vector<vector<double>>& u, vector<vector<double>>& v) {
    for(auto& row : psi)
        fill(row.begin(), row.end(), 0.0);
    for(auto& row : zeta)
        fill(row.begin(), row.end(), 0.0);
    for(auto& row : u)
        fill(row.begin(), row.end(), 0.0);
    for(auto& row : v)
        fill(row.begin(), row.end(), 0.0);
}

void calculateSpeeds(const Params& params, vector<vector<double>>& psi, vector<vector<double>>& zeta, vector<vector<double>>& u, vector<vector<double>>& v) {
    for(int i = 1; i < params.nx; i++)
    {
        for(int j = 1; j < params.ny; j++)
        {
            if((i <= params.i1 && j > params.j1) || (i > params.i1))
            {
                u[i][j] = (psi[i][j+1] - psi[i][j-1])/(2.*params.delta);
                v[i][j] = (-psi[i+1][j] + psi[i-1][j])/(2.*params.delta);
            }
            else
            {
                u[i][j] = 0.0;
                v[i][j] = 0.0;
            }
        }
    }
}


int main() {
    Params params{0.01, 1., 1., 200, 90, 50, 55, 20000};
    vector<vector<double>> psi;
    vector<vector<double>> zeta;
    vector<vector<double>> u;
    vector<vector<double>> v;
    string filename;
    
    psi.resize(params.nx + 1);
    for(auto& x : psi) {
        x.resize(params.ny + 1, 0.0);
    }

    zeta.resize(params.nx + 1);
    for(auto& x : zeta) {
        x.resize(params.ny + 1, 0.0);
    }

    u.resize(params.nx + 1);
    for(auto& x : u) {
        x.resize(params.ny + 1, 0.0);
    }

    v.resize(params.nx + 1);
    for(auto& x : v) {
        x.resize(params.ny + 1, 0.0);
    }

    double Q = -1000.;
    filename = "-Q1k.txt";

    clearTables(params, psi, zeta, u, v);
    setBoundaryConditionsPsi(params, psi, Q);
    setBoundaryConditionsZeta(params, psi, zeta, Q);
    relaxation(params, psi, zeta, Q);
    calculateSpeeds(params, psi, zeta, u, v);
    printMatrix("psi" + filename, params, psi);
    printMatrix("zeta" + filename, params, zeta);
    printMatrix("u" + filename, params, u);
    printMatrix("v" + filename, params, v);

    Q = -4000.;
    filename = "-Q4k.txt";
    clearTables(params, psi, zeta, u, v);
    setBoundaryConditionsPsi(params, psi, Q);
    setBoundaryConditionsZeta(params, psi, zeta, Q);
    relaxation(params, psi, zeta, Q);
    calculateSpeeds(params, psi, zeta, u, v);
    printMatrix("psi" + filename, params, psi);
    printMatrix("zeta" + filename, params, zeta);
    printMatrix("u" + filename, params, u);
    printMatrix("v" + filename, params, v);

    Q = 4000.;
    filename = "+Q4k.txt";
    clearTables(params, psi, zeta, u, v);
    setBoundaryConditionsPsi(params, psi, Q);
    setBoundaryConditionsZeta(params, psi, zeta, Q);
    relaxation(params, psi, zeta, Q);
    calculateSpeeds(params, psi, zeta, u, v);
    printMatrix("psi" + filename, params, psi);
    printMatrix("zeta" + filename, params, zeta);
    printMatrix("u" + filename, params, u);
    printMatrix("v" + filename, params, v);

    return 0;
}

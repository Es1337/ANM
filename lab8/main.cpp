#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

#define OUT

using namespace std;


struct DataLine
{
    DataLine() {};
    DataLine(int _i, int _j, double _psi) : i{_i}, j{_j}, psi{_psi} {};
    
    int i, j;
    double psi;
};

ostream& operator<<(ostream& os, DataLine dl)
{
    os << dl.i << ", " << dl.j << ", " << dl.psi;
    return os;
}

struct Params
{
    Params(int _nx, int _ny, int _i1, int _i2, int _j1, double _delta) : nx{_nx}, ny{_ny}, i1{_i1}, i2{_i1}, j1{_j1}, delta{_delta} {};
    int nx;
    int ny;
    int i1;
    int i2;
    int j1;
    double delta;
    double sigma = 10*delta;
    double xa = 0.45;
    double ya = 0.45;
    double PI = 3.14159265359;
};

void printMatrix(const string& filename, const Params& params, vector<vector<double>>& matrix) 
{
    fstream fs;
    fs.open(filename, fstream::out);
    for(int i = 0; i < params.nx + 1; i++){
        for(int j = 0; j < params.ny + 1; j++){
            fs << params.delta*i << ", " << params.delta*j << ", " << matrix[i][j] << endl;
        }
    }
    fs.close();
}

void setupVelocities(const Params& p, vector<vector<double>>& vx, vector<vector<double>>& vy, vector<DataLine>& psiV)
{
    for(int i = 1; i < p.nx; i++)
    {
        for(int j = 1; j < p.ny; j++)
        {
            vx[i][j] = (psiV[i*91+j+1].psi - psiV[i*91+j-1].psi)/(2.*p.delta);
            vy[i][j] = (psiV[(i+1)*91+j].psi - psiV[(i-1)*91+j-1].psi)/(2.*p.delta);
        }
    }

    for(int i = p.i1; i < p.i2+1; i++)
    {
        for(int j = 0; j < p.j1+1; j++)
        {
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }

    for(int i = 0; i < p.nx; i++)
    {
        vx[i][0] = 0.0;
        vy[i][p.ny] = 0.0;
    }

    for(int j = 0; j < p.ny+1; j++)
    {
        vx[0][j] = vx[1][j];
        vx[p.nx][j] = vx[p.nx-1][j];
    }
}

double maxV(const Params& p, vector<vector<double>>& vx, vector<vector<double>>& vy)
{
    double max = 0.0;
    double tmp = 0.0;

    for(int i = 0; i < p.nx; i++)
    {
        for(int j = 0; j < p.ny; j++)
        {
            tmp = sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2));
            if(tmp > max)
                max = tmp;
        }
    }

    return max;
}

// double crankNicolson(const Params& p, vector<vector<double>>& u0, vector<vector<double>>& u1, double mu, int i, int j, double dt)
// {
//     double V = maxV(p, vx, vy);
//     double dt = p.delta / 4.0*;
//     double result = 0.0;
//     result += u0[i][j];
//     result -= 
// }

void AD(const Params& p, double D, vector<DataLine>& psiV, vector<vector<double>>& vx, vector<vector<double>>& vy)
{
    vector<vector<double>> u0, u1;
    u0.resize(p.nx+1);
    for(auto& item : u0)
        item.resize(p.ny, 0.0);

    u1.resize(p.nx+1);
    for(auto& item : u1)
        item.resize(p.ny, 0.0);

    for(int i = 0; i < p.nx+1; i++)
    {
        double x = i*p.delta;
        for(int j = 0; j < p.ny; j++)
        {
            double y = j*p.delta;
            u0[i][j] = (1/(2*p.PI*p.sigma*p.sigma))*exp(-(pow(x-p.xa, 2) + pow(y - p.ya, 2))/(2*p.sigma*p.sigma));
        }
    }

    for(int it = 1; it <= 60; it++)
    {
        u1 = u0;
        for(int k = 1; k <= 20; k++)
        {
            for(int i = 0; i < p.nx+1; i++)
            {
                for(int j = 1; j <p.ny; j++)
                {
                    if(i >= p.i1 && i <= p.i2 && j <= p.j1)
                    {
                        continue;
                    }
                    else if(i == 0 || i == p.nx)
                    {
                        
                    }
                    else 
                    {

                    }
                }
            }
        }
        u0 = u1;
    }
}

int main()
{
    Params params{400, 90, 200, 210, 50, 0.01};

    string ls;
    ifstream is("psi.dat");
    double a, b, c;
    vector<string> tmp;
    vector<DataLine> psiV;
    vector<vector<double>> vx ,vy; 
    psiV.resize(401*91); 

    vx.resize(params.nx + 1);
    for(auto& item : vx)
        item.resize(params.ny + 1, 0.0);

    vy.resize(params.nx + 1);
    for(auto& item : vy)
        item.resize(params.ny + 1, 0.0);
    
    for(DataLine& line : psiV)
    {
        is >> a >> b >> c;
        line.i = a;
        line.j = b;
        line.psi = c;
        // cout << line << "\n";
    }
    is.close();

    psiV.insert(psiV.begin(), *psiV.end());
    psiV.insert(psiV.end(), *(psiV.begin()+1));

    setupVelocities(params, vx, vy, psiV);

    printMatrix("vx.txt", params, vx);
    printMatrix("vy.txt", params, vy);


    return 0;
}
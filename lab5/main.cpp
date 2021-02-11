// Compile with -O2 for fastest execution
#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>


using namespace std;

vector<int> ks{16, 8, 4, 2, 1};
const double delta = 0.2;
const int n_x = 128;
const int n_y = 128;
const double pi = M_PI;

double S(vector<vector<double>> &V, const int& k)
{
    double sum = 0.0;

    for(int i = 0; i <= n_x-k; i+=k)
    {
        for(int j = 0; j <= n_y-k; j+=k)
        {
            sum += 0.5*pow(k*delta,2) * (pow(((V[i+k][j] - V[i][j]) + (V[i+k][j+k] - V[i][j+k]))/(2*k*delta), 2) 
                                    + pow(((V[i][j+k] - V[i][j]) + (V[i+k][j+k] - V[i+k][j]))/(2*k*delta), 2));
        }
    }

    return sum;
}

int main()
{
    fstream fs;
    double tol = 1e-8;

    string filename;
    vector<vector<double>> V(n_x + 1);

    double sum = 0.;
    int iteration = 0;
    double prev_sum = 1.;

    for(int i = 0; i < n_x + 1; i++)
    {
        V[i].resize(n_y + 1, 0.);

        for(int j = 0; j < n_y + 1; j++)
        {
            V[i][j] = 0.;
        }
    }

    for(int i = 0; i < n_x + 1; i++)
    {
        V[i][0] = sin(2*pi * i / n_x );
        V[i][n_y] = -sin(2*pi * i / n_x );

    }
    for(int j = 0; j < n_y + 1; j++)
    {
        V[0][j] = sin(pi * j / n_y );
        V[n_x][j] = sin(pi * j / n_y );
    }

    for(const auto& k : ks)
    {
        sum = 0.0;
        filename = "rel_" + to_string(k) + "_";
        fs.open(filename + "s.txt", fstream::out);
        do{
            for(int i = k; i <= n_x-k; i+=k)
            {
                for(int j = k; j <= n_y-k; j+=k)
                {
                    V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
                }
            }

            prev_sum = sum;
            sum = S(V, k);
            iteration++;

            fs << iteration << ", " << sum << "\n";
        
        }while(fabs((sum - prev_sum)/prev_sum) > tol);
        fs.close();

        fs.open(filename + "v.txt", fstream::out);
        for(int i = 0; i <= n_x-k; i+=k)
        {
            for(int j = 0; j <= n_y-k; j+=k)
            {
                fs << delta*i << ", " << delta*j << ", " << V[i][j] << endl;
            }
        }
        fs.close();

        if(k != 1)
        {
            for(int i = 0; i < n_x-k; i+=k)
            {   
                for(int j = 0; j < n_y-k; j+=k)
                {
                    V[i+k/2][j+k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
                    V[i+k][j+k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]);
                    V[i+k/2][j+k] = 0.5*(V[i][j+k] + V[i+k][j+k]);
                    V[i+k/2][j] = 0.5*(V[i][j] + V[i+k][j]);
                    V[i][j+k/2] = 0.5*(V[i][j] + V[i][j+k]);
                }
        }
        }
    }

    return 0;
}
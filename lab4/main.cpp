// Compile with -O2 for fastest execution
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

const double delta = 0.1;
const int n_x = 150;
const int n_y = 100;

double ro(double x, double y)
{
    double sigma_x_sqr = pow(0.1 * n_x * delta, 2);
    double sigma_y_sqr = pow(0.1 * n_y * delta, 2);

    double numerator_x = pow(x - 0.35 * n_x * delta, 2);
    double numerator_y = pow(y - 0.5 * n_y * delta, 2);

    double ro_1 = exp(-(numerator_x/sigma_x_sqr) - (numerator_y/sigma_y_sqr));

    numerator_x = pow(x - 0.65 * n_x * delta, 2);

    double ro_2 = exp(-(numerator_x/sigma_x_sqr) - (numerator_y/sigma_y_sqr));

    return ro_1 - ro_2;
}

double S(vector<vector<double>> &V)
{
    double sum = 0.0;

    for(int i = 0; i < n_x; i++)
    {
        for(int j = 0; j < n_y; j++)
        {
            sum += pow(delta,2) * (0.5 * pow((V[i+1][j] - V[i][j])/delta, 2) 
                                  + 0.5 * pow((V[i][j+1] - V[i][j])/delta, 2)
                                  - (ro(i*delta, j*delta)) * V[i][j]);
        }
    }

    return sum;
}

int main()
{
    fstream fs;
    vector<double> omega_glob{0.6, 1.0};
    vector<double> omega_loc{1.0, 1.4, 1.8, 1.9};
    double epsilon = 1.;
    double delta_sq = 0.01;
    double v_1 = 10.;

    double tol = 1e-8;

    string filename;

///////////////////// Global relaxation ///////////////////////
    for(auto& omega : omega_glob)
    {
        vector<vector<double>> V(n_x + 1);
        vector<vector<double>> V_prev(n_x + 1);
        vector<vector<double>> ro_v(n_x + 1);
        double sum = 0.;
        int iteration = 0;
        double prev_sum = 0.;
        filename = "glob_" + to_string(omega) + "_";
        fs.open(filename + "s.txt", fstream::out);
        #pragma omp parallel for
        for(int i = 0; i < n_x + 1; i++)
        {
            V[i].resize(n_y + 1, 0.);
            V_prev[i].resize(n_y + 1, 0.);
            ro_v[i].resize(n_y + 1, 0.);

            V[i][0] = v_1;
            V_prev[i][0] = v_1;

            for(int j = 0; j < n_y + 1; j++)
            {
                ro_v[i][j] = ro(i*delta, j*delta);
            }
        }

        do{
            #pragma omp parallel for
            for(int i = 0; i < n_x+1; i++)
            {
                for(int j = 0; j < n_y+1; j++)
                {
                    V_prev[i][j] = V[i][j];
                }
            }

            #pragma omp parallel for
            for(int i = 1; i < n_x; i++)
            {
                for(int j = 1; j < n_y; j++)
                {
                    V[i][j] = 0.25*(V_prev[i+1][j] + V_prev[i-1][j] + V_prev[i][j+1] + V_prev[i][j-1] + delta_sq*ro_v[i][j]);
                }
            }

            #pragma omp parallel for
            for(int j = 1; j < n_y + 1; j++)
            {
                V[0][j] = V[1][j];
                V[n_x][j] = V[n_x - 1][j];
            }

            #pragma omp parallel for
            for(int i = 0; i < n_x+1; i++)
            {   
                #pragma omp parallel for
                for(int j = 0; j < n_y+1; j++)
                {
                    V[i][j] = (1- omega)*V_prev[i][j] + omega*V[i][j];
                }
            }

            prev_sum = sum;
            sum = S(V);
            iteration++;

            fs << iteration << ", " << sum << "\n";
            
        }while(fabs((sum - prev_sum)/prev_sum) > tol);
        fs.close();

        fs.open(filename + "err.txt", fstream::out);
        #pragma omp parallel for
        for(int i = 1; i < n_x; i++)
        {
            #pragma omp parallel for
            for(int j = 1; j < n_y; j++)
            {
                fs << ((V[i+1][j] - 2*V[i][j] + V[i-1][j])/delta_sq) + ((V[i][j+1] -2*V[i][j] + V[i][j-1])/delta_sq) + ro_v[i][j] << ", " << i*delta << ", " << j*delta << "\n";
            }
        }
        fs.close();

        fs.open(filename + "v.txt", fstream::out);
        #pragma omp parallel for
        for(int i = 1; i < n_x; i++)
        {
            #pragma omp parallel for
            for(int j = 1; j < n_y; j++)
            {
                fs << V[i][j] << ", " << i*delta << ", " << j*delta << "\n";
            }
        }
        fs.close();
    }

///////////////////// Local relaxation ///////////////////////
    for(auto& omega : omega_loc)
    {
        vector<vector<double>> V(n_x + 1);
        vector<vector<double>> V_prev(n_x + 1);
        vector<vector<double>> ro_v(n_x + 1);
        double sum = 0.;
        int iteration = 0;
        double prev_sum = 0.;

        filename = "loc_" + to_string(omega) + "_";
        fs.open(filename + "s.txt", fstream::out);
        #pragma omp parallel for
        for(int i = 0; i < n_x + 1; i++)
        {
            V[i].resize(n_y + 1, 0.);
            ro_v[i].resize(n_y + 1, 0.);

            V[i][0] = v_1;

            for(int j = 0; j < n_y + 1; j++)
            {
                ro_v[i][j] = ro(i*delta, j*delta);
            }
        }

        do{
            #pragma omp parallel for
            for(int i = 1; i < n_x; i++)
            {
                for(int j = 1; j < n_y; j++)
                {
                    V[i][j] = (1 - omega)*V[i][j] + omega*0.25*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + delta_sq*ro_v[i][j]);
                }
            }

            #pragma omp parallel for
            for(int j = 1; j < n_y + 1; j++)
            {
                V[0][j] = V[1][j];
                V[n_x][j] = V[n_x - 1][j];
            }

            prev_sum = sum;
            sum = S(V);
            iteration++;

            fs << iteration << ", " << sum << "\n";
            
        }while(fabs((sum - prev_sum)/prev_sum) > tol);
        fs.close();

        fs.open(filename + "err.txt", fstream::out);
        #pragma omp parallel for
        for(int i = 1; i < n_x; i++)
        {
            #pragma omp parallel for
            for(int j = 1; j < n_y; j++)
            {
                fs << ((V[i+1][j] - 2*V[i][j] + V[i-1][j])/delta_sq) + ((V[i][j+1] -2*V[i][j] + V[i][j-1])/delta_sq) + ro_v[i][j] << ", " << i*delta << ", " << j*delta << "\n";
            }
        }
        fs.close();

        fs.open(filename + "v.txt", fstream::out);
        #pragma omp parallel for
        for(int i = 1; i < n_x; i++)
        {
            #pragma omp parallel for
            for(int j = 1; j < n_y; j++)
            {
                fs << V[i][j] << ", " << i*delta << ", " << j*delta << "\n";
            }
        }
        fs.close();
    }

    return 0;
}
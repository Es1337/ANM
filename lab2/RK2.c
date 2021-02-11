#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t, double u)
{
    double beta = 0.001;
    double N = 500;
    double gamma = 0.1;

    return (beta*N - gamma)*u - beta*pow(u,2);
}

int main(void)
{
    FILE *file;
    double beta = 0.001;
    double N = 500;
    double gamma = 0.1;
    double tMax = 100;
    double dt = 0.1;
    double u0 = 1;
    double tol = 0.000001;
    double mi = 0;
    double n = 0;
    double alpha = beta * N - gamma;

    double u, u_mi;
    double uPrev, u_miPrev;
    double time = 0.0;
    u = u0;
    u_mi = u0;

//////////////////// IMPLICIT RK2 ////////////////////

//////////////////// Butcher tables ////////////////////

    double a[2][2] = {{1/4.,              1/4. - sqrt(3.)/6.},
                     {1/4. + sqrt(3.)/6., 1/4.}};
    double c[2] = {1/2. - sqrt(3.)/6., 1/2. + sqrt(3.)/6.};
    double b = 1/2.;

    double F1, F2;
    double U1, U2;
    double U1_mi, U2_mi;
    double U1_miPrev, U2_miPrev;
    double dU1, dU2;
    double m[2][2];

    file = fopen("ImplicitRK2.txt", "w");

    n=0;

    U1 = u0;
    U2 = u0;

    U1_mi = u0;
    U2_mi = u0;

    while(time <= tMax)
    {
        mi = 0;

        while(mi <= 20.0)
        {
            F1 = U1 - u - dt*( a[0][0] * (alpha*U1 - beta*pow(U1, 2)) + a[0][1] * (alpha*U2 - beta*pow(U2, 2)));
            F2 = U2 - u - dt*( a[1][0] * (alpha*U1 - beta*pow(U1, 2)) + a[1][1] * (alpha*U2 - beta*pow(U2, 2)));

            m[0][0] = 1 - dt*a[0][0]*(alpha - 2*beta*U1);
            m[0][1] = -dt*a[0][1]*(alpha - 2*beta*U2);
            m[1][0] = -dt*a[1][0]*(alpha - 2*beta*U1);
            m[1][1] = 1 - dt*a[1][1]*(alpha - 2*beta*U2);

            dU1 = (F2*m[0][1] - F1*m[1][1]) / (m[0][0]*m[1][1] - m[0][1]*m[1][0]);
            dU2 = (F1*m[1][0] - F2*m[0][0]) / (m[0][0]*m[1][1] - m[0][1]*m[1][0]);

            U1_miPrev = U1_mi;
            U2_miPrev = U2_mi;

            U1_mi = U1_mi + dU1;
            U2_mi = U2_mi + dU2;

            mi++;
            if(abs(dU1 < tol) || (abs(dU2 < tol)))
            {
                break;
            }
        }

        U1 = U1_mi;
        U2 = U2_mi;

        u = u + dt*(b * f(time + c[0]*dt, U1) + b * f(time + c[1]*dt, U2));
        fprintf(file, "%f, %f, %f\n", time, u, (N-u));
        n++;
        time = n*dt;
    }

    fclose(file);

    return 0;
}
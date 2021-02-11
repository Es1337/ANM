#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

//////////////////// TRAPEZOID ////////////////////

//////////////////// Newton iteration ////////////////////

    file = fopen("Newton.txt", "w");

    while(time <= tMax)
    {
        u_mi = u;
        mi = 0;

        while(mi <= 20.0)
        {
            u_miPrev = u_mi;
            u_mi = u_mi - ( u_mi - u - dt/2.*( (alpha*u - beta*pow(u, 2)) + (alpha*u_mi - beta*pow(u_mi, 2)) ) )
                                            /( 1 - dt/2.*(alpha - 2*beta*u_mi) ); 
            mi++;
            if(abs(u_mi - u_miPrev) < tol)
            {
                break;
            }
        }

        uPrev = u;
        u = u_mi;
        fprintf(file, "%f, %f, %f\n", u, (N-u), time);
        n++;
        time = n*dt;
    }

    fclose(file);
    return 0;
}
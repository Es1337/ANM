#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define OUT

double maximum(double a, double b)
{
    if(a > b)
        return a;
    else
        return b;
}

void rk2(OUT double* xPointer, OUT double* vPointer, double dt, double alpha)
{
    double x, v;
    double k1x, k1v, k2x, k2v;

    x = *xPointer;
    v = *vPointer;

    k1x = v;
    k1v = alpha * v * (1 - pow(x,2)) - x;

    double t = x + dt * k1x;
    k2x = v + dt * k1v;
    k2v = alpha * k2x * (1 - pow(t,2)) - t;

    *xPointer += (dt/2.) * (k1x + k2x);
    *vPointer += (dt/2.) * (k1v + k2v);
}

void trapezoid(OUT double* xPointer, OUT double* vPointer, double dt, double alpha)
{
    double delta = 1e-10;
    double F, G;
    double xTmp, vTmp;
    double dx, dv;
    double a[2][2];
    double x, v;
    int counter = 0;

    x = *xPointer;
    v = *vPointer;

    xTmp = *xPointer;
    vTmp = *vPointer;

    while(counter < 20)
    {
        F = xTmp - x - (dt/2.)*(v + vTmp);
        G = vTmp - v - (dt/2.)*(alpha * v * (1 - pow(x,2) - x)
                               + alpha * vTmp * (1 - pow(xTmp,2)) - xTmp);

        a[0][0] = 1;
        a[0][1] = -(dt/2.);
        a[1][0] = -(dt/2.)*(-2 * alpha * xTmp * vTmp - 1);
        a[1][1] = 1 - (dt/2.)* alpha *(1 - pow(xTmp,2));

        dx = (-F * a[1][1] + G * a[0][1])
            /(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
        dv = (F * a[1][0] - G * a[0][0])
            /(a[0][0]*a[1][1] - a[0][1]*a[1][0]);

        xTmp += dx;
        vTmp += dv;
        counter++;

        if(fabs(dx) < delta && fabs(dv) < delta)
        {
            break;
        }
    }

    *xPointer = xTmp;
    *vPointer = vTmp;
}

int main(void){
    FILE* file;
    double x0 = 0.01;
    double v0 = 0.;
    double dt0 = 1.;
    double s = 0.75;
    double tMax = 40.;
    double alpha = 5.;
    double tol1 = 1e-2;
    double tol2 = 1e-5;
    double p = 2.;

    double time = 0.;
    double dt;
    double x, v;
    double x1, v1; 
    double x2, v2;
    double ex, ev;

    // TOL = 10^-2 RK2
    file = fopen("rk2tol1.txt", "w");
    x = x0;
    v = v0;
    x1 = x0;
    v2 = v0;
    x2 = x0;
    v2 = v0;
    dt = dt0;
    while(time < tMax)
    {
        x2 = x;
        x1 = x;
        v2 = v;
        v1 = v;
        rk2(&x2, &v2, dt, alpha);
        rk2(&x2, &v2, dt, alpha);

        rk2(&x1, &v1, 2*dt, alpha);

        ex = fabs((x2 - x1)/3.);
        ev = fabs((v2 - v1)/3.);

        if(maximum(ex,ev) < tol1)
        {
            x = x2;
            v = v2;
            time += 2*dt;
            fprintf(file, "%f, %f, %f, %f\n", time, dt, x, v); 
        }

        double k = (s*tol1)/(maximum(ex,ev));
        dt =  pow(k, 1/3.) * dt;
    }
    fclose(file);

    // TOL = 10^-5 RK2
    file = fopen("rk2tol2.txt", "w");
    time = 0;
    x = x0;
    v = v0;
    x1 = x0;
    v2 = v0;
    x2 = x0;
    v2 = v0;
    dt = dt0;
    while(time < tMax)
    {
        x2 = x;
        x1 = x;
        v2 = v;
        v1 = v;
        rk2(&x2, &v2, dt, alpha);
        rk2(&x2, &v2, dt, alpha);

        rk2(&x1, &v1, 2*dt, alpha);

        ex = fabs((x2 - x1)/3.);
        ev = fabs((v2 - v1)/3.);

        if(maximum(ex,ev) < tol2)
        {
            x = x2;
            v = v2;
            time += 2*dt;
            fprintf(file, "%f, %f, %f, %f\n", time, dt, x, v); 
        }

        double k = (s*tol2)/(maximum(ex,ev));
        dt =  pow(k, 1/3.) * dt;
    }
    fclose(file);

    // TOL = 10^-2 Trapezoid
    time = 0;
    x = x0;
    v = v0;
    x1 = x0;
    v2 = v0;
    x2 = x0;
    v2 = v0;
    dt = dt0;
    file = fopen("trapeztol1.txt", "w");
    while(time < tMax)
    {
        x2 = x;
        x1 = x;
        v2 = v;
        v1 = v;
        trapezoid(&x2, &v2, dt, alpha);
        trapezoid(&x2, &v2, dt, alpha);

        trapezoid(&x1, &v1, 2*dt, alpha);

        ex = fabs((x2 - x1)/3.);
        ev = fabs((v2 - v1)/3.);

        if(maximum(ex,ev) < tol1)
        {
            x = x2;
            v = v2;
            time += 2*dt;
            fprintf(file, "%f, %f, %f, %f\n", time, dt, x, v); 
        }

        double k = (s*tol1)/(maximum(ex,ev));
        dt =  pow(k, 1/3.) * dt;
    }
    fclose(file);

    // TOL = 10^-5 Trapezoid
    file = fopen("trapeztol2.txt", "w");
    time = 0;
    x = x0;
    v = v0;
    x1 = x0;
    v2 = v0;
    x2 = x0;
    v2 = v0;
    dt = dt0;
    while(time < tMax)
    {
        x2 = x;
        x1 = x;
        v2 = v;
        v1 = v;
        trapezoid(&x2, &v2, dt, alpha);
        trapezoid(&x2, &v2, dt, alpha);

        trapezoid(&x1, &v1, 2*dt, alpha);

        ex = fabs((x2 - x1)/3.);
        ev = fabs((v2 - v1)/3.);

        if(maximum(ex,ev) < tol2)
        {
            x = x2;
            v = v2;
            time += 2*dt;
            fprintf(file, "%f, %f, %f, %f\n", time, dt, x, v); 
        }

        double k = (s*tol2)/(maximum(ex,ev));
        dt =  pow(k, 1/3.) * dt;
    }
    fclose(file);

    return 0;
}
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

struct Omega {
    private:
        double omega0 = 1 / sqrt(0.0001);
    
    public:
        double omega1 = 0.5*omega0;
        double omega2 = 0.8*omega0;
        double omega3 = omega0;
        double omega4 = 1.2*omega0;
};

double potential(double omega, double t)
{
    return 10*sin(omega*t);
}

struct Params1 {
    Params1(double _y0, double _lambda, double _maxTime, double _t1, double _t2, double _t3)
         : y0{_y0}, lambda{_lambda}, maxTime{_maxTime}, t1{_t1}, t2{_t2}, t3{_t3} {};

    double y0;
    double lambda;
    double maxTime;
    double t1;
    double t2;
    double t3;

    double e = 2.71828;
};

struct Params2 {
    Params2(double _Q0, double _I0, double _dt, double _R, double _L, double _C) 
            : Q0{_Q0}, I0{_I0}, dt{_dt}, R{_R}, L{_L}, C{_C} {};

    double dt;
    double R;
    double L;
    double C;
    double Q0;
    double I0;

    double PI = 3.14159265359;

    double omega0 = 1/sqrt(L*C);
    double T0 = 2*PI / omega0;
    double maxTime = 4 * T0;
};

void euler(string filename, Params1 p, double dt) 
{
    fstream fs;
    fs.open(filename, fstream::out);

    double time = 0.0;
    int n = 0;

    double y_0 = p.y0;
    double y;

    double analytical = pow(p.e, p.lambda*time);
    double err = y_0 - analytical;

    fs << p.y0 << ", " << err << ", " << time << endl;

    y = y_0 + dt * p.lambda * y_0;

    n++;
    time = dt * n;

    analytical = pow(p.e, p.lambda*time);
    err = y - analytical;

    fs << y << ", " << err << ", " << time << endl;

    while(time < p.maxTime)
    {
        y_0 = y;
        y = y + dt * p.lambda * y;

        n++;
        time = dt*n;

        analytical = pow(p.e, p.lambda*time);
        err = y - analytical;

        fs << y << ", " << err << ", " << time << endl;
    }

    fs.close();
}

void rk2(string filename, Params1 p, double dt) 
{
    fstream fs;
    fs.open(filename, fstream::out);

    double time = 0.0;
    int n = 0;

    double y_0 = p.y0;
    double y;
    double k1, k2;
    k1 = p.lambda * y_0;
    k2 = p.lambda * (y_0 + dt*k1);  

    double analytical = pow(p.e, p.lambda*time);
    double err = y_0 - analytical;

    fs << p.y0 << ", " << err << ", " << time << endl;

    y = y_0 + (dt/2.)*(k1 + k2);

    n++;
    time = dt * n;

    analytical = pow(p.e, p.lambda*time);
    err = y - analytical;

    fs << y << ", " << err << ", " << time << endl;

    while(time < p.maxTime)
    {
        k1 = p.lambda * y;
        k2 = p.lambda * (y + dt*k1);  

        y_0 = y;
        y = y + (dt/2.)*(k1 + k2);

        n++;
        time = dt*n;

        analytical = pow(p.e, p.lambda*time);
        err = y - analytical;

        fs << y << ", " << err << ", " << time << endl;
    }

    fs.close();
}

void rk4(string filename, Params1 p, double dt) 
{
    fstream fs;
    fs.open(filename, fstream::out);

    double time = 0.0;
    int n = 0;

    double y_0 = p.y0;
    double y;
    double k1, k2, k3, k4;
    k1 = p.lambda * y_0;
    k2 = p.lambda * (y_0 + (dt/2.)*k1);
    k3 = p.lambda * (y_0 + (dt/2.)*k2);
    k4 = p.lambda * (y_0 + dt*k3);  

    double analytical = pow(p.e, p.lambda*time);
    double err = y_0 - analytical;

    fs << p.y0 << ", " << err << ", " << time << endl;

    y = y_0 + (dt/6.)*(k1 + 2*k2 + 2*k3 + k4);

    n++;
    time = dt * n;

    analytical = pow(p.e, p.lambda*time);
    err = y - analytical;

    fs << y << ", " << err << ", " << time << endl;

    while(time < p.maxTime)
    {
        k1 = p.lambda * y;
        k2 = p.lambda * (y_0 + (dt/2.)*k1);
        k3 = p.lambda * (y_0 + (dt/2.)*k2);
        k4 = p.lambda * (y_0 + dt*k3);  

        y_0 = y;
        y = y_0 + (dt/6.)*(k1 + 2*k2 + 2*k3 + k4);

        n++;
        time = dt*n;

        analytical = pow(p.e, p.lambda*time);
        err = y - analytical;

        fs << y << ", " << err << ", " << time << endl;
    }

    fs.close();
}

void ode2(string filename, double omega, Params2 p)
{
    fstream fs;
    fs.open(filename, fstream::out);

    double k1Q, k2Q, k3Q, k4Q;
    double k1I, k2I, k3I, k4I;
    double Q, I;
    double Q0, I0;
    Q0 = p.Q0;
    I0 = p.I0;

    double Vn = 0.0;
    int n = 0;
    double time = 0.0;

    fs << Q0 << ", " << I0 << ", " << time << endl;

    k1Q = I0;
    k1I = Vn/p.L - Q0/(p.L*p.C) - (p.R/p.L)*I0;

    n++;
    time = 0.5*p.dt * n;
    Vn = potential(omega, 0.5*p.dt * n);
    k2Q = I0 + (0.5*p.dt)*k1I;
    k2I = Vn/p.L - (Q0 + (0.5*p.dt)*k1Q )/(p.L*p.C) - (p.R/p.L)*k2Q;

    k3Q = I0 + (0.5*p.dt)*k2I;
    k3I = Vn/p.L - (Q0 + (0.5*p.dt)*k2Q )/(p.L*p.C) - (p.R/p.L)*k3Q;

    n++;
    time = 0.5*p.dt*n;
    Vn = potential(omega, 0.5*p.dt*n);
    k4Q = I0 + p.dt*k3I;
    k4I = Vn/p.L - (Q0 + p.dt*k3Q )/(p.L*p.C) - (p.R/p.L)*k4Q;


    Q = Q0 + (p.dt/6.)*(k1Q + 2*k2Q + 2*k3Q + k4Q);
    I = I0 + (p.dt/6.)*(k1I + 2*k2I + 2*k3I + k4I);

    fs << Q << ", " << I << ", " << time << endl;

    while(time < p.maxTime)
    {
        k1Q = I;
        k1I = Vn/p.L - Q/(p.L*p.C) - (p.R/p.L)*I;

        n++;
        time = 0.5*p.dt * n;
        Vn = potential(omega, 0.5*p.dt * n);
        k2Q = I + (0.5*p.dt)*k1I;
        k2I = Vn/p.L - (Q + (0.5*p.dt)*k1Q )/(p.L*p.C) - (p.R/p.L)*k2Q;

        k3Q = I0 + (0.5*p.dt)*k2I;
        k3I = Vn/p.L - (Q + (0.5*p.dt)*k2Q )/(p.L*p.C) - (p.R/p.L)*k3Q;

        n++;
        time = 0.5*p.dt*n;
        Vn = potential(omega, 0.5*p.dt*n);
        k4Q = I + p.dt*k3I;
        k4I = Vn/p.L - (Q + p.dt*k3Q )/(p.L*p.C) - (p.R/p.L)*k4Q;
    
        Q0 = Q;
        I0 = I;
        Q += (p.dt/6.)*(k1Q + 2*k2Q + 2*k3Q + k4Q);
        I += (p.dt/6.)*(k1I + 2*k2I + 2*k3I + k4I);

        fs << Q << ", " << I << ", " << time << endl;
    }

    fs.close();
}

int main() {
    Omega omega;
    Params1 params1{1., -1., 5., 0.01, 0.1, 1.};
    Params2 params2{0.0, 0.0, 0.0001, 100., 0.1, 0.001};

    euler("t1_euler.txt", params1, params1.t3);
    euler("t01_euler.txt", params1, params1.t2);
    euler("t001_euler.txt", params1, params1.t1);

    rk2("t1_rk2.txt", params1, params1.t3);
    rk2("t01_rk2.txt", params1, params1.t2);
    rk2("t001_rk2.txt", params1, params1.t1);

    rk4("t1_rk4.txt", params1, params1.t3);
    rk4("t01_rk4.txt", params1, params1.t2);
    rk4("t001_rk4.txt", params1, params1.t1);

    ode2("omega1_ode2.txt", omega.omega1, params2);
    ode2("omega2_ode2.txt", omega.omega2, params2);
    ode2("omega3_ode2.txt", omega.omega3, params2);
    ode2("omega4_ode2.txt", omega.omega4, params2);

    return 0;
}
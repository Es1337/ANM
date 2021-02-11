#include <stdio.h>
#include <math.h>

double potential(double t)
{
    double omega_0 = 1 / (sqrt(0.0001));
    double omega_V1 = 0.5*omega_0;
    double omega_V2 = 0.8*omega_0;
    double omega_V3 = 1.0*omega_0;
    double omega_V4 = 1.2*omega_0;

    // SWAP OUT OMEGAS
    return 10*sin(omega_V4 * t); 
}

int main(void) {
  
    double y_0 = 1.;
    double lambda = -1.;
    double max_time = 5.;
    double t1 = 0.01;
    double t2 = 0.1;
    double t3 = 1.;
    double time = 0.;
    double y;
    int n = 0;
    double tmp;
    double analytical;
    double err;
    double e = 2.71828;

/// UNCOMMENT ONE SEGMENT AT A TIME //

///////////////////////////// Euler /////////////////////////////

    // GNUPLOT COMMANDS ! uncomment err and y printfs separately ! 
    // plot "t1_result.txt" using 2:1 with points, "t01_result.txt" using 2:1 with points, "t001_result.txt" using 2:1 with points, 2.71828**(-x) 
    // plot "err1.txt" using 2:1 with lines, "err01.txt" using 2:1 with lines, "err001.txt" using 2:1 with lines    

///////////////// dt = 0.01 /////////////////
    // printf("%f, %f\n", y_0, time);
    analytical = pow(e, lambda*time);
    err = y_0 - analytical;
    printf("%f, %f\n", err, time);

    y = y_0 + t1 * lambda * y_0;
    n++;
    time = t1 * n;
    analytical = pow(e, lambda*time);
    err = y - analytical;
    printf("%f, %f\n", err, time);
    // printf("%f, %f\n", y, time);

    while(time < max_time)
    {
        y_0 = y;
        y = y + t1 * lambda * y;
        n++;
        time = t1 * n;
        analytical = pow(e, lambda*time);
        err = y - analytical;
        printf("%f, %f\n", err, time);
        // printf("%f, %f\n", y, time);
    }

///////////////// dt = 0.1 /////////////////
    // // printf("%f, %f\n", y_0, time);
    // analytical = pow(e, lambda*time);
    // err = y_0 - analytical;
    // printf("%f, %f\n", err, time);

    // y = y_0 + t2 * lambda * y_0;
    // n++;
    // time = t2 * n;
    // analytical = pow(e, lambda*time);
    // err = y - analytical;
    // printf("%f, %f\n", err, time);
    // // printf("%f, %f\n", y, time);

    // while(time < max_time)
    // {
    //     y_0 = y;
    //     y = y + t2 * lambda * y;
    //     n++;
    //     time = t2 * n;
    //     analytical = pow(e, lambda*time);
    //     err = y - analytical;
    //     printf("%f, %f\n", err, time);
    //     // printf("%f, %f\n", y, time);
    // }

///////////////// dt = 1.0 /////////////////
    // // printf("%f, %f\n", y_0, time);
    // analytical = pow(e, lambda*time);
    // err = y_0 - analytical;
    // printf("%f, %f\n", err, time);
    // y = y_0 + t3 * lambda * y_0;
    // n++;
    // analytical = pow(e, lambda*time);
    // err = y - analytical;
    // time = t3 * n;
    // printf("%f, %f\n", err, time);
    // // printf("%f, %f\n", y, time);

    // while(time < max_time)
    // {
    //     y_0 = y;
    //     y = y + t3 * lambda * y;
    //     analytical = pow(e, lambda*time);
    //     err = y - analytical;
    //     n++;
    //     time = t3 * n;
    //     printf("%f, %f\n", err, time);
    //     // printf("%f, %f\n", y, time);
    // }

///////////////////////////// RK2 /////////////////////////////

    // GNUPLOT COMMANDS ! uncomment err and y printfs separately !
    // plot "t1_resultRK2.txt" using 2:1 with points, "t01_resultRK2.txt" using 2:1 with points, "t001_resultRK2.txt" using 2:1 with points, 2.71828**(-x)
    // plot "err1RK2.txt" using 2:1 with lines, "err01RK2.txt" using 2:1 with lines, "err001RK2.txt" using 2:1 with lines         

    double k_1;
    double k_2;

///////////////// dt = 0.01 /////////////////
    k_1 = lambda * y_0;
    k_2 = lambda * (y_0 + t1*k_1);
    // printf("%f, %f\n", y_0, time);
    analytical = pow(e, lambda*time);
    err = y_0 - analytical;
    printf("%f, %f\n", err, time);

    y = y_0 + (t1/2.)*(k_1 + k_2);
    n++;
    time = t1 * n;
    // printf("%f, %f\n", y, time);
    
    analytical = pow(e, lambda*time);
    err = y - analytical;
    printf("%f, %f\n", err, time);

    while(time < max_time)
    {
        k_1 = lambda * y;
        k_2 = lambda * (y + t1*k_1);

        y_0 = y;
        y = y + (t1/2.)*(k_1 + k_2);
        n++;
        time = t1 * n;
        // printf("%f, %f\n", y, time);

        analytical = pow(e, lambda*time);
        err = y - analytical;
        printf("%f, %f\n", err, time);
    }

///////////////// dt = 0.1 /////////////////
    // k_1 = lambda * y_0;
    // k_2 = lambda * (y_0 + t2*k_1);
    // // printf("%f, %f\n", y_0, time);
    // analytical = pow(e, lambda*time);
    // err = y_0 - analytical;
    // printf("%f, %f\n", err, time);

    // y = y_0 + (t2/2.)*(k_1 + k_2);
    // n++;
    // time = t2 * n;
    // // printf("%f, %f\n", y, time);
    
    // analytical = pow(e, lambda*time);
    // err = y - analytical;
    // printf("%f, %f\n", err, time);

    // while(time < max_time)
    // {
    //     k_1 = lambda * y;
    //     k_2 = lambda * (y + t2*k_1);

    //     y_0 = y;
    //     y = y + (t2/2.)*(k_1 + k_2);
    //     n++;
    //     time = t2 * n;
    //     // printf("%f, %f\n", y, time);

    //     analytical = pow(e, lambda*time);
    //     err = y - analytical;
    //     printf("%f, %f\n", err, time);
    // }

///////////////// dt = 1.0 /////////////////
    // k_1 = lambda * y_0;
    // k_2 = lambda * (y_0 + t3*k_1);
    // printf("%f, %f\n", y_0, time);
    // analytical = pow(e, lambda*time);
    // err = y_0 - analytical;
    // // printf("%f, %f\n", err, time);

    // y = y_0 + (t3/2.)*(k_1 + k_2);
    // n++;
    // time = t3 * n;
    // printf("%f, %f\n", y, time);
    
    // analytical = pow(e, lambda*time);
    // err = y - analytical;
    // // printf("%f, %f\n", err, time);

    // while(time < max_time)
    // {
    //     k_1 = lambda * y;
    //     k_2 = lambda * (y + t3*k_1);

    //     y_0 = y;
    //     y = y + (t3/2.)*(k_1 + k_2);
    //     n++;
    //     time = t3 * n;
    //     printf("%f, %f\n", y, time);

    //     analytical = pow(e, lambda*time);
    //     err = y - analytical;
    //     // printf("%f, %f\n", err, time);
    // }

///////////////////////////// RK4 /////////////////////////////

    // GNUPLOT COMMANDS ! uncomment err and y printfs separately !
    // plot "t1_resultRK4.txt" using 2:1 with points, "t01_resultRK4.txt" using 2:1 with points, "t001_resultRK4.txt" using 2:1 with points, 2.71828**(-x)  
    // plot "err1RK4.txt" using 2:1 with lines, "err01RK4.txt" using 2:1 with lines, "err001RK4.txt" using 2:1 with lines 

    double k_3;
    double k_4;

///////////////// dt = 0.01 /////////////////
    // k_1 = lambda * y_0;
    // k_2 = lambda * (y_0 + (t1/2.)*k_1);
    // k_3 = lambda * (y_0 + (t1/2.)*k_2);
    // k_4 = lambda * (y_0 + t1*k_3);

    // // printf("%f, %f\n", y_0, time);

    // analytical = pow(e, lambda*time);
    // err = y_0 - analytical;
    // printf("%f, %f\n", err, time);

    // y = y_0 + (t1/6.)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    // n++;
    // time = t1 * n;
    // // printf("%f, %f\n", y, time);
    
    // analytical = pow(e, lambda*time);
    // err = y - analytical;
    // printf("%f, %f\n", err, time);

    // while(time < max_time)
    // {
    //     k_1 = lambda * y;
    //     k_2 = lambda * (y + (t1/2.)*k_1);
    //     k_3 = lambda * (y + (t1/2.)*k_2);
    //     k_4 = lambda * (y + t1*k_3);

    //     y_0 = y;
    //     y = y + (t1/6.)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    //     n++;
    //     time = t1 * n;
    //     // printf("%f, %f\n", y, time);

    //     analytical = pow(e, lambda*time);
    //     err = y - analytical;
    //     printf("%f, %f\n", err, time);
    // }

///////////////// dt = 0.1 /////////////////
    // k_1 = lambda * y_0;
    // k_2 = lambda * (y_0 + (t2/2.)*k_1);
    // k_3 = lambda * (y_0 + (t2/2.)*k_2);
    // k_4 = lambda * (y_0 + t2*k_3);

    // printf("%f, %f\n", y_0, time);

    // analytical = pow(e, lambda*time);
    // err = y_0 - analytical;
    // // printf("%f, %f\n", err, time);

    // y = y_0 + (t2/6.)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    // n++;
    // time = t2 * n;
    // printf("%f, %f\n", y, time);
    
    // analytical = pow(e, lambda*time);
    // err = y - analytical;
    // // printf("%f, %f\n", err, time);

    // while(time < max_time)
    // {
    //     k_1 = lambda * y;
    //     k_2 = lambda * (y + (t2/2.)*k_1);
    //     k_3 = lambda * (y + (t2/2.)*k_2);
    //     k_4 = lambda * (y + t2*k_3);

    //     y_0 = y;
    //     y = y + (t2/6.)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    //     n++;
    //     time = t2 * n;
    //     printf("%f, %f\n", y, time);

    //     analytical = pow(e, lambda*time);
    //     err = y - analytical;
    //     // printf("%f, %f\n", err, time);
    // }

///////////////// dt = 1.0 /////////////////
    // k_1 = lambda * y_0;
    // k_2 = lambda * (y_0 + (t3/2.)*k_1);
    // k_3 = lambda * (y_0 + (t3/2.)*k_2);
    // k_4 = lambda * (y_0 + t3*k_3);

    // printf("%f, %f\n", y_0, time);

    // analytical = pow(e, lambda*time);
    // err = y_0 - analytical;
    // // printf("%f, %f\n", err, time);

    // y = y_0 + (t3/6.)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    // n++;
    // time = t3 * n;
    // printf("%f, %f\n", y, time);
    
    // analytical = pow(e, lambda*time);
    // err = y - analytical;
    // // printf("%f, %f\n", err, time);

    // while(time < max_time)
    // {
    //     k_1 = lambda * y;
    //     k_2 = lambda * (y + (t3/2.)*k_1);
    //     k_3 = lambda * (y + (t3/2.)*k_2);
    //     k_4 = lambda * (y + t3*k_3);

    //     y_0 = y;
    //     y = y + (t3/6.)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    //     n++;
    //     time = t3 * n;
    //     printf("%f, %f\n", y, time);

    //     analytical = pow(e, lambda*time);
    //     err = y - analytical;
    //     // printf("%f, %f\n", err, time);
    // }

///////////////////////////// ODE 2nd degree /////////////////////////////

    // GNUPLOT COMMANDS
    // plot "omega1.txt" using 3:1 with lines, "omega2.txt" using 3:1 with lines, "omega3.txt" using 3:1 with lines, "omega4.txt" using 3:1 with lines
    // plot "omega1.txt" using 3:2  with lines, "omega2.txt" using 3:2 with lines, "omega3.txt" using 3:2 with lines, "omega4.txt" using 3:2 with lines 

    // double PI = 3.14159265359;
    // double dt = 0.0001/2.;
    // double R = 100.;
    // double L = 0.1;
    // double C = 0.001;

    // double omega_0 = 1 / (sqrt(L*C));
    // double T_0 = 2*PI / omega_0;
    // max_time = 4*T_0;

    // double Q_0 = 0;
    // double I_0 = 0;

    // double k_1Q, k_2Q, k_3Q, k_4Q;
    // double k_1I, k_2I, k_3I, k_4I;
    // double Q, I;

    // double Vn = 0;

    // n = 0;
    // time = 0;
    // printf("%f, %f, %f\n", Q_0, I_0, time);

    // k_1Q = I_0;
    // k_1I = Vn/L - Q_0/(L*C) - (R/L)*I_0;

    // n++;
    // time = dt*n;
    // Vn = potential(dt*n);
    // k_2Q = I_0 + (dt)*k_1I;
    // k_2I = Vn/L - (Q_0 + (dt)*k_1Q )/(L*C) - (R/L)*k_2Q;

    // k_3Q = I_0 + (dt)*k_2I;
    // k_3I = Vn/L - (Q_0 + (dt)*k_2Q )/(L*C) - (R/L)*k_3Q;

    // n++;
    // time = dt*n;
    // Vn = potential(dt*n);
    // k_4Q = I_0 + 2*dt*k_3I;
    // k_4I = Vn/L - (Q_0 + 2*dt*k_3Q )/(L*C) - (R/L)*k_4Q;


    // Q = Q_0 + (dt/3.)*(k_1Q + 2*k_2Q + 2*k_3Q + k_4Q);
    // I = I_0 + (dt/3.)*(k_1I + 2*k_2I + 2*k_3I + k_4I);
    // printf("%f, %f, %f\n", Q, I, time);
    
    // while(time < max_time)
    // {
    //     k_1Q = I;
    //     k_1I = Vn/L - Q/(L*C) - (R/L)*I;

    //     n++;
    //     time = dt*n;
    //     Vn = potential(dt*n);
    //     k_2Q = I + (dt)*k_1I;
    //     k_2I = Vn/L - (Q + (dt)*k_1Q )/(L*C) - (R/L)*k_2Q;

    //     k_3Q = I + (dt)*k_2I;
    //     k_3I = Vn/L - (Q + (dt)*k_2Q )/(L*C) - (R/L)*k_3Q;

    //     n++;
    //     time = dt*n;
    //     Vn = potential(dt*n);
    //     k_4Q = I + 2*dt*k_3I;
    //     k_4I = Vn/L - (Q + 2*dt*k_3Q )/(L*C) - (R/L)*k_4Q;

    //     Q_0 = Q;
    //     I_0 = I;
    //     Q = Q + (dt/3.)*(k_1Q + 2*k_2Q + 2*k_3Q + k_4Q);
    //     I = I + (dt/3.)*(k_1I + 2*k_2I + 2*k_3I + k_4I);
    //     printf("%f, %f, %f\n", Q, I, time);
    // }

    return 0;
}
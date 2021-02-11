reset
set term jpeg

set xrange [-0.1:5.1]
set yrange [-0.1:1.1]
set ylabel "y(t)"
set xlabel "time"
set out "euler.jpg"
set title "Metoda Eulera - rozw"
plot "t1_euler.txt" using 3:1 with points title "dt = 1", "t01_euler.txt" using 3:1 with points title "dt = 0.1", "t001_euler.txt" using 3:1 with points title "dt = 0.01", 2.71828**(-x) title "analityczne"

set out "rk2.jpg"
set title "Metoda RK2 - rozw"
plot "t1_rk2.txt" using 3:1 with points title "dt = 1", "t01_rk2.txt" using 3:1 with points title "dt = 0.1", "t001_rk2.txt" using 3:1 with points title "dt = 0.01", 2.71828**(-x) title "analityczne"

set out "rk4.jpg"
set title "Metoda RK4 - rozw"
plot "t1_rk4.txt" using 3:1 with points title "dt = 1", "t01_rk4.txt" using 3:1 with points title "dt = 0.1", "t001_rk4.txt" using 3:1 with points title "dt = 0.01", 2.71828**(-x) title "analityczne"

set yrange [-0.4:0.05]
set ylabel "err(t)"

set out "euler_err.jpg"
set title "Metoda Eulera - błąd"
plot "t1_euler.txt" using 3:2 with lines title "dt = 1", "t01_euler.txt" using 3:2 with lines title "dt = 0.1", "t001_euler.txt" using 3:2 with lines title "dt = 0.01"

set yrange [-0.005 : 0.14]
set out "rk2_err.jpg"
set title "Metoda RK2 - błąd"
plot "t1_rk2.txt" using 3:2 with lines title "dt = 1", "t01_rk2.txt" using 3:2 with lines title "dt = 0.1", "t001_rk2.txt" using 3:2 with lines title "dt = 0.01"

set yrange [-0.5 : 0.008]
set out "rk4_err.jpg"
set title "Metoda RK4 - błąd"
plot "t1_rk4.txt" using 3:2 with lines title "dt = 1", "t01_rk4.txt" using 3:2 with lines title "dt = 0.1", "t001_rk4.txt" using 3:2 with lines title "dt = 0.01"

set xrange [-0.01:0.26]
set yrange [-0.002:0.0035]
set xlabel "time"
set ylabel "Q(t)"
set out "Q.jpg"
set title "Q(t) Metoda RK4"
plot "omega1_ode2.txt" using 3:1 with lines title "0.5ω_0", "omega2_ode2.txt" using 3:1 with lines title "0.8ω_0", "omega3_ode2.txt" using 3:1 with lines title "ω_0", "omega4_ode2.txt" using 3:1 with lines title "1.2 ω_0"

set yrange [-0.11:0.11]
set ylabel "I(t)"
set out "I.jpg"
set title "I(t) Metoda RK4"
plot "omega1_ode2.txt" using 3:2 with lines title "0.5ω_0", "omega2_ode2.txt" using 3:2 with lines title "0.8ω_0", "omega3_ode2.txt" using 3:2 with lines title "ω_0", "omega4_ode2.txt" using 3:2 with lines title "1.2 ω_0"
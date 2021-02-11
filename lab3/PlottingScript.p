set terminal png
set output 'trapez_dt_t.png'
set xlabel "time"
set ylabel "dt(t)"
plot "trapeztol1.txt" using 1:2 with lines title "tol = 10^(-2)", "trapeztol2.txt" using 1:2 with lines title "tol = 10^(-5)"
set terminal png
set output 'trapez_x_t.png'
set xlabel "time"
set ylabel "x(t)"
plot "trapeztol1.txt" using 1:3 with lines title "tol = 10^(-2)", "trapeztol2.txt" using 1:3 with lines title "tol = 10^(-5)"
set terminal png
set output 'trapez_v_t.png'
set xlabel "time"
set ylabel "v(t)"
plot "trapeztol1.txt" using 1:4 with lines title "tol = 10^(-2)", "trapeztol2.txt" using 1:4 with lines title "tol = 10^(-5)"
set terminal png
set output 'trapez_v_x.png'
set xlabel "x"
set ylabel "v(x)"
plot "trapeztol1.txt" using 3:4 with lines title "tol = 10^(-2)", "trapeztol2.txt" using 3:4 with lines title "tol = 10^(-5)"

set terminal png
set output 'rk2_dt_t.png'
set xlabel "time"
set ylabel "dt(t)"
plot "rk2tol1.txt" using 1:2 with lines title "tol = 10^(-2)", "rk2tol2.txt" using 1:2 with lines title "tol = 10^(-5)"
set terminal png
set output 'rk2_x_t.png'
set xlabel "time"
set ylabel "x(t)"
plot "rk2tol1.txt" using 1:3 with lines title "tol = 10^(-2)", "rk2tol2.txt" using 1:3 with lines title "tol = 10^(-5)"
set terminal png
set output 'rk2_v_t.png'
set xlabel "time"
set ylabel "v(t)"
plot "rk2tol1.txt" using 1:4 with lines title "tol = 10^(-2)", "rk2tol2.txt" using 1:4 with lines title "tol = 10^(-5)"
set terminal png
set output 'rk2_v_x.png'
set xlabel "x"
set ylabel "v(x)"
plot "rk2tol1.txt" using 3:4 with lines title "tol = 10^(-2)", "rk2tol2.txt" using 3:4 with lines title "tol = 10^(-5)" 
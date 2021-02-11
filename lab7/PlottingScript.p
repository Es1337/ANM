reset
set encoding utf8
set yrange [-200:200]
set xlabel "Iterations"
set ylabel ' Γ '
set term jpeg
set out "err.jpeg"
plot "-1000.000000err.txt" u 1:2 w l t "err Q = -1000", "-4000.000000err.txt" u 1:2 w l t "err Q=-4000", "4000.000000err.txt" u 1:2 w l t "err Q=4000"

reset
set xlabel "x"
set ylabel "y"
set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )
set title "Q=-1000, u(x,y)"
set cbrange [-2:16]
set out "map_u-Q1k.jpeg"
p "u-Q1k.txt" u 1:2:3 w p pt 7 palette

set title "Q=-1000, v(x,y)"
set cbrange [-6:1]
set out "map_v-Q1k.jpeg"
p "v-Q1k.txt" u 1:2:3 w p pt 7 palette

set title "Q=-4000, u(x,y)"
set cbrange [-10:70]
set out "map_u-Q4k.jpeg"
p "u-Q4k.txt" u 1:2:3 w p pt 7 palette

set title "Q=-4000, v(x,y)"
set cbrange [-14:4]
set out "map_v-Q4k.jpeg"
p "v-Q4k.txt" u 1:2:3 w p pt 7 palette

set title "Q=4000, u(x,y)"
set cbrange [-70:10]
set out "map_u+Q4k.jpeg"
p "u+Q4k.txt" u 1:2:3 w p pt 7 palette

set title "Q=4000, v(x,y)"
set cbrange [-5:35]
set out "map_v+Q4k.jpeg"
p "v+Q4k.txt" u 1:2:3 w p pt 7 palette

reset
set terminal postscript eps enhanced
set xlabel "x"
set ylabel "y"
set palette defined (0 "violet", 1 "red", 2 "yellow")

set title "Q=-1000, ψ(x,y)"
set cbrange [-55:-50]
set term jpeg
set out "map_psi-Q1k.jpeg"
p "psi-Q1k.txt" u 1:2:3 w p pt 7 palette

set terminal postscript eps enhanced 
set title "Q=-4000, ψ(x,y)"
set cbrange [-218:-202]
set term jpeg
set out "map_zeta-Q1k.jpeg"
p "psi-Q4k.txt" u 1:2:3 w p pt 7 palette

set terminal postscript eps enhanced 
set title "Q=4000, ψ(x,y)"
set cbrange [202:218]
set term jpeg
set out "map_psi-Q4k.jpeg"
p "psi+Q4k.txt" u 1:2:3 w p pt 7 palette

set terminal postscript eps enhanced 
set title "Q=-1000, ζ(x,y)"
set cbrange [-200:350]
set term jpeg
set out "map_zeta-Q4k.jpeg"
p "zeta-Q1k.txt" u 1:2:3 w p pt 7 palette

set terminal postscript eps enhanced 
set title "Q=-4000, ζ(x,y)"
set cbrange [-800:1350]
set term jpeg
set out "map_psi+Q4k.jpeg"
p "zeta-Q4k.txt" u 1:2:3 w p pt 7 palette

set terminal postscript eps enhanced 
set title "Q=4000, ζ(x,y)"
set cbrange [-800:1000]
set term jpeg
set out "map_zeta+Q4k.jpeg"
p "zeta+Q4k.txt" u 1:2:3 w p pt 7 palette
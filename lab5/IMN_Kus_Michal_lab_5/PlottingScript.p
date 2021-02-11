set term jpeg

set title "S(it)"
set xlabel "it"
set ylabel "S"
set out "rel_s.jpeg"
plot "rel_16_s.txt" u 1:2 w l, "rel_8_s.txt" u 1:2 w l, "rel_4_s.txt" u 1:2 w l, "rel_2_s.txt" u 1:2 w l, "rel_1_s.txt" u 1:2 w l

set title "v(k=16)"
set xlabel "x"
set ylabel "y"
set out "rel_v_16.jpeg"
plot "rel_16_v.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "v(k=8)"
set xlabel "x"
set ylabel "y"
set out "rel_v_8.jpeg"
plot "rel_8_v.txt" u 1:2:3 w p pt 7 palette t 'V'  

set title "v(k=4)"
set xlabel "x"
set ylabel "y"
set out "rel_v_4.jpeg"
plot "rel_4_v.txt" u 1:2:3 w p pt 7 palette t 'V'  

set title "v(k=2)"
set xlabel "x"
set ylabel "y"
set out "rel_v_2.jpeg"
plot "rel_2_v.txt" u 1:2:3 w p pt 7 palette t 'V'  

set title "v(k=1)"
set xlabel "x"
set ylabel "y"
set out "rel_v_1.jpeg"
plot "rel_1_v.txt" u 1:2:3 w p pt 7 palette t 'V'  
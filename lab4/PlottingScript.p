set term jpeg

set logscale x
set xlabel "nr iteracji"
set ylabel "S"

set title "Relaksacja globalna - S"
set out "glob_s.jpeg"
plot "glob_0.600000_s.txt" u 1:2 w l title "{/Symbol w} = 0.6", "glob_1.000000_s.txt" u 1:2 w l title "{/Symbol w} = 1.0"

set title "Relaksacja lokalna - S"
set out "loc_s.jpeg"
plot "loc_1.000000_s.txt" u 1:2 w l title "{/Symbol w} = 1.0", "loc_1.400000_s.txt" u 1:2 w l title "{/Symbol w} = 1.4", "loc_1.800000_s.txt" u 1:2 w l title "{/Symbol w} = 1.8", "loc_1.900000_s.txt" u 1:2 w l title "{/Symbol w} = 1.9"

unset logscale x
set xlabel "X"
set ylabel "Y"

set title "Relaksacja globalna - V"
set out "glob_v_0.6.jpeg"
plot "glob_0.600000_v.txt" u 2:3:1 w p pt 7 palette t 'V'

set out "glob_v_1.0.jpeg"
plot "glob_1.000000_v.txt" u 2:3:1 w p pt 7 palette t 'V'

set title "Relaksacja globalna - err"
set out "glob_err_0.6.jpeg"
plot "glob_0.600000_err.txt" u 2:3:1 w p pt 7 palette t 'err'

set out "glob_err_1.0.jpeg"
plot "glob_1.000000_err.txt" u 2:3:1 w p pt 7 palette t 'err'
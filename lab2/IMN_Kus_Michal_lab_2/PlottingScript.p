set xlabel "time"
set ylabel "u(t), z(t)"
plot "Picard.txt" using 3:1 title "u(t)" with lines, "" using 3:2 title "z(t)" with lines 
set xlabel "time"
set ylabel "u(t), z(t)"  
plot "Newton.txt" using 3:1 title "u(t)" with lines, "" using 3:2 title "z(t)" with lines 
set xlabel "time"
set ylabel "u(t), z(t)"  
plot "ImplicitRK2.txt" using 1:2 title "u(t)" with lines, "" using 1:3 title "z(t)" with lines
set term jpeg
set cbrange[-10:10]
set palette defined(0 "red", 1 "white", 2 "blue")

set title "Map nx=ny=4"
set xlabel "x"
set ylabel "y"
set out "map_4.jpeg"
plot "map_4.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "Map nx=ny=50"
set xlabel "x"
set ylabel "y"
set out "map_50.jpeg"
plot "map_50.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "Map nx=ny=100"
set xlabel "x"
set ylabel "y"
set out "map_100.jpeg"
plot "map_100.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "Map nx=ny=200"
set xlabel "x"
set ylabel "y"
set out "map_200.jpeg"
plot "map_200.txt" u 1:2:3 w p pt 7 palette t 'V'

set cbrange[-0.8: 0.8]

set title 'Map nx=ny=100 {/Symbol e}2 = 1'
set xlabel "x"
set ylabel "y"
set out "map_100_ro1.jpeg"
plot "map_100_ro1.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "Map nx=ny=100 {/Symbol e}2 = 2"
set xlabel "x"
set ylabel "y"
set out "map_100_ro2.jpeg"
plot "map_100_ro2.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "Map nx=ny=100 {/Symbol e}2 = 10"
set xlabel "x"
set ylabel "y"
set out "map_100_ro10.jpeg"
plot "map_100_ro10.txt" u 1:2:3 w p pt 7 palette t 'V'
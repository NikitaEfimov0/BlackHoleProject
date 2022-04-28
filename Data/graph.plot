#! /usr/bin/gnuplot -persist
    set terminal postscript  enhanced color solid
    set terminal pdf
    set output "data.pdf"
    #set xrange [-100:100]
    #set yrange [-100:100]
    #set zrange [-100:100]
    set grid

    f(x)=x
splot 'S2.dat' using 2:3:4 w l, 'S55.dat' using 2:3:4 w l, 'S38.dat' using 2:3:4 w l, 'BlackHole.dat' u 1:2:3 w p ps 0.1
plot 'S2.dat' using 3:4 w l, 'S55.dat' using 3:4 w l, 'S38.dat' using 3:4 w l, 'BlackHole.dat' u 1:2:3 w p ps 0.1
plot 'S2.dat' using 2:4 w l, 'S55.dat' using 2:4 w l, 'S38.dat' using 2:4 w l, 'BlackHole.dat' u 1:2:3 w p ps 0.1
plot 'S2.dat' using 2:3 w l, 'S55.dat' using 2:3 w l, 'S38.dat' using 2:3 w l, 'BlackHole.dat' u 1:2:3 w p ps 0.1

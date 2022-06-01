#! /usr/bin/gnuplot -persist
    set terminal postscript  enhanced color solid
    set terminal pdf
    set output "data.pdf"
    set xrange [] reverse
    #set yrange [] reverse
    #set zrange [-100:100]
    set ylabel "Δ‎Dec. ('')"
    set xlabel "Δ‎R.A.('')"
    set grid

#splot 'S2.dat' using 2:3:4 w l, 'S55.dat' using 2:3:4 w l, 'S38.dat' using 2:3:4 w l, 'BlackHole.dat' u 1:2:3 w p ps 0.1

plot 'S2.dat' using 2:3 w l, 'S55.dat' using 2:3 w l, 'S38.dat' using 2:3 w l, 'BlackHole.dat' u 1:2:3 w p ps 0.1, 'S2Original.dat' u 2:3 w p ps 0.1, 'S38Original.dat' u 1:2 w p ps 0.1, 'S55Original.dat' u 1:2 w p ps 0.1

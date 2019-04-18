# Elegimos el tipo de archivo de salida, letra, y tamaño de la letra
set terminal postscript eps color enhanced "Helvetica-bold" 16
set output 'EvsZ.eps'

set xlabel 'Tip-surface distance (Ang)'
# Z limit of the plot 
# set xrange [0:* ]
set xrange [0:5 ]

set ylabel 'Energy (eV)'
# set ylabel '{/Symbol \141} (10^5/K)'
# set yrange [0:2.6]

# set format x "%.3f"
set format y "%.2f"

# Contrlamos la leyenda de lo que es cada linea de la grafica
# unset key
set key bottom right  font "Helvetica,12"
# set tics font "Helvetica,14" offset -0.25,-0.25,0

set grid

# Elegimos que queremos plotear y como (tipo de linea, color....)
set style line 1 lt 1 lc rgb 'blue' pt 2
set style line 2 lt 1 lc rgb 'red' pt 6

stats 'grid_point_1/EvsZ.dat' using 2
plot 'grid_point_1/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 1'  ,\
     'grid_point_2/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 2'  ,\
     'grid_point_3/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 3'  ,\
     'grid_point_4/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 4'  ,\
     'grid_point_5/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 5'  ,\
     'grid_point_6/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 6'  ,\
     'grid_point_7/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 7'  ,\
     'grid_point_8/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 8'  ,\
     'grid_point_9/EvsZ.dat'  u 1:($2-STATS_max) w linespoints   t 'GP 9'  ,\
     'grid_point_10/EvsZ.dat' u 1:($2-STATS_max) w linespoints   t 'GP 10' 

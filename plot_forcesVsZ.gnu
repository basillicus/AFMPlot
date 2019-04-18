# Elegimos el tipo de archivo de salida, letra, y tamaño de la letra
set terminal postscript eps color enhanced "Helvetica-bold" 16
set output 'FvsZ.eps'

set xlabel 'Tip-surface distance (Ang)'
# Z limit of the plot 
#set xrange [0:* ]
set xrange [0:5 ]

set ylabel 'Force (eV/Ang)'
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

plot 'forces_GP_1.dat'  u 1:2 w linespoints  t 'GP 1'  ,\
     'forces_GP_2.dat'  u 1:2 w linespoints  t 'GP 2'  ,\
     'forces_GP_3.dat'  u 1:2 w linespoints  t 'GP 3'  ,\
     'forces_GP_4.dat'  u 1:2 w linespoints  t 'GP 4'  ,\
     'forces_GP_5.dat'  u 1:2 w linespoints  t 'GP 5'  ,\
     'forces_GP_6.dat'  u 1:2 w linespoints  t 'GP 6'  ,\
     'forces_GP_7.dat'  u 1:2 w linespoints  t 'GP 7'  ,\
     'forces_GP_8.dat'  u 1:2 w linespoints  t 'GP 8'  ,\
     'forces_GP_9.dat'  u 1:2 w linespoints  t 'GP 9'  ,\
     'forces_GP_10.dat' u 1:2 w linespoints  t 'GP 10' 

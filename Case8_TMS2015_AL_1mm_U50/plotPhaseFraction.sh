#!/bin/bash

DATA1=/home/saeedipour/OpenFOAM/saeedipour-2.1.1/run/September_2014/Case8_TMS_1mm_50/log.HPDCinterFoam
DATA2=/home/saeedipour/OpenFOAM/saeedipour-2.1.1/run/September_2014/Case8_TMS_1mm_interFoam/log.interFoam

grep 'Phase-1 volume fraction =' $DATA1 | \
  cut -d' ' -f5 > alpha1.tmp

grep 'Phase-1 volume fraction =' $DATA2 | \
  cut -d' ' -f5 > alpha2.tmp

grep 'Time = ' $DATA1 | \
  grep -v Execution | cut -d' ' -f3 > time1.tmp

grep 'Time = ' $DATA2 | \
  grep -v Execution | cut -d' ' -f3 > time2.tmp

paste time1.tmp alpha1.tmp > alpha1.dat
paste time2.tmp alpha2.tmp > alpha2.dat


rm *.tmp

gnuplot << EOF

reset


set term png
set output "phaseFraction.png"

set xlabel "Time (s)"
set xtics 0.002

set ylabel "alpha (-)"

set key top right Left

plot \
    'alpha1.dat' using 1:2 with lines title "Simulation",\
    'alpha2.dat' using 1:2 with lines title "Liquid flow rate"

EOF

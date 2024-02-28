#!/usr/bin/env gnuplot

#
# file.gplot
#
# Copyright (c) 2020 Carlos Braga
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.
#

# -----------------------------------------------------------------------------
## graphics output
#set terminal aqua persist size 640,480               enhanced font 'helvetica,14'
#set terminal x11 persist size 640,480 position 0,24 enhanced font 'helvetica,14' linewidth 2.0
set terminal wxt  persist size 640,480 position 0,24 enhanced font 'helvetica,14' linewidth 2.0 dashlength 0.1
#set terminal epslatex size 8.0cm,6.0cm color solid colortext rounded font 'Times-Roman,8' standalone header \
#"\\usepackage[utf8]{inputenc}\n \
#\\usepackage[T1]{fontenc}\n \
#\\renewcommand{\\familydefault}{\\rmdefault}\n \
#\\usepackage{mathptmx}\n \
#\\usepackage{helvet}\n \
#\\newcommand{\\ft}[0]{\\footnotesize}"
# set output 'figure.tex'
reset



# -----------------------------------------------------------------------------
## plot uniform distribution
set size    1.0, 1.0
set origin  0.0, 0.0
set samples 1000

set xrange [0:1.0];  set xtics 0.5; set mxtics 5
set yrange [0:2.0];  set ytics 0.5; set mytics 5

set format x '%.1f'; set xlabel '$\mathrm{x}$'
set format y '%.1f'; set ylabel '$\mathrm{\rho(X)}$'

set key at -1.5, 0.4
set key bottom right

unset key

set term wxt 0
plot '/tmp/out.hist1d_kiss_urand' u 1:3 w l lw 1.0, \
     '/tmp/out.hist1d_cmwc_urand' u 1:3 w l lw 1.0


# -----------------------------------------------------------------------------
## plot normal distribution
set size    1.0, 1.0
set origin  0.0, 0.0
set samples 1000

set xrange [-10:10.0];  set xtics 2.0; set mxtics 5
set yrange [0:0.5];     set ytics 0.1; set mytics 5

set format x '%.1f'; set xlabel '$\mathrm{x}$'
set format y '%.1f'; set ylabel '$\mathrm{\rho(X)}$'

set key at -1.5, 0.4
set key bottom right

unset key

set term wxt 1
plot '/tmp/out.hist1d_kiss_nrand' u 1:3 w l lw 1.0, \
     '/tmp/out.hist1d_cmwc_nrand' u 1:3 w l lw 1.0

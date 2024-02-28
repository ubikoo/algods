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
# graphics output
#
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
# plot ode trajectories
#
set size    1.0, 1.0
set origin  0.0, 0.0
set samples 1000

#set xrange [0:1.0];
set xtics 0.5; set mxtics 5
#set yrange [0:2.0];
set ytics 0.5; set mytics 5


set format x '%.1e'; set xlabel '$\mathrm{x}$'
set format y '%.1e'; set ylabel '$\mathrm{\rho(X)}$'

set key at -1.5, 0.4
set key bottom right

unset key


#
# set data files prefix
#
prefix = "/tmp"


# -----------------------------------------------------------------------------
#  t  x(0)  x(1)  e1  x  v  e2
#
set term wxt 0
set title "non oscillatory spring - all"
show title
plot sprintf("%s/out.euler%d", prefix, 1) u ($2):($3) w l lw 0.5 lc rgb 'dark-red', \
     sprintf("%s/out.gauss%d", prefix, 1) u ($2):($3) w p ps 0.1 lc rgb 'dark-green', \
     sprintf("%s/out.kutta%d", prefix, 1) u ($2):($3) w p ps 0.1 lc rgb 'dark-blue', \
     sprintf("%s/out.gear%d",  prefix, 1) u ($2):($3) w p ps 0.1 lc rgb 'dark-yellow'

set term wxt 1
set title "damped oscillatory spring - all"
show title
plot sprintf("%s/out.euler%d", prefix, 2) u ($2):($3) w l lw 0.5 lc rgb 'dark-red', \
     sprintf("%s/out.gauss%d", prefix, 2) u ($2):($3) w p ps 0.1 lc rgb 'dark-green', \
     sprintf("%s/out.kutta%d", prefix, 2) u ($2):($3) w p ps 0.1 lc rgb 'dark-blue', \
     sprintf("%s/out.gear%d",  prefix, 2) u ($2):($3) w p ps 0.1 lc rgb 'dark-yellow'

set term wxt 2
set title "oscillatory spring - all"
show title
plot sprintf("%s/out.euler%d", prefix, 3) u ($2):($3) w l lw 0.5 lc rgb 'dark-red', \
     sprintf("%s/out.gauss%d", prefix, 3) u ($2):($3) w p ps 0.1 lc rgb 'dark-green', \
     sprintf("%s/out.kutta%d", prefix, 3) u ($2):($3) w p ps 0.1 lc rgb 'dark-blue', \
     sprintf("%s/out.gear%d",  prefix, 3) u ($2):($3) w p ps 0.1 lc rgb 'dark-yellow'



# -----------------------------------------------------------------------------
#  t  x(0)  x(1)  e1  x  v  e2
#
set term wxt 3
set title "non oscillatory spring - gauss, kutta, gear"
show title
plot sprintf("%s/out.gauss%d", prefix, 1) u ($2):($3) w p lw 0.5 lc rgb 'red', \
     sprintf("%s/out.kutta%d", prefix, 1) u ($2):($3) w p ps 0.1 lc rgb 'dark-green', \
     sprintf("%s/out.gear%d",  prefix, 1) u ($2):($3) w p ps 0.1 lc rgb 'dark-yellow'

set term wxt 4
set title "damped oscillatory spring - gauss, kutta, gear"
show title
plot sprintf("%s/out.gauss%d", prefix, 2) u ($2):($3) w p lw 0.5 lc rgb 'red', \
     sprintf("%s/out.kutta%d", prefix, 2) u ($2):($3) w p ps 0.1 lc rgb 'dark-green', \
     sprintf("%s/out.gear%d",  prefix, 2) u ($2):($3) w p ps 0.1 lc rgb 'dark-yellow'

set term wxt 5
set title "oscillatory spring - gauss, kutta, gear"
show title
plot sprintf("%s/out.gauss%d", prefix, 3) u ($2):($3) w p lw 0.5 lc rgb 'red', \
     sprintf("%s/out.kutta%d", prefix, 3) u ($2):($3) w p ps 0.1 lc rgb 'dark-green', \
     sprintf("%s/out.gear%d",  prefix, 3) u ($2):($3) w p ps 0.1 lc rgb 'dark-yellow'



# -----------------------------------------------------------------------------
#  t  x(0)  x(1)  e1  x  v  e2
#
set term wxt 6
set title "non oscillatory spring - gauss"
show title
plot sprintf("%s/out.gauss%d", prefix, 1) u ($2):($3) w p ps 0.1 lc rgb 'dark-green'

set term wxt 7
set title "damped oscillatory spring - gauss"
show title
plot sprintf("%s/out.gauss%d", prefix, 2) u ($2):($3) w p ps 0.1 lc rgb 'dark-green'

set term wxt 8
set title "oscillatory spring - gauss"
show title
plot sprintf("%s/out.gauss%d", prefix, 3) u ($2):($3) w p ps 0.1 lc rgb 'dark-green'

set term postscript eps enhanced color "Helvetica,22" 
set output "images/RadRes.eps"

set size 0.9,1.25

color="blue" 
set style line 11 dt 1 lc rgb color    lw 4.5 pt 1.0 ps 5
set style line 1 dt 1 lc rgb color     lw 3.0 ps 1.0 pt 5 #11

set style line 2 dt 1 lc rgb "black"   lw 2.0 pt 4  ps 1.5
set style line 3 dt 1 lc rgb "black"   lw 2.0 pt 8  ps 1.7
set style line 4 dt 1 lc rgb "black"   lw 2.0
set style line 5 dt 1 lc rgb "black"   lw 2.0 pt 12 ps 1.7
set style line 6 dt 1 lc rgb "black"   lw 2.0 pt 2  ps 1.5
set style line 7 dt 1 lc rgb "black"   lw 2.0 pt 13 ps 1.7
set style line 8 dt 1 lc rgb "black"   lw 2.0 pt 6  ps 1.5
set style line 9 dt 1 lc rgb "black"   lw 2.0 pt 10 ps 1.7
set style line 10 dt 1 lc rgb "black"  lw 2.0 pt 9  ps 1.7

set style line 20 dt 1 lw 2 lc rgb "white"

set multiplot 

set object rectangle from graph 0,0 to graph 1,1 behind fillcolor rgb "#EAEAF4" fillstyle solid noborder
set grid xtics ytics ls 20

set key font "Helvetica, 12"
set key spacing 1.2

set logscale x
set format x "10^{%L}"
set border lw 1.5
set tics scale 1.5
set xtics ("1" 1, "10^{2}" 1e2, "10^{4}" 1e4, "10^{6}" 1e6)
set xtics 100
set xlabel "" #"Time (ps)"
set xrange[1:1.5e6]

dir="/Users/jose.r/TOPAS_FORUM_NBIO/IRT/regression/"

set format x ""
set size 0.5,0.6
set origin 0.0,0.65
set rmargin 0
set bmargin 0

f=1e-1/0.1035
set yrange[2.1:5.2]
set ylabel "G value [molec. / 100 eV]" offset 1.0, -8
set key bottom left 
tcut=0 
set label "^{\267}OH" at graph 0.75,0.85
plot "<echo 7.0 4.9 0.2" u 1:($2*f):($3*f) w errorbars ls 2 notitle,\
     "<echo 10.0 4.8 0.12" u 1:($2*f):($3*f) w errorbars ls 2 notitle,\
     "<echo 296276, 2.72, 0.14" u 1:2:3 w errorbar ls 3 notitle "LaVerne 2000",\
     "<echo 493692, 2.47, 0.12" u 1:2:3 w errorbar ls 3 notitle,\
     "<echo 605751, 2.44, 0.12" u 1:2:3 w errorbar ls 3 notitle,\
     "<echo 1e6,    2.49, 0.12" u 1:2:3 w errorbar ls 3 notitle,\
     "<echo 1e6,    2.49, 0.12" u 1:2 w p ls 3 notitle "LaVerne 2000",\
     "hydroxyl_long.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 lc rgb "gray50" notitle,\
     "hydroxyl_short.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle,\
     sprintf("%sgvalue-pure-water.phsp","./") u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 11 notitle,\



set size 0.5,0.6
set origin 0.5,0.65
set lmargin 0
set rmargin 8.5
set format y ""
set xrange[2:]
set ylabel ""
tcut=1.8e1
unset label
set label "e^@{/Symbol \055}_{aq}" at graph 0.75, 0.85
plot "<echo 7.0, 4.4, 0.2" u 1:($2*f):($3*f) w errorbars ls 4 notitle,\
     "<echo 7.0, 4.4, 0.2" u 1:($2*f) w p ls 4 notitle "Wang 2018",\
     "<echo 20.0, 4.2, 0.2" u 1:($2*f):($3*f) w errorbars ls 4 notitle,\
     "<echo 70.0e3  2.93 0.2"  w errorbars ls 5 notitle "Shirashi 1988",\
     "<echo 70.0e3  2.93 0.2"  w p ls 5 notitle "Shirashi 1988",\
     "<echo 300.0e3 2.67 0.15" w errorbars ls 5 notitle,\
     "<echo 1e5,    2.7" u 1:2:($2*0.03) w errorbars pt 5 ps 1.3 lc 0 lw 4 notitle "Bartels 2000",\
     "electron_short.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle, "electron_long.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle,\
     sprintf("%sgvalue-pure-water.phsp","./") u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 11 notitle "TOPAS-nBio",\

set size 0.5,0.6
set origin 0.0,0.05
set lmargin 8.2
set rmargin 0
set bmargin 2
set tmargin 0
set format y "%g"
set format x "10^{%L}"
set xtics ("1" 1, "10^{2}" 1e2, "10^{4}" 1e4, "10^{6}" 1e6)
set xtics 100
set xlabel "Time (ps)" 

set xrange[1:]
set yrange[0.05:0.85] 
kobs=9.7e8
histories=1
unset label
set label "H_{2}O_{2}" at graph 0.05, 0.85

plot "../HydrogenPeroxide/Hiroki2020.csv" using (1e12/($1*kobs)):(1e12/($1*kobs)>0.e2 ? $2 : 1/0):($2*0.02) w errorbars ls 8 notitle "Hiroki et al., 2002",\
"../HydrogenPeroxide/Hiroki2020.csv" using (1e12/($1*kobs)):(1e12/($1*kobs)>0.e2 ? $2 : 1/0) w p ls 8 notitle "Hiroki 2002",\
     "tempH2O2Ritchie.csv" index 0 u (1e12/($1*kobs)):((1e12/($1*kobs)>0.e2 ? $2 : 1/0)):($3*sqrt(histories)) w errorbars ls 1 notitle,\
     "tempH2O2Ritchie.csv" index 0 u (1e12/($1*kobs)):(1e12/($1*kobs)>0.e2 ? $2 : 1/0) with     lp   ls 1 dt 1 notitle,\
#     sprintf("%sgvalue-pure-water.phsp","./") u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l ls 11 notitle,\


set size 0.5,0.6
set origin 0.5,0.05
set lmargin 0 #8.2
set rmargin 8.5 
set bmargin 2
set tmargin 0
set format y ""
set xrange[2:]
k1=1.1e10
k2=2.9e10
set autoscale y 
unset label
set label "H_{2}" at graph 0.05, 0.85
plot "../MolecularHydrogen/Pastina1999.csv" index 0 using (1e12/($1*k1*$2)):(1e12/($1*k1*$2) > 5e2 ? $3 : 1/0):($3*0.1) w errorbars ls 6 notitle "H2O2",\
     "../MolecularHydrogen/Pastina1999.csv" index 1 using (1e12/($1*k2*$2)):(1e12/($1*k2*$2) > 5e2 ? $3 : 1/0):($3*0.1) w errorbars ls 7 notitle "Pastina1999 K2Cr2O",\
     "../MolecularHydrogen/Pastina1999.csv" index 0 using (1e12/($1*k1*$2)):(1e12/($1*k1*$2) > 5e2 ? $3 : 1/0) w p ls 6 notitle "Pastina 1999 (H_{2}O_{2})",\
     "../MolecularHydrogen/Pastina1999.csv" index 1 using (1e12/($1*k2*$2)):(1e12/($1*k2*$2) > 5e2 ? $3 : 1/0) w p ls 7 notitle "Pastina 1999 (K_{2}Cr_{2}O_{7})",\
     sprintf("%sgvalue-pure-water.phsp","./") u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 11 notitle,\



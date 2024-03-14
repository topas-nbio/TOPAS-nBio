set term postscript eps enhanced color "Helvetica,16" 
set output "images/Figure.eps"

set multiplot layout 2,3
set ylabel offset 1.5
set lmargin 6
set rmargin 1.1

set key font "Helvetica, 11"
set key spacing 1.2

set logscale x
set format x "10^{%L}"
set xtics ("1" 1, "10^{2}" 1e2, "10^{4}" 1e4, "10^{6}" 1e6)
set xtics 100
set xlabel "" #"Time (ps)"
set xrange[1:1.5e6]

color="blue" 

set style line 1 dt 1 lc rgb color     lw 3.0 ps 1.0 pt 5 #11
set style line 2 dt 1 lc rgb "black"   lw 1.2 pt 4 ps 1.2
set style line 3 dt 1 lc rgb "black"   lw 1.2 pt 8 ps 1.5
set style line 4 dt 1 lc rgb "black"   lw 1.0
set style line 5 dt 1 lc rgb "black"   lw 1.2 pt 12 ps 1.5
set style line 6 dt 1 lc rgb "black"   lw 1.2 pt 2 ps 1.2
set style line 7 dt 1 lc rgb "black"   lw 1.2 pt 13 ps 1.5
set style line 8 dt 1 lc rgb "black"   lw 1.2 pt 6 ps 1.2
set style line 9 dt 1 lc rgb "black"   lw 1.2 pt 10 ps 1.5
set style line 10 dt 1 lc rgb "black"  lw 1.2 pt 9 ps 1.5

set style line 20 dt 1 lw 2 lc rgb "white"

set object rectangle from graph 0,0 to graph 1,1 behind fillcolor rgb "#EAEAF4" fillstyle solid noborder
set grid xtics ytics ls 20

dir="/Users/jose.r/TOPAS_FORUM_NBIO/IRT/regression/"

set ylabel "G(^{\267}OH) (molec. / 100 eV)" offset 1.5
f=1e-1/0.1035
set yrange[2:5.2]
set key bottom left #box opaque
tcut=0 

plot sprintf("%selectron_999.999_keV.gp","./") index 0 u 1:2 w l ls 1 notitle,\
     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "OH^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "black" notitle,\

     #sprintf("%selectron_999.999_keV.gp","/Users/ramosj/Projects/RadiolysisValidation/gvalueLET/e_999.999_keV/") index 0 u 1:2 w l ls 1 notitle,\

set ylabel "G(e^@{/Symbol \055}_{aq}) (molec. / 100 eV)"
set ytics 0.4
set yrange[2.4:4.5]
set key at 0.5,2.45
tcut=1.8e1

plot sprintf("%selectron_999.999_keV.gp","./") index 1 u 1:2 w l ls 1 notitle,\
     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "e_aq^-1" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "black" notitle,\

     #sprintf("%selectron_999.999_keV.gp","/Users/ramosj/Projects/RadiolysisValidation/gvalueLET/e_999.999_keV/") index 1 u 1:2 w l ls 1 notitle,\

set yrange[0.1:0.55] 
set ytics 0.1
set label "Time (ps)" at graph 0.3, graph -0.2

set ylabel "G(H_{2}) (molec. / 100 eV)"
set key at 0.5, 0.12
plot sprintf("%selectron_999.999_keV.gp","./") index 4 u 1:2 w l ls 1 notitle,\
     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "H_2^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "black" notitle,\

unset label
set tmargin 0.3
set ylabel "G(H_{2}O_{2}) (molec. / 100 eV)"
kobs=9.7e8
set format x "10^{%L}"
set xlabel "Time (ps)"
set yrange[0:0.85]
set ytics 0.2
set key at 0.5, 0.75
histories=1
plot sprintf("%selectron_999.999_keV.gp","./") index 3 u 1:2 w l ls 1 notitle,\
     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "H2O2^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "black" notitle,\

set ylabel "G(H^{\267}) (molec. / 100 eV)"
kobs=3.2e9

set autoscale xy
set xrange[1:1.5e6]
set yrange[0.2:0.85]

MC=0.422305
eMC=0.00176205
MCNoSpin=0.396928   
eMCNoSpin=0.00168968
set key at 1.0e0,0.65
plot sprintf("%selectron_999.999_keV.gp","./") index 5 u 1:2 w l ls 1 title "TsIRT",\
     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "H^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "black" title "G4Step-by-step",\




set term postscript eps enhanced color "Helvetica,16" 
set output "images/Figure_2.eps"

set multiplot layout 2,3
set ylabel offset 1.5
set lmargin 6
set rmargin 1.1

set key font "Helvetica, 12"
set key spacing 1.2

set logscale x
set format x "10^{%L}"
set xtics ("1" 1, "10^{2}" 1e2, "10^{4}" 1e4, "10^{6}" 1e6)
set xtics 100
set xlabel "" #"Time (ps)"
set xrange[1:1.5e6]

color="blue" 

set style line 11 dt 1 lc rgb color    lw 3.0 pt 1.0 ps 5
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
set key bottom left 
tcut=0 
set label "A)" at graph -0.3, graph 1.02

plot "<echo 7.0 4.9 0.2" u 1:($2*f):($3*f) w errorbars ls 2 notitle,\
     "<echo 10.0 4.8 0.12" u 1:($2*f):($3*f) w errorbars ls 2 notitle,\
     "<echo 296276, 2.72, 0.14" u 1:2:3 w errorbar ls 3 notitle "LaVerne 2000",\
     "<echo 493692, 2.47, 0.12" u 1:2:3 w errorbar ls 3 notitle,\
     "<echo 605751, 2.44, 0.12" u 1:2:3 w errorbar ls 3 notitle,\
     "<echo 1e6,    2.49, 0.12" u 1:2:3 w errorbar ls 3 notitle,\
     "<echo 1e6,    2.49, 0.12" u 1:2 w p ls 3 notitle "LaVerne 2000",\
     "hydroxyl_long.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 lc rgb "gray50" notitle,\
     "hydroxyl_short.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle,\
     sprintf("%sgvalue-pure-water-temp.phsp","./") u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 11 notitle,\

     #"<echo 1e6,    2.53" u 1:2 w p pt 5 ps 1 lc 0 notitle "Elliot 2009",\
#     sprintf("%sgvalue-pure-water-G4.phsp","./") u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 1 notitle,\
#     sprintf("%sgvalue-pure-water-opt-old.phsp","./") u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 11 notitle,\
#     sprintf("%sgvalue-pure-water-opt-noSpin.phsp","./") u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 1 notitle,\
#     sprintf("%sgvalue-pure-water-opt-noSpin-g4dev.phsp","./") u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\
#     sprintf("%selectron_Gvalue_SBS_Fix_0.01ps.phsp","./") u 3:(strcol(4) eq "OH^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\
     #sprintf("%sgvalue-pure-water-opt.phsp","./") u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 1 notitle,\
#     "/Users/ramosj/Projects/RadiolysisValidation/gvalueLET/e_999.999_keV/electron_999.999_keV.gp" index 0 u 1:2 w l ls 1 dt 4 notitle,\
#     "LaVerneModel.csv" index 1 u 1:2 w l dt 4 lw 4 lc rgb "green" notitle,\

unset label
set label "B)" at graph -0.3, graph 1.02
set ylabel "G(e^@{/Symbol \055}_{aq}) (molec. / 100 eV)"
set ytics 0.4
set yrange[2.4:4.5]
#set key at 0.5,2.45
tcut=1.8e1

plot "<echo 7.0, 4.4, 0.2" u 1:($2*f):($3*f) w errorbars ls 4 notitle,\
     "<echo 7.0, 4.4, 0.2" u 1:($2*f) w p ls 4 notitle "Wang 2018",\
     "<echo 20.0, 4.2, 0.2" u 1:($2*f):($3*f) w errorbars ls 4 notitle,\
     "<echo 70.0e3  2.93 0.2"  w errorbars ls 5 notitle "Shirashi 1988",\
     "<echo 70.0e3  2.93 0.2"  w p ls 5 notitle "Shirashi 1988",\
     "<echo 300.0e3 2.67 0.15" w errorbars ls 5 notitle,\
     "<echo 1e5,    2.7" u 1:2:($2*0.03) w errorbars pt 5 ps 1.3 lc 0 notitle "Bartels 2000",\
     "electron_short.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle, "electron_long.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle,\
     sprintf("%sgvalue-pure-water-temp.phsp","./") u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 11 notitle "TOPAS-nBio",\

#     sprintf("%sgvalue-pure-water-G4.phsp","./") u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 1 notitle,\
#
#     "<echo 1e6,    2.64" u 1:2 w p pt 5 ps 1 lc 0 notitle "Elliot 2009",\

#     sprintf("%sgvalue-pure-water-opt-old.phsp","./") u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 11 title "G4DNA-IRT",\
#     sprintf("%sgvalue-pure-water-opt-noSpin.phsp","./") u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 1 title "This work",\
#     sprintf("%sgvalue-pure-water-opt-noSpin-g4dev.phsp","./") u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 1 lc rgb "green" title "This work-g4dev",\

#     sprintf("%selectron_Gvalue_SBS_Fix_0.01ps.phsp","./") u 3:(strcol(4) eq "e_aq^-1" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\
     #sprintf("%sgvalue-pure-water-opt.phsp","./") u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 1 notitle,\
#     "/Users/ramosj/Projects/RadiolysisValidation/gvalueLET/e_999.999_keV/electron_999.999_keV.gp" index 1 u 1:2 w l ls 1 dt 4 notitle,\

     #"Fanning_1975.csv" u 1:2 w p pt 9 ps 1 lc 0 notitle "Fanning 1975",\
#     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "e_aq^-1" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "red" notitle,\

k1=1.1e10
k2=2.9e10
set yrange[0.1:0.55] 
set ytics 0.1

unset label
set label "C)" at graph -0.3, graph 1.02
set label "Time (ps)" at graph 0.3, graph -0.2

set ylabel "G(H_{2}) (molec. / 100 eV)"
set key at 0.5, 0.12
plot "../MolecularHydrogen/Pastina1999.csv" index 0 using (1e12/($1*k1*$2)):(1e12/($1*k1*$2) > 5e2 ? $3 : 1/0):($3*0.1) w errorbars ls 6 notitle "H2O2",\
     "../MolecularHydrogen/Pastina1999.csv" index 1 using (1e12/($1*k2*$2)):(1e12/($1*k2*$2) > 5e2 ? $3 : 1/0):($3*0.1) w errorbars ls 7 notitle "Pastina1999 K2Cr2O",\
     "../MolecularHydrogen/Pastina1999.csv" index 0 using (1e12/($1*k1*$2)):(1e12/($1*k1*$2) > 5e2 ? $3 : 1/0) w p ls 6 notitle "Pastina 1999 (H_{2}O_{2})",\
     "../MolecularHydrogen/Pastina1999.csv" index 1 using (1e12/($1*k2*$2)):(1e12/($1*k2*$2) > 5e2 ? $3 : 1/0) w p ls 7 notitle "Pastina 1999 (K_{2}Cr_{2}O_{7})",\
     sprintf("%sgvalue-pure-water-temp.phsp","./") u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 11 notitle,\

#plot "<echo 1e6,    0.42" u 1:2 w p pt 5 ps 1 lc 0 notitle "Elliot 2009",\
#     sprintf("%sgvalue-pure-water-G4.phsp","./") u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 1 notitle,\

#     sprintf("%sgvalue-pure-water-opt-old.phsp","./") u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 11 notitle,\
#     sprintf("%sgvalue-pure-water-opt-noSpin.phsp","./") u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 1 notitle,\
#     sprintf("%sgvalue-pure-water-opt-noSpin-g4dev.phsp","./") u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\

#     sprintf("%selectron_Gvalue_SBS_Fix_0.01ps.phsp","./") u 3:(strcol(4) eq "H_2^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\

     #sprintf("%sgvalue-pure-water-opt.phsp","./") u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 1 notitle,\
#     "/Users/ramosj/Projects/RadiolysisValidation/gvalueLET/e_999.999_keV/electron_999.999_keV.gp" index 4 u 1:2 w l ls 1 dt 4 notitle,\
#     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "H_2^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "red" notitle,\
    "H2_vs_H2O2_scav.csv" u 1:2:3 w errorbars ls 1 notitle,\
    "H2_vs_H2O2_scav_var.csv" u 1:2:3 w errorbars ls 1 lc rgb "red" notitle,\

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
unset label
set label "D)" at graph -0.3, graph 1.02
plot "../HydrogenPeroxide/Hiroki2020.csv" using (1e12/($1*kobs)):(1e12/($1*kobs)>0.e2 ? $2 : 1/0):($2*0.02) w errorbars ls 8 notitle "Hiroki et al., 2002",\
     "../HydrogenPeroxide/Hiroki2020.csv" using (1e12/($1*kobs)):(1e12/($1*kobs)>0.e2 ? $2 : 1/0) w p ls 8 notitle "Hiroki 2002",\
     "tempH2O2Ritchie.csv" index 0 u (1e12/($1*kobs)):((1e12/($1*kobs)>0.e2 ? $2 : 1/0)):($3*sqrt(histories)) w errorbars ls 1 notitle,\
     "tempH2O2Ritchie.csv" index 0 u (1e12/($1*kobs)):(1e12/($1*kobs)>0.e2 ? $2 : 1/0) with     lp   ls 1 dt 3 notitle,\
     sprintf("%sgvalue-pure-water-temp.phsp","./") u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l ls 11 notitle,\

     #"<echo 1e6,    0.75" u 1:2 w p pt 5 ps 1 lc 0 notitle "Elliot 2008",\
#     sprintf("%sgvalue-pure-water-G4.phsp","./") u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l ls 11 notitle,\
#     "../HydrogenPeroxide/Pastina1999.csv" using (1e12/$1):2:($2*0.05) w errorbars ls 9 notitle "Pastina 1999",\
#     "../HydrogenPeroxide/Pastina1999.csv" using (1e12/$1):2           w p ls 9 notitle "Pastina 1999",\
#     sprintf("%sgvalue-pure-water-Meesugnoen.phsp","./") u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l ls 1 notitle,\
#     sprintf("%sgvalue-pure-water-opt-noSpin.phsp","./") u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l ls 1 notitle,\
#     sprintf("%sgvalue-pure-water-opt-noSpin-g4dev.phsp","./") u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\
#     "tempH2O2Ritchie.csv" u (1e12/($1*kobs)):((1e12/($1*kobs)>0.e2 ? $2 : 1/0)):($3*sqrt(histories)) w errorbars ls 1 notitle,\
#     "tempH2O2Ritchie.csv" u (1e12/($1*kobs)):(1e12/($1*kobs)>0.e2 ? $2 : 1/0) with     lp   ls 1 dt 3 notitle,\
#     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "H2O2^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "red" notitle,\

     #sprintf("%sgvalue-pure-water-opt.phsp","./") u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l ls 1 notitle,\
     #sprintf("%selectron_Gvalue_SBS_Fix_0.01ps.phsp","./") u 3:(strcol(4) eq "H2O2^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\
     #"/Users/ramosj/Projects/RadiolysisValidation/gvalueLET/e_999.999_keV/electron_999.999_keV.gp" index 3 u 1:2 w l ls 1 dt 4 notitle,\

set ylabel "G(H^{\267}) (molec. / 100 eV)"
#set label "1 mM" at 0.4e2, 0.8
#set label "24 mM" at 0.5e3, 0.3
kobs=3.2e9
#set label "Huerta-Parajon 2008" at 1.5, 0.1 font "Helvetica, 11"

set autoscale xy
set xrange[1:1.5e6]
set yrange[0.2:0.85]
#set yrange[0:0.99]
unset label
set label "E)" at graph -0.3, graph 1.02

MC=0.422305
eMC=0.00176205
MCNoSpin=0.396928   
eMCNoSpin=0.00168968
set key at 1.0e0,0.65

plot sprintf("%sgvalue-pure-water-temp.phsp","./") u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 11 notitle "No-spin",\
     "../HydrogenAtom/HuertaParajon.csv" index 0 u (1e12/($1*2.1e8 + 1e-3 * 1.4e6)):($2-4.60526e-1):($2*0.05) w errorbars ls 10 notitle,\
     "H2_minus_X_equal_to_H_noSpin.csv" index 1 u 1:($2-MCNoSpin):(sqrt($3*$3+eMCNoSpin*eMCNoSpin)) w errorbars ls 1 notitle,\
     "H2_minus_X_equal_to_H_noSpin.csv" index 1 u 1:($2-MCNoSpin) w lp ls 1 dt 3 notitle,\

     #"<echo 1e6,    0.56" u 1:2 w p pt 5 ps 1 lc rgb "black" notitle,\
#plot "<echo 1e6,    0.56" u 1:2 w p pt 5 ps 1 lc rgb "black" notitle,\
#     sprintf("%sgvalue-pure-water-G4.phsp","./") u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 1 notitle "No-spin",\
     #sprintf("%sgvalue-pure-water-opt-noSpin.phsp","./") u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 1 notitle "No-spin",\
     #sprintf("%sgvalue-pure-water-opt-noSpin-g4dev.phsp","./") u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 1 lc rgb "green" notitle "No-spin",\
#     sprintf("%selectron_Gvalue_SBS_Fix_0.01ps.phsp","./") u 3:(strcol(4) eq "H^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "green" notitle,\

     #"H2_minus_X_equal_to_H.csv" index 1 u 1:($2-MC):(sqrt($3*$3+eMC*eMC)) w errorbars ls 1 notitle,\
     #"H2_minus_X_equal_to_H.csv" index 1 u 1:($2-MC) w lp ls 1 dt 3 notitle,\
     #sprintf("%sgvalue-pure-water-opt.phsp","./") u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 1 notitle "Spin",\

     #"/Users/ramosj/Projects/RadiolysisValidation/gvalueLET/e_999.999_keV/electron_999.999_keV.gp" index 5 u 1:2 w l ls 1 dt 4 notitle,\
     #"H2_minus_X_equal_to_H.csv" index 0 u 1:($2-MC):(sqrt($3*$3+eMC*eMC)) w errorbars ls 1 notitle,\
     #"H2_minus_X_equal_to_H.csv" index 0 u 1:($2-MC) w lp ls 1 dt 3 notitle,\

     #"../HydrogenAtom/HuertaParajon.csv" index 1 u (1e12/($1)):($2-4.60526e-1):($2*0.05) w errorbars ls 10 pt 28 notitle,\
#     sprintf("%selectron_Gvalue_SBS_Fix_0.01.phsp","./") u 3:(strcol(4) eq "H^0" && $3 < 1e6? $1 : 1/0) w l ls 1 lc rgb "red" notitle,\





set term postscript eps enhanced color "Helvetica,16" 
set output "Figure_2b1.eps"

set multiplot layout 2,1
set ylabel offset 1.5
#set lmargin 7
#set rmargin 1.1

set key font "Helvetica, 11"
set key spacing 1.2

set logscale x
set format x "10^{%L}"
set xtics 100

set xlabel "Time (ps)"
#set xrange[1:1.5e6]

color="royalblue"
set style line 1 dt 1 lc rgb "black"   lw 2.2 ps 0.5 pt 4
set style line 2 dt 2 lc rgb "black"   lw 2.2 pt 4 ps 0.8
set style line 3 dt 3 lc rgb "black"   lw 2.2 pt 6 ps 0.8
set style line 4 dt 4 lc rgb "black"   lw 2.2 pt 4 ps 0.8
set style line 5 dt 5 lc rgb "black"   lw 2.2 pt 10 ps 1.2
set style line 6 dt 6 lc rgb "black"   lw 2.2 pt 12 ps 1.2
set style line 7 dt 7 lc rgb "black"   lw 2.2 pt 14 ps 1.2
set style line 8 dt 8 lc rgb "black"   lw 2.2 pt 7 ps 1.2
set style line 9 dt 9 lc rgb "black"   lw 2.2 pt 3 ps 1.2
set style line 10 dt 1 lc rgb "black"  lw 2.2 pt 22 ps 1.2

dir="/Users/jose.r/TOPAS_FORUM_NBIO/IRT/regression/"

set ylabel "G(^{\267}OH) (molec. / 100 eV)"
set ytics nomirror 
set y2tics 0.2
set key bottom left
#set yrange[-0.5:5]
#set y2range[-0.1:0.85]

plot "temporal/HCO2m_1e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 1 title "Xe-2",\
     "temporal/HCO2m_1e-2_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 2 title "Xe-1",\
     "temporal/HCO2m_300e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 3 title "300e-3",\
     "temporal/HCO2m_1e-1_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 4 title "Xe-1",\
     "temporal/HCO2m_1e-0_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls 4 title "Xe+0",\
     "temporal/HCO2m_1e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 1 lc 2 axis x1y2 notitle "Xe-2",\
     "temporal/HCO2m_1e-2_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 2 lc 2 axis x1y2 notitle "Xe-1",\
     "temporal/HCO2m_300e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 3 lc 2 axis x1y2 notitle "Xe-0",\
     "temporal/HCO2m_1e-1_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 4 lc 2 axis x1y2 notitle "Xe+1",\
     "temporal/HCO2m_1e-0_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls 5 lc 2 axis x1y2 notitle "Xe+1",\

set ylabel "G(OH) (molec. / 100 eV)"
set y2label "G(e-) (molec. / 100 eV)"
#kobs=9.7e8

#set yrange[0:1] 
set ytics 1
set y2tics 1
#set key at 0.5, 0.75

plot "temporal/HCO2m_1e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 1 title "Xe-2",\
     "temporal/HCO2m_1e-2_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 2 title "Xe-1",\
     "temporal/HCO2m_300e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 3 title "300e-3",\
     "temporal/HCO2m_1e-1_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 4 title "Xe-1",\
     "temporal/HCO2m_1e-0_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls 4 title "Xe+0",\
     "temporal/HCO2m_1e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 1 lc 2 axis x1y2 notitle "Xe-2",\
     "temporal/HCO2m_1e-2_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 2 lc 2 axis x1y2 notitle "Xe-1",\
     "temporal/HCO2m_300e-3_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 3 lc 2 axis x1y2 notitle "Xe-0",\
     "temporal/HCO2m_1e-1_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 4 lc 2 axis x1y2 notitle "Xe+1",\
     "temporal/HCO2m_1e-0_M_Br_1e-3_M_NO3m_1e-3_M.phsp" u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls 5 lc 2 axis x1y2 notitle "Xe+1",\



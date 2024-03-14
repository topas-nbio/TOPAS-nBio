set term postscript enhanced color "Helvetica"
set output "images/deltaG.eps"
set multiplot layout 2,3

set logscale yx
set format xy "10^{%L}"

set style line 1 lw 2 lc 0
set style line 2 lw 2 lc 1
set style line 3 lw 2 lc 2
set style line 4 lw 2 lc 3
set style line 5 lw 2 lc 4
set style line 6 lw 2 lc 5

set yrange[1e-6:1]

plot "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==1 ? $8 : 1/0) w l ls 1 title "1",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==2 ? $8 : 1/0) w l ls 2 title "2",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==3 ? $8 : 1/0) w l ls 3 title "3",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==4 ? $8 : 1/0) w l ls 4 title "4",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==5 ? $8 : 1/0) w l ls 5 title "5",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==6 ? $8 : 1/0) w l ls 6 title "6",\

plot "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==7 ? $8 : 1/0) w l ls 1 title "7",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==8 ? $8 : 1/0) w l ls 2 title "8",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==9 ? $8 : 1/0) w l ls 3 title "9",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==10 ? $8 : 1/0) w l ls 4 title "10",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==11 ? $8 : 1/0) w l ls 5 title "11",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==12 ? $8 : 1/0) w l ls 6 title "12",\

plot "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==13 ? $8 : 1/0) w l ls 1 title "13",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==14 ? $8 : 1/0) w l ls 2 title "14",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==15 ? $8 : 1/0) w l ls 3 title "15",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==16 ? $8 : 1/0) w l ls 4 title "16",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==17 ? $8 : 1/0) w l ls 5 title "17",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==18 ? $8 : 1/0) w l ls 6 title "18",\

plot "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==19 ? $8 : 1/0) w l ls 1 title "19",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==20 ? $8 : 1/0) w l ls 2 title "20",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==21 ? $8 : 1/0) w l ls 3 title "21",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==22 ? $8 : 1/0) w l ls 4 title "22",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==23 ? $8 : 1/0) w l ls 5 title "23",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==24 ? $8 : 1/0) w l ls 6 title "24",\

plot "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==25 ? $8 : 1/0) w l ls 1 title "25",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==26 ? $8 : 1/0) w l ls 2 title "26",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==27 ? $8 : 1/0) w l ls 3 title "27",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==28 ? $8 : 1/0) w l ls 4 title "28",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==29 ? $8 : 1/0) w l ls 5 title "29",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==30 ? $8 : 1/0) w l ls 6 title "30",\
     "gvalue-pure-water-opt_DeltaG.phsp" u 7:($1==31 ? $8 : 1/0) w l ls 1 title "31",\

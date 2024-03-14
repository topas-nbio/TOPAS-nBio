import os
def plot(fileNames, titles):
    out = open('Figure_2_Automatic.gp','w')
    out.write('set terminal pngcairo enhanced size 2800, 1950 font "Helvetica, 44"\n')
    out.write('set output "images/Figure_2_Automatic.png"\n')
    out.write('\n')
    out.write('set ylabel offset 1.5\n')
    out.write('set lmargin 7\n')
    out.write('set rmargin 1.1\n')
    out.write('\n')
    out.write('set key font "Helvetica, 34"\n')
    out.write('set key spacing 1.5\n')
    out.write('set key top right\n')
    out.write('\n')
    out.write('set logscale x\n')
    out.write('set format x "10^{%L}"\n')
    out.write('set xtics 100\n')
    out.write('\n')
    out.write('set xlabel "Time (ps)"\n')
    out.write('set xrange[1:2.5e6]\n')
    out.write('\n')
    out.write('color="black"\n')
    out.write('set style line 1 dt 1 lc rgb color     lw 4.0 ps 0.5 pt 4\n')
    out.write('set style line 2 dt 1 lc rgb "blue"   lw 4.0 pt 4 ps 1.\n')
    out.write('set style line 3 dt 1 lc rgb "red"   lw 4.0 pt 8 ps 1.\n')
    out.write('set style line 4 dt 1 lc rgb "green"   lw 4.0\n')
    out.write('set style line 5 dt 1 lc 0 lw 4.0 pt 12 ps 1.\n')
    out.write('set style line 6 dt 2 lc 1   lw 4.0 pt 2 ps 1.\n')
    out.write('set style line 7 dt 3 lc 2   lw 4.0 pt 14 ps 1.\n')
    out.write('set style line 8 dt 4 lc 3   lw 4.0 pt 6 ps 1.\n')
    out.write('set style line 9 dt 5 lc 4   lw 4.0 pt 10 ps 1.\n')
    out.write('set style line 10 dt 1 lc rgb "black"  lw 4.0 pt 22 ps 1.\n')

    out.write('set style line 20 dt 1 lw 2 lc rgb "white"\n')
    out.write('set object rectangle from graph 0,0 to graph 1,1 behind fillcolor rgb "#EAEAF4" fillstyle solid noborder\n')
    out.write('set grid xtics ytics ls 20\n')

    out.write('\n')
    out.write('set grid\n')
    out.write('\n')
    out.write('set multiplot layout 2,3\n')
    out.write('set ylabel "G(^{\267}OH) (molec. / 100 eV)"\n')
    out.write('f=1e-1/0.1035\n')
    out.write('set yrange[2:5.2]\n')
    out.write('tcut=0\n')
    out.write('plot "<echo 7.0 4.9 0.2" u 1:($2*f):($3*f) w errorbars ls 2 notitle,\\\n')
    out.write('     "<echo 296276, 2.72, 0.14" u 1:2:3 w errorbar ls 3 notitle "LaVerne 2000",\\\n')
    out.write('     "<echo 493692, 2.47, 0.12" u 1:2:3 w errorbar ls 3 notitle,\\\n')
    out.write('     "<echo 605751, 2.44, 0.12" u 1:2:3 w errorbar ls 3 notitle,\\\n')
    out.write('     "<echo 1e6,    2.49, 0.12" u 1:2:3 w errorbar ls 3 notitle,\\\n')
    out.write('     "<echo 1e6,    2.49, 0.12" u 1:2 w p ls 3 title "LaVerne 2000",\\\n')
    out.write('     "<echo 1e6,    2.53" u 1:2 w p pt 5 ps 1 lc 0 title "Elliot 2009",\\\n')
    out.write('     "hydroxyl_long.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 lc rgb "gray50" notitle,\\\n')
    out.write('     "hydroxyl_short.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle,\\\n')
    out.write('     "LaVerneModel.csv" index 1 u 1:2 w l dt 4 lw 4 lc rgb "gray60" notitle,\\\n')
#    out.write('     "gvalue-pure-water-opt.phsp" u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l dt 1 lw 4 lc rgb "blue" notitle,\\\n')
    marker = 1
    for name, title in zip(fileNames, titles):
        out.write('     "%s" u 3:(strcol(4) eq "OH^0" ? $1 : 1/0) w l ls %d notitle "%s",\\\n' % (name, marker, title))
        marker += 1
    out.write('    \n')
    out.write('\n\n')
    out.write('set autoscale y\n')
    out.write('tcut=1.8e1\n')
    out.write('plot ')
    out.write('     "<echo 7.0, 4.4, 0.2" u 1:($2*f):($3*f) w errorbars ls 4 notitle,\\\n')
    out.write('     "<echo 7.0, 4.4, 0.2" u 1:($2*f) w p ls 4 notitle "Wang 2018",\\\n')
    out.write('     "<echo 20.0, 4.2, 0.2" u 1:($2*f):($3*f) w errorbars ls 4 notitle,\\\n')
    out.write('     "<echo 70.0e3  2.93 0.2"  w errorbars ls 5 notitle "Shirashi 1988",\\\n')
    out.write('     "<echo 70.0e3  2.93 0.2"  w p ls 5 notitle "Shirashi 1988",\\\n')
    out.write('     "<echo 300.0e3 2.67 0.15" w errorbars ls 5 notitle,\\\n')
    out.write('     "<echo 1e5,    2.7" u 1:2:($2*0.03) w errorbars pt 5 ps 1 lc 1 notitle "Elliot 2009",\\\n')
    out.write('     "<echo 1e6,    2.64" u 1:2 w p pt 5 ps 1 lc 0 notitle "Elliot 2009",\\\n')
    out.write('     "electron_short.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle, "electron_long.csv" u 1:($1>tcut ? $2 : 1/0) w l ls 4 notitle,\\\n')
#    out.write('     "gvalue-pure-water-opt.phsp" u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l dt 1 lw 4 lc rgb "blue" notitle,\\\n')
    marker = 1
    for name, title in zip(fileNames, titles):
        out.write('     "%s" u 3:(strcol(4) eq "e_aq^-1" ? $1 : 1/0) w l ls %d title "%s",\\\n' % (name, marker, title))
        marker += 1
    out.write('\n\n')

    out.write('plot ')
#    out.write('     "gvalue-pure-water-opt.phsp" u 3:(strcol(4) eq "H3O^1" ? $1 : 1/0) w l dt 1 lw 4 lc rgb "blue" notitle,\\\n')
    marker = 1
    for name, title in zip(fileNames, titles):
        out.write('     "%s" u 3:(strcol(4) eq "H3O^1" ? $1 : 1/0) ls %d notitle "%s",\\\n' % (name, marker, title))
        marker += 1
    out.write('\n\n')

    out.write('plot ')
#    out.write('     "gvalue-pure-water-opt.phsp" u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) w l dt 1 lw 4 lc rgb "blue" notitle,\\\n')
    marker = 1
    for name, title in zip(fileNames, titles):
        out.write('     "%s" u 3:(strcol(4) eq "H2O2^0" ? $1 : 1/0) ls %d notitle "%s",\\\n' % (name, marker, title))
        marker += 1
    out.write('\n\n')
    out.write('k1=1.1e10\n')
    out.write('k2=2.9e10\n')

    out.write('plot ')
    out.write('     "<echo 1e6,    0.42" u 1:2 w p pt 5 ps 1 lc 0 notitle "Elliot 2009",\\\n')
    out.write('     "../MolecularHydrogen/Pastina1999.csv" index 0 using (1e12/($1*k1*$2)):(1e12/($1*k1*$2) > 5e2 ? $3 : 1/0):($3*0.1) w errorbars ls 6 lc 0 notitle "H2O2",\\\n')
    out.write('     "../MolecularHydrogen/Pastina1999.csv" index 1 using (1e12/($1*k2*$2)):(1e12/($1*k2*$2) > 5e2 ? $3 : 1/0):($3*0.1) w errorbars ls 7 lc 0 notitle "Pastina1999 K2Cr2O",\\\n')
    out.write('     "../MolecularHydrogen/Pastina1999.csv" index 0 using (1e12/($1*k1*$2)):(1e12/($1*k1*$2) > 5e2 ? $3 : 1/0) w p ls 6 lc 0 notitle "Pastina 1999 (H_{2}O_{2})",\\\n')
    out.write('     "../MolecularHydrogen/Pastina1999.csv" index 1 using (1e12/($1*k2*$2)):(1e12/($1*k2*$2) > 5e2 ? $3 : 1/0) w p ls 7 lc 0 notitle "Pastina 1999 (K_{2}Cr_{2}O_{7})",\\\n')
#    out.write('     "gvalue-pure-water-opt.phsp" u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l dt 1 lw 4 lc rgb "blue" notitle,\\\n')
    marker = 1
    for name, title in zip(fileNames, titles):
        out.write('     "%s" u 3:(strcol(4) eq "H_2^0" ? $1 : 1/0) w l ls %d notitle "%s",\\\n' % (name, marker, title))
        marker += 1
    out.write('\n\n')

    out.write('plot ')
#    out.write('     "gvalue-pure-water-opt.phsp" u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l dt 1 lw 4 lc rgb "blue" notitle,\\\n')
    marker = 1
    for name, title in zip(fileNames, titles):
        out.write('     "%s" u 3:(strcol(4) eq "H^0" ? $1 : 1/0) w l ls %d notitle "%s",\\\n' % (name, marker, title))
        marker += 1
    out.write('\n\n')

plot([\
       'gvalue-pure-water-opt-noSpin.phsp',\
       'gvalue-pure-water-Meesugnoen.phsp',\
       'Reference.phsp',\
       'gvalue-pure-water-wholeEtrack.phsp',\
#Geant4-Dev.phsp',\
      ],\
      [\
      'Old',\
      'Mees',\
      'Ref',\
      'long',\
      ])

os.system('gnuplot Figure_2_Automatic.gp')



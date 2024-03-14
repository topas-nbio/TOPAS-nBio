set ylabel "G(H^{\267}) (molec. / 100 eV)"
set label "1 mM" at 0.5e3, 0.9
set label "24 mM" at 0.5e3, 0.3
kobs=3.2e9
set label "Huerta-Parajon 2008" at 1.5, 0.1 font "Helvetica, 11"

set logscale x
plot "<echo 1e6,    0.56" u 1:2 w p pt 5 ps 1 lc rgb "black" notitle,\
     "HuertaParajon.csv" index 0 u (1e12/($1*2.1e8 + 1e-3 * 1.4e6)):($2-4.60526e-1):($2*0.05) w errorbars ls 10 notitle,\
     "HuertaParajon.csv" index 1 u (1e12/($1)):($2-4.60526e-1):($2*0.05) w errorbars ls 10 pt 28 notitle,\
     "H2_minus_X_equal_to_H.csv" index 0 u 1:($2-0.398095):3 w errorbars ls 2 notitle,\
     "H2_minus_X_equal_to_H.csv" index 1 u 1:($2-0.398095):3 w errorbars ls 1 notitle,\
     "H2_minus_X_equal_to_H.csv" index 0 u 1:($2-0.398095) w lp ls 2 notitle,\
     "H2_minus_X_equal_to_H.csv" index 1 u 1:($2-0.398095) w lp ls 1 notitle,\

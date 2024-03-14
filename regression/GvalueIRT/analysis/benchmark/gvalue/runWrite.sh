gtopas gvalueIRTElectron.txt 
gnuplot Figure_2.gp
convert -density 300 images/Figure_2.{eps,png}

################## Hydrogen peroxide
concentration=(e-2 e-1 e-0 e+1 e+2 e+3 e+4)
for i in 0 1 2 3 4
 do 
   c=${concentration[$i]}
   sed 's/theConcentration/0.1031'$c'/g' HydrogenPeroxide.txt > runme.txt
   gtopas runme.txt > /dev/null
   sed 's/theConcentration/0.1031'$c'/g' HydrogenPeroxide0.25e-3M.txt  > runme.txt
   gtopas runme.txt > /dev/null
   echo "################################ " $i
done

concentration=(e-2 e-1 e-0 e+1 e+2 e+3 e+4)
for i in 0 1 2 3 4 
  do 
    c=${concentration[$i]}
    echo 0.1031$c `sed -n '200,200p' temporal/Methanol_0.1031"$c"_M_Nitrate_25e-3_M.phsp` | awk '{print $1, $2, $3, $3*sqrt(2605), $4, $5}' 
done > tempH2O2Ritchie.csv
echo "" >> tempH2O2Ritchie.csv
echo "" >> tempH2O2Ritchie.csv
for i in 0 1 2 3 4 
  do 
    c=${concentration[$i]}
    echo 0.1031$c `sed -n '200,200p' temporal/Methanol_0.1031"$c"_M_Nitrate_0.25e-3_M.phsp` | awk '{print $1, $2, $3, $3*sqrt(2605), $4, $5}' 
done >> tempH2O2Ritchie.csv

gnuplot Figure_2.gp
#gnuplot Figure_2a.gp
convert -density 300 images/Figure_2.{eps,png}

################### Hydrogen atom
concentration=(1e-2 1e-1 300e-3 1e-0)
for i in 0 1 2 3 
 do
   c=${concentration[$i]}
   sed 's/theConcentration/'$c'/g' HydrogenAtom.txt > runme.txt 
   gtopas runme.txt  #> /dev/null
   sed 's/theConcentration/'$c'/g' HydrogenAtom24mM.txt > runme.txt 
   gtopas runme.txt  #> /dev/null
   echo $i
done 

kobsBr=1.1e10 # c=1e-3
kobsHCO2m=2.1e8 #3.2e9 # c=variable  Formate + H -> H2 + CO2-
echo "# 24 mM" > H2_minus_X_equal_to_H.csv 
for i in 0 1 2 3 
  do
    c=${concentration[$i]}
    echo $c `sed -n '500,500p' temporal/HCO2m_"$c"_M_Br_1e-3_M_NO3m_24e-3_M.phsp` | awk '{print 1e12/($1*2.1e8+1e-3*1.4e6), $2, $3, $4, $5}'
done >> H2_minus_X_equal_to_H.csv
echo "" >> H2_minus_X_equal_to_H.csv
echo "" >> H2_minus_X_equal_to_H.csv
echo "# 1 mM" >> H2_minus_X_equal_to_H.csv
for i in 0 1 2 3
  do
    c=${concentration[$i]}
    echo $c `sed -n '500,500p' temporal/HCO2m_"$c"_M_Br_1e-3_M_NO3m_1e-3_M.phsp` | awk '{print 1e12/($1*2.1e8+1e-3*1.4e6), $2, $3, $4, $5}'
done  >> H2_minus_X_equal_to_H.csv
gtopas HydrogenAtom_H20_1mM_Br_1mM_NO3m.txt
gnuplot Figure_2.gp
convert images/Figure_2.{eps,png}


################## Molecular hydrogen
#concentration=(0.000356 0.009 0.1435 1.0 4.1 9.225)
#factor=(1.004 1.018 1.070 1.175 1.307 1.396) 
##factor=(1.0   1.0   1.0   1.0   1.0   1.0)
#for i in `seq 0 5`
# do
#   c=${concentration[$i]}
#   f=${factor[$i]}
#   sed 's/theConcentration/'$c'/g' MolecularHydrogen.txt | sed 's/theReactionRate/'$f'/g' > runme.txt
#   gtopas runme.txt 
#   echo "####################" $c " and " $f
#done
#concentration=(0.000356 0.009 0.1435 1.0 4.1 9.225)
#factor=(1.004 1.018 1.070 1.175 1.307 1.396) 
##factor=(1.0   1.0   1.0   1.0   1.0   1.0)
#echo "#Scav" > H2_vs_H2O2_scav_var.csv 
#for i in `seq 0 5`
# do
#   c=${concentration[$i]}
#   f=${factor[$i]}
#   echo $c `sed -n '300,300p' temporal/H2O2_"$c"_M_KBr_1e-3_M_Factor_"$f".phsp` $f | awk '{print 1e12/($1*1.1e10*$6), $2, $3, $4, $5}' 
#done >> H2_vs_H2O2_scav_var.csv
#


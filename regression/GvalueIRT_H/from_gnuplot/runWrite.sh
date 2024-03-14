################### Hydrogen atom
concentration=(1e-2 1e-1 300e-3 1e-0)

kobsBr=1.1e10 # c=1e-3
kobsHCO2m=2.1e8 #3.2e9 # c=variable  Formate + H -> H2 + CO2-
echo "# 24 mM" > H2_minus_X_equal_to_H.csv 
for i in 0 1 2 3 
  do
    c=${concentration[$i]}
    echo $c `sed -n '500,500p' /Applications/TOPAS/openTOPAS/openTOPAS_nBio_dev/TOPAS-nBio-dev/regression/GvalueIRT_H/run/2024Mar8/mainPython/0/HCO2m_"$c"_M_Br_1e-3_M_NO3m_24e-3_M.phsp` | awk '{print 1e12/($1*2.1e8+1e-3*1.4e6), $2, $3, $4, $5}'
done >> H2_minus_X_equal_to_H.csv
echo "" >> H2_minus_X_equal_to_H.csv
echo "" >> H2_minus_X_equal_to_H.csv
echo "# 1 mM" >> H2_minus_X_equal_to_H.csv
for i in 0 1 2 3
  do
    c=${concentration[$i]}
    echo $c `sed -n '500,500p' /Applications/TOPAS/openTOPAS/openTOPAS_nBio_dev/TOPAS-nBio-dev/regression/GvalueIRT_H/run/2024Mar8/mainPython/0/HCO2m_"$c"_M_Br_1e-3_M_NO3m_1e-3_M.phsp` | awk '{print 1e12/($1*2.1e8+1e-3*1.4e6), $2, $3, $4, $5}'
done  >> H2_minus_X_equal_to_H.csv



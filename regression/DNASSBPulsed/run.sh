#!/bin/bash
# CONV
#
# Considerations:
# Dose rate:
#   0.12-18.9 Gy/min -> Milligan et al.
#   Higher dose rate in a FFF beam is 2400 MU/Min -> 24 Gy/Min at calibration point
#   Using a value in between: 10 Gy/Min 
# Dose:
#   Use 5 Gy dose, value used by Peter Wardman calculation, clincial consideration and
#   simulation statistics. 
# Time:
#   5 Gy -> 30 s @ 10 Gy/Min

np=(1)
time=(5) #125) # in seconds. Upper limit +5 seg
dose=(20) # in Gy
mean=(2) #0.5e-6) #60) # in seconds.
fwhm=(4) #1.0e-6) #120) # in seconds.

for i in 0 
  do
    echo "Pulse: " ${np[$i]} "Time limit (s): " ${time[$i]} "Dose (Gy): " ${dose[$i]} " mean (s): " ${mean[$i]} " FWHM (s): " ${fwhm[$i]}
    tcsh submitLocally.sh 0 10 ${dose[$i]}"Gy" ${dose[$i]} ${np[$i]} ${time[$i]} ${mean[$i]} ${fwhm[$i]}
done

############################
### Import Libraries

import numpy as np
import sys
from os import system as command

############################
### Define Water Density Function

def GetDensity(T_Celsius):
    return 0.999 + 1.094e-4  * T_Celsius \
                 - 7.397e-6  * T_Celsius**2 \
                 + 2.693e-8  * T_Celsius**3 \
                 - 4.714e-11 * T_Celsius**4

############################
### Define Main

def Main():
    TOPAS = sys.argv[1]
    Temperatures = np.array([0, 10, 20, 30, 40, 50, 60, 70 ,80, 90, 100])
    Densities    = GetDensity(Temperatures)

    for i in range(len(Temperatures)):
        command("sed 's/fDensity/'%s'/g' depFile1.tps > run1.tps"%(Densities[i]))
        command("sed 's/fTemperature/'%s'/g' run1.tps > run2.tps"%(Temperatures[i]))

        command("%s run2.tps"%(TOPAS))
        command("rm run1.tps")
        command("rm run2.tps")

############################
### Run Main

if __name__=="__main__":
    Main()

############################
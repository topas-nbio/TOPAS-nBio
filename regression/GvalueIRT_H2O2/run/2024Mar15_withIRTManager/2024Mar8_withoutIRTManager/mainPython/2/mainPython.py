############################
### Import Libraries

import numpy as np
import sys 
from os import system as command

############################
### Define Main

def Main():
    TOPAS = sys.argv[1]
    Concentrations = ["1e-3", "1e-2", "1e-1", "1e-0", "1e+1"]
    #Concentrations = ["e-4",             , "e-1", "e-0", "e+1", "e+2"]

    for c in Concentrations:
        command("sed 's/theConcentration/'%s'/g' mainTopas.txt > run1.txt" % c)
        command("%s run1.txt"%(TOPAS))
        command("mv run1.txt run_%s.txt"%(c))

############################
### Run Main

if __name__=="__main__":
    Main()

############################


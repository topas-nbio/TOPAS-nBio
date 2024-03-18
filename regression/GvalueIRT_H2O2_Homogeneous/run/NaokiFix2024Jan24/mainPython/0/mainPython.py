############################
### Import Libraries

import numpy as np
import sys 
from os import system as command

############################
### Define Main

def Main():
    TOPAS = sys.argv[1]
    Concentrations = ["e-4", "e-3", "e-2", "e-1", "e-0", "e+1", "e+2"]
    Concentrations = ["e-4",                      "e-0", "e+1", "e+2"]

    for c in Concentrations:
        command("sed 's/theConcentration/0.1031'%s'/g' mainTopas.txt > run1.txt" % c)

        command("%s run1.txt"%(TOPAS))

############################
### Run Main

if __name__=="__main__":
    Main()

############################


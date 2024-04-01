############################
### Import Libraries

import numpy as np
import sys
from os import system as command

############################
### Define Main

def Main():
    TOPAS = sys.argv[1]

    DMSO = ["1e-4", "1e-3", "1e-2", "1e-1", "1e-0"]
    repeats = [200, 200, 200, 200, 200, 200]
    kObs    = [ 1.32e7 * (7.1e9 * float(c) )**0.29 for c in DMSO]
     
    for k, c, r in zip(kObs, DMSO, repeats):
        command("sed 's/fOHDNAKobs/'%1.3e'/g' depFile1.txt | sed 's/fDMSOConcentration/'%s'/g' | sed 's/fNumberOfRepeats/'%d'/g' > run.txt " % (k, c, r))
        command("%s run.txt"%(TOPAS))

############################
### Run Main

if __name__=="__main__":
    Main()

############################

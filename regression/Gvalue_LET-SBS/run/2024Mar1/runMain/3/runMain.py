###############################
### Import Libraries

import os
import numpy

###############################
### Declare Globals

Protons = {"Energies"   :[5E2, 1E4, 1E5],\
           "Percentages":[1.0, 1.0, 1.0],\
           "TrackLength":[2.5,   5,  20],\
           "Histories"  :[  5,   5,   5]}

Electrons = {"Energies"   :[  2.0,  30.0, 999.999],\
             "Percentages":[  0.6,  0.20,    0.01],\
             "TrackLength":[100.0, 100.0,   100.0],\
             "Histories"  :[  100,   100,     100]}

Alphas = {"Energies"   :[4E3, 2E4, 5E4],\
		  "Percentages":[1.0, 1.0, 1.0],\
          "TrackLength":[2.5, 2.5, 2.5],\
          "Histories"  :[  3,   3,   3]}

###############################
### Declare Function

def ChangeFileAndRun(Particle, Energy, Histories, TrackLength, Percentage):
	os.system("sed 's/TheEnergy/'%s'/g' depFile1.txt > runstep1.txt"%(Energy))
	os.system("sed 's/TheParticle/'%s'/g' runstep1.txt > runstep2.txt"%(Particle))
	os.system("sed 's/TheHistories1/'%s'/g' runstep2.txt > runstep3.txt"%(int(Histories)))
	os.system("sed 's/TheHistories2/'%s'/g' runstep3.txt > runstep4.txt"%(int(Histories/10)))
	os.system("sed 's/ThePercentage/'%s'/g' runstep4.txt > runstep5.txt"%(Percentage))
	os.system("sed 's/TheTrackLength/'%s'/g' runstep5.txt > runstep6.txt"%(TrackLength))

	if (Percentage < 1):
		os.system("sed 's/TheKillMethod/'%s'/g' runstep6.txt > RUN.txt"%("True"))

	else:
		os.system("sed 's/TheKillMethod/'%s'/g' runstep6.txt > RUN.txt"%("False"))

	os.system("rm runstep1.txt")
	os.system("rm runstep2.txt")
	os.system("rm runstep3.txt")
	os.system("rm runstep4.txt")
	os.system("rm runstep5.txt")
	os.system("rm runstep6.txt")
	os.system("nBio_dev RUN.txt")
	os.system("rm RUN.txt")

###############################
### Declare Main

def Main():
	for i in range(len(Electrons["Energies"])):
		Energy      = Electrons["Energies"][i]
		Histories   = Electrons["Histories"][i]
		TrackLength = Electrons["TrackLength"][i]
		Percentage  = Electrons["Percentages"][i]
		ChangeFileAndRun("e-", Energy, Histories, TrackLength, Percentage)

	for i in range(len(Protons["Energies"])):
		Energy      = Protons["Energies"][i]
		Histories   = Protons["Histories"][i]
		TrackLength = Protons["TrackLength"][i]
		Percentage  = Protons["Percentages"][i]
		ChangeFileAndRun("proton", Energy, Histories, TrackLength, Percentage)

	for i in range(len(Alphas["Energies"])):
		Energy      = Alphas["Energies"][i]
		Histories   = Alphas["Histories"][i]
		TrackLength = Alphas["TrackLength"][i]
		Percentage  = Alphas["Percentages"][i] 
		ChangeFileAndRun("alpha", Energy, Histories, TrackLength, Percentage)

###############################
### Run Main

if __name__=="__main__":
	Main()

###############################
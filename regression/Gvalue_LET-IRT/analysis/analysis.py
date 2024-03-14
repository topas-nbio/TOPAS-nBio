#!/usr/bin/env python

####################################################
### Import Libraries

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os, subprocess
from os.path import isdir, join, split
from glob import glob
from functools import reduce
from pylab import rcParams

####################################################
### Declare Globals

rcParams['legend.numpoints'] = 1

ElectronInfo = ["2.0","30.0","999.999"]
ElectronLET  = np.array([[9.28044,0.215123],[1.37702,0.0391928],[0.195033,0.0121145]])

ProtonInfo   = ["500.0","10000.0","100000.0"]
ProtonLET    = np.array([[40.2692,0.32277],[4.49355,0.109095],[0.695069,0.0417797]])

AlphaInfo    = ["4000.0","20000.0","50000.0"]
AlphaLET     = np.array([[97.3568,0.868983],[28.3289,0.336693],[13.4607,0.32049]])

####################################################
### Define Read Functions

def ReadExperimental(File):
    Results = []
    f = open(File,"r")
    for z in f.readlines():
        A = z.split()
        if A[0] == "#":
            continue

        LET = float(A[0])
        GV  = float(A[1])
        GE  = float(A[2])

        Results.append(np.array([LET,GV,GE]))
    return np.array(Results)

def GetGValues(file,molecules):
    f = open(file)
    Result = {}
    for z in f.readlines():
        A = z.split()

        GValue = float(A[0])
        STD    = float(A[1])
        Time   = float(A[2])
        Name   = A[3]

        if Name in molecules and Time == 1E6:
            Result[Name] = [[GValue, STD]]

    return Result

def GetMeanGValue(Dir,Repeats,Molecules):
    ElectronRes = {}
    ProtonRes   = {}
    AlphaRes    = {}

    for i in range(len(ElectronInfo)):
        GValues_Holder = {}
        Energy = ElectronInfo[i]
        for j in range(len(Repeats)):
            File    = Dir+"/"+Repeats[j]+"/e-_"+Energy+"_GValue.phsp"
            GValues = GetGValues(File,Molecules)
            
            if GValues_Holder == {}:
                GValues_Holder = GValues

            else:
                for k in GValues:
                    GValues_Holder[k].append(GValues[k][0])

        for j in Molecules:
            DATA = np.array(GValues_Holder[j])
            MEAN = np.mean(DATA[:,0])
            STD  = sum(np.sqrt(DATA[:,1]**2))/len(DATA)

            if j not in ElectronRes:
                ElectronRes[j] = [[float(Energy),MEAN,STD]]

            else:
                ElectronRes[j].append([float(Energy),MEAN,STD])

    for i in range(len(ProtonInfo)):
        GValues_Holder = {}
        Energy = ProtonInfo[i]
        for j in range(len(Repeats)):
            File    = Dir+"/"+Repeats[j]+"/proton_"+Energy+"_GValue.phsp"
            GValues = GetGValues(File,Molecules)

            if GValues_Holder == {}:
                GValues_Holder = GValues

            else:
                for k in GValues:
                    GValues_Holder[k].append(GValues[k][0])

        for j in Molecules:
            DATA = np.array(GValues_Holder[j])
            MEAN = np.mean(DATA[:,0])
            STD  = sum(np.sqrt(DATA[:,1]**2))/len(DATA)

            if j not in ProtonRes:
                ProtonRes[j] = [[float(Energy),MEAN,STD]]

            else:
                ProtonRes[j].append([float(Energy),MEAN,STD])

    for i in range(len(AlphaInfo)):
        GValues_Holder = {}
        Energy = AlphaInfo[i]
        for j in range(len(Repeats)):
            File    = Dir+"/"+Repeats[j]+"/alpha_"+Energy+"_GValue.phsp"
            GValues = GetGValues(File,Molecules)

            if GValues_Holder == {}:
                GValues_Holder = GValues

            else:
                for k in GValues:
                    GValues_Holder[k].append(GValues[k][0])

        for j in Molecules:
            DATA = np.array(GValues_Holder[j])
            MEAN = np.mean(DATA[:,0])
            STD  = sum(np.sqrt(DATA[:,1]**2))/len(DATA)

            if j not in AlphaRes:
                AlphaRes[j] = [[float(Energy),MEAN,STD]]

            else:
                AlphaRes[j].append([float(Energy),MEAN,STD])

    for i in Molecules:
        ElectronRes[i] = np.array(ElectronRes[i])
        ProtonRes[i]   = np.array(ProtonRes[i])
        AlphaRes[i]    = np.array(AlphaRes[i])

    return ElectronRes, ProtonRes, AlphaRes


def GetMeanTime(Dir, Repeats):
    RealTimes = []
    UserTimes = []
    SysTimes  = []
    for  i in Repeats:
        File = Dir+"/"+i+"/log.out"
        Times = subprocess.check_output("tail -1 "+File, shell=True)
        Times  = Times.split()
        RealTimes.append(float(Times[0]))
        UserTimes.append(float(Times[2]))
        SysTimes.append(float(Times[4]))

    RealTimes = np.array(RealTimes)
    UserTimes = np.array(UserTimes)
    SysTimes  = np.array(SysTimes)

    return [np.mean(RealTimes),np.std(RealTimes),np.mean(UserTimes),np.std(UserTimes),np.mean(SysTimes),np.std(SysTimes)]

####################################################
### Define Plot Function

def plot_results(sut_dir, ref_dir, args):
    tag = split(sut_dir.rstrip(os.sep))[-1]

    molecules = ['OH^0','e_aq^-1','H3O^1','H2O2^0','H_2^0','H^0']
    moleculesName = ['$^{.}OH$','$e-_{aq}$','$H_{3}O^{+}$','$H_{2}O_{2}$','$H_{2}$','$H^{.}$']

    Ref_Repeats = os.listdir(ref_dir)
    Sut_Repeats = os.listdir(sut_dir)

    ElectronRef, ProtonRef, AlphaRef = GetMeanGValue(ref_dir,Ref_Repeats,molecules)
    ElectronSut, ProtonSut, AlphaSut = GetMeanGValue(sut_dir,Sut_Repeats,molecules)

    RefTimes = GetMeanTime(ref_dir,Ref_Repeats)
    SutTimes = GetMeanTime(sut_dir,Sut_Repeats)

    Experimental = {"electron":{},"proton":{},"alpha":{}}

    Experimental["electron"]["OH"]   = {}
    Experimental["electron"]["eaq"]  = {}
    Experimental["electron"]["H"]    = {}
    Experimental["electron"]["H2"]   = {}
    Experimental["electron"]["H2O2"] = {}

    Experimental["proton"]["OH"]   = {}
    Experimental["proton"]["eaq"]  = {}
    Experimental["proton"]["H"]    = {}
    Experimental["proton"]["H2"]   = {}
    Experimental["proton"]["H2O2"] = {}

    Experimental["alpha"]["OH"]   = {}
    Experimental["alpha"]["eaq"]  = {}
    Experimental["alpha"]["H"]    = {}
    Experimental["alpha"]["H2"]   = {}
    Experimental["alpha"]["H2O2"] = {}

    Experimental["electron"]["OH"]["Appleby"]    = ReadExperimental("./analysis/benchmark//Electron/OH_Appleby1969.txt")
    Experimental["electron"]["OH"]["Burns"]      = ReadExperimental("./analysis/benchmark//Electron/OH_Burns1981.txt")
    Experimental["electron"]["eaq"]["Appleby"]   = ReadExperimental("./analysis/benchmark//Electron/eaq_Appleby1969.txt")
    Experimental["electron"]["H"]["Appleby"]     = ReadExperimental("./analysis/benchmark//Electron/H_Appleby1969.txt")
    Experimental["electron"]["H"]["Elliot"]      = ReadExperimental("./analysis/benchmark//Electron/H_Elliot1993.txt")
    Experimental["electron"]["H2"]["Appleby"]    = ReadExperimental("./analysis/benchmark//Electron/H2_Appleby1969.txt")
    Experimental["electron"]["H2"]["Crumiere"]   = ReadExperimental("./analysis/benchmark//Electron/H2_Crumiere2013.txt")
    Experimental["electron"]["H2O2"]["Appleby"]  = ReadExperimental("./analysis/benchmark//Electron/H2O2_Appleby1969.txt")
    Experimental["electron"]["H2O2"]["Wasselin"] = ReadExperimental("./analysis/benchmark//Electron/H2O2_Wasselin2002.txt")

    Experimental["proton"]["H2O2"]["Appleby"]  = ReadExperimental("./analysis/benchmark//Proton/H2O2_Appleby1969.txt")
    Experimental["proton"]["H2O2"]["Pastina"]  = ReadExperimental("./analysis/benchmark//Proton/H2O2_Pastina1999.txt")
    Experimental["proton"]["H2O2"]["Wasselin"] = ReadExperimental("./analysis/benchmark//Proton/H2O2_Wasselin2002.txt")
    Experimental["proton"]["H2"]["Appleby"]    = ReadExperimental("./analysis/benchmark//Proton/H2_Appleby1969.txt")
    Experimental["proton"]["H"]["Appleby"]     = ReadExperimental("./analysis/benchmark//Proton/H_Appleby1969.txt")
    Experimental["proton"]["OH"]["Anderson"]   = ReadExperimental("./analysis/benchmark//Proton/OH_Anderson1969.txt")
    Experimental["proton"]["OH"]["Appleby"]    = ReadExperimental("./analysis/benchmark//Proton/OH_Appleby1969.txt")
    Experimental["proton"]["OH"]["Burns"]      = ReadExperimental("./analysis/benchmark//Proton/OH_Burns1981.txt")
    Experimental["proton"]["eaq"]["Appleby"]   = ReadExperimental("./analysis/benchmark//Proton/eaq_Appleby1969.txt")
    Experimental["proton"]["eaq"]["Sauer"]     = ReadExperimental("./analysis/benchmark//Proton/eaq_Sauer1977.txt")

    Experimental["alpha"]["H2O2"]["Appleby"]  = ReadExperimental("./analysis/benchmark//Alpha/H2O2_Appleby1969.txt")
    Experimental["alpha"]["H2O2"]["Pastina"]  = ReadExperimental("./analysis/benchmark//Alpha/H2O2_Pastina1993.txt")
    Experimental["alpha"]["H2O2"]["Wasselin"] = ReadExperimental("./analysis/benchmark//Alpha/H2O2_Wasselin2002.txt")
    Experimental["alpha"]["H2"]["Appleby"]    = ReadExperimental("./analysis/benchmark//Alpha/H2_Appleby1969.txt")
    Experimental["alpha"]["H2"]["Crumiere"]   = ReadExperimental("./analysis/benchmark//Alpha/H2_Crumiere2013.txt")
    Experimental["alpha"]["H"]["Appleby"]     = ReadExperimental("./analysis/benchmark//Alpha/H_Appleby1969.txt")
    Experimental["alpha"]["OH"]["Anderson"]   = ReadExperimental("./analysis/benchmark//Alpha/OH_Anderson1969.txt")
    Experimental["alpha"]["OH"]["Appleby"]    = ReadExperimental("./analysis/benchmark//Alpha/OH_Appleby1969.txt")
    Experimental["alpha"]["OH"]["Burns"]      = ReadExperimental("./analysis/benchmark//Alpha/OH_Burns1981.txt")
    Experimental["alpha"]["eaq"]["Anderson"]  = ReadExperimental("./analysis/benchmark//Alpha/eaq_Anderson1969.txt")
    Experimental["alpha"]["eaq"]["Burns"]     = ReadExperimental("./analysis/benchmark//Alpha/eaq_Burns1981.txt")

    plt.rcParams['legend.title_fontsize'] = 15

    fig  = plt.figure(figsize=(33,27))
    ax1  = plt.subplot2grid((3,3),(0,0))
    ax2  = plt.subplot2grid((3,3),(0,1))
    ax3  = plt.subplot2grid((3,3),(0,2))
    ax4  = plt.subplot2grid((3,3),(1,0))
    ax5  = plt.subplot2grid((3,3),(1,1))
    ax6  = plt.subplot2grid((3,3),(1,2))
    ax7  = plt.subplot2grid((3,3),(2,0),colspan=3,rowspan=1)


    ## H2O2
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.tick_params(axis='both', which='minor', labelsize=20)
    ax1.set_xlabel(r"$LET$ $keV/\mu m$",fontsize=24)
    ax1.set_ylabel(r"GValue $(Molecules/100eV)$",fontsize=24)
    ax1.errorbar(ElectronLET[:,0],ElectronSut["H2O2^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronSut["H2O2^0"][:,2], linewidth=2, fmt="b--", label="TOPAS-Sut e-")
    ax1.errorbar(ProtonLET[:,0],  ProtonSut["H2O2^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonSut["H2O2^0"][:,2],   linewidth=2, fmt="b-.", label="TOPAS-Sut Proton")
    ax1.errorbar(AlphaLET[:,0],   AlphaSut["H2O2^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaSut["H2O2^0"][:,2],    linewidth=2, fmt="b:", label="TOPAS-Sut Alpha")
    ax1.errorbar(ElectronLET[:,0],ElectronRef["H2O2^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronRef["H2O2^0"][:,2], linewidth=2, fmt="r--", label="TOPAS-Ref e-")
    ax1.errorbar(ProtonLET[:,0],  ProtonRef["H2O2^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonRef["H2O2^0"][:,2],   linewidth=2, fmt="r-.", label="TOPAS-Ref Proton")
    ax1.errorbar(AlphaLET[:,0],   AlphaRef["H2O2^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaRef["H2O2^0"][:,2],    linewidth=2, fmt="r:", label="TOPAS-Ref Alpha")
    ax1.errorbar(Experimental["electron"]["H2O2"]["Appleby"][:,0], Experimental["electron"]["H2O2"]["Appleby"][:,1], yerr=Experimental["electron"]["H2O2"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax1.errorbar(Experimental["electron"]["H2O2"]["Wasselin"][:,0],Experimental["electron"]["H2O2"]["Wasselin"][:,1],yerr=Experimental["electron"]["H2O2"]["Wasselin"][:,2],fmt="p",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax1.errorbar(Experimental["proton"]["H2O2"]["Appleby"][:,0],   Experimental["proton"]["H2O2"]["Appleby"][:,1],   yerr=Experimental["proton"]["H2O2"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax1.errorbar(Experimental["proton"]["H2O2"]["Pastina"][:,0],   Experimental["proton"]["H2O2"]["Pastina"][:,1],   yerr=Experimental["proton"]["H2O2"]["Pastina"][:,2],fmt="v",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax1.errorbar(Experimental["proton"]["H2O2"]["Wasselin"][:,0],  Experimental["proton"]["H2O2"]["Wasselin"][:,1],  yerr=Experimental["proton"]["H2O2"]["Wasselin"][:,2],fmt="p",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax1.errorbar(Experimental["alpha"]["H2O2"]["Appleby"][:,0],    Experimental["alpha"]["H2O2"]["Appleby"][:,1],    yerr=Experimental["alpha"]["H2O2"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax1.errorbar(Experimental["alpha"]["H2O2"]["Pastina"][:,0],    Experimental["alpha"]["H2O2"]["Pastina"][:,1],    yerr=Experimental["alpha"]["H2O2"]["Pastina"][:,2],fmt="v",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax1.errorbar(Experimental["alpha"]["H2O2"]["Wasselin"][:,0],   Experimental["alpha"]["H2O2"]["Wasselin"][:,1],   yerr=Experimental["alpha"]["H2O2"]["Wasselin"][:,2],fmt="p",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax1.set_xlim(0.1,120)
    ax1.set_ylim(0.5,1.3)
    ax1.set_xscale("log")
    ax1.grid(True,dashes=[5,5])
    ax1.legend(loc=2,fontsize=20,title=r"$H_{2}O_{2}$")

    # eaq
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='minor', labelsize=20)
    ax2.set_xlabel(r"$LET$ $keV/\mu m$",fontsize=24)
    ax2.set_ylabel(r"GValue $(Molecules/100eV)$",fontsize=24)
    ax2.errorbar(ElectronLET[:,0],ElectronSut["e_aq^-1"][:,1],xerr=ElectronLET[:,1],yerr=ElectronSut["e_aq^-1"][:,2], linewidth=2, fmt="b--", label="TOPAS-Sut e-")
    ax2.errorbar(ProtonLET[:,0],  ProtonSut["e_aq^-1"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonSut["e_aq^-1"][:,2],   linewidth=2, fmt="b-.", label="TOPAS-Sut Proton")
    ax2.errorbar(AlphaLET[:,0],   AlphaSut["e_aq^-1"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaSut["e_aq^-1"][:,2],    linewidth=2, fmt="b:", label="TOPAS-Sut Alpha")
    ax2.errorbar(ElectronLET[:,0],ElectronRef["e_aq^-1"][:,1],xerr=ElectronLET[:,1],yerr=ElectronRef["e_aq^-1"][:,2], linewidth=2, fmt="r--", label="TOPAS-Ref e-")
    ax2.errorbar(ProtonLET[:,0],  ProtonRef["e_aq^-1"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonRef["e_aq^-1"][:,2],   linewidth=2, fmt="r-.", label="TOPAS-Ref Proton")
    ax2.errorbar(AlphaLET[:,0],   AlphaRef["e_aq^-1"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaRef["e_aq^-1"][:,2],    linewidth=2, fmt="r:", label="TOPAS-Ref Alpha")
    ax2.errorbar(Experimental["electron"]["eaq"]["Appleby"][:,0],Experimental["electron"]["eaq"]["Appleby"][:,1],yerr=Experimental["electron"]["eaq"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax2.errorbar(Experimental["proton"]["eaq"]["Appleby"][:,0],Experimental["proton"]["eaq"]["Appleby"][:,1],yerr=Experimental["proton"]["eaq"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax2.errorbar(Experimental["proton"]["eaq"]["Sauer"][:,0],Experimental["proton"]["eaq"]["Sauer"][:,1],yerr=Experimental["proton"]["eaq"]["Sauer"][:,2],fmt="d",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax2.errorbar(Experimental["alpha"]["eaq"]["Anderson"][:,0],Experimental["alpha"]["eaq"]["Anderson"][:,1],yerr=Experimental["alpha"]["eaq"]["Anderson"][:,2],fmt="^",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax2.errorbar(Experimental["alpha"]["eaq"]["Burns"][:,0],Experimental["alpha"]["eaq"]["Burns"][:,1],yerr=Experimental["alpha"]["eaq"]["Burns"][:,2],fmt="s",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax2.set_xlim(0.1,120)
    ax2.set_ylim(0,3)
    ax2.set_xscale("log")
    ax2.grid(True,dashes=[5,5])
    ax2.legend(loc=1,fontsize=20,title=r"$e_{aq}^{-}$")

    # OH
    ax3.tick_params(axis='both', which='major', labelsize=20)
    ax3.tick_params(axis='both', which='minor', labelsize=20)
    ax3.set_xlabel(r"$LET$ $keV/\mu m$",fontsize=24)
    ax3.set_ylabel(r"GValue $(Molecules/100eV)$",fontsize=24)
    ax3.errorbar(ElectronLET[:,0],ElectronSut["OH^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronSut["OH^0"][:,2], linewidth=2, fmt="b--", label="TOPAS-Sut e-")
    ax3.errorbar(ProtonLET[:,0],  ProtonSut["OH^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonSut["OH^0"][:,2],   linewidth=2, fmt="b-.", label="TOPAS-Sut Proton")
    ax3.errorbar(AlphaLET[:,0],   AlphaSut["OH^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaSut["OH^0"][:,2],    linewidth=2, fmt="b:", label="TOPAS-Sut Alpha")
    ax3.errorbar(ElectronLET[:,0],ElectronRef["OH^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronRef["OH^0"][:,2], linewidth=2, fmt="r--", label="TOPAS-Ref e-")
    ax3.errorbar(ProtonLET[:,0],  ProtonRef["OH^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonRef["OH^0"][:,2],   linewidth=2, fmt="r-.", label="TOPAS-Ref Proton")
    ax3.errorbar(AlphaLET[:,0],   AlphaRef["OH^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaRef["OH^0"][:,2],    linewidth=2, fmt="r:", label="TOPAS-Ref Alpha")
    ax3.errorbar(Experimental["electron"]["OH"]["Appleby"][:,0],Experimental["electron"]["OH"]["Appleby"][:,1],yerr=Experimental["electron"]["OH"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax3.errorbar(Experimental["electron"]["OH"]["Burns"][:,0],Experimental["electron"]["OH"]["Burns"][:,1],yerr=Experimental["electron"]["OH"]["Burns"][:,2],fmt="s",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax3.errorbar(Experimental["proton"]["OH"]["Anderson"][:,0],Experimental["proton"]["OH"]["Anderson"][:,1],yerr=Experimental["proton"]["OH"]["Anderson"][:,2],fmt="^",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax3.errorbar(Experimental["proton"]["OH"]["Appleby"][:,0],Experimental["proton"]["OH"]["Appleby"][:,1],yerr=Experimental["proton"]["OH"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax3.errorbar(Experimental["proton"]["OH"]["Burns"][:,0],Experimental["proton"]["OH"]["Burns"][:,1],yerr=Experimental["proton"]["OH"]["Burns"][:,2],fmt="s",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax3.errorbar(Experimental["alpha"]["OH"]["Anderson"][:,0],Experimental["alpha"]["OH"]["Anderson"][:,1],yerr=Experimental["alpha"]["OH"]["Anderson"][:,2],fmt="^",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax3.errorbar(Experimental["alpha"]["OH"]["Appleby"][:,0],Experimental["alpha"]["OH"]["Appleby"][:,1],yerr=Experimental["alpha"]["OH"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax3.errorbar(Experimental["alpha"]["OH"]["Burns"][:,0],Experimental["alpha"]["OH"]["Burns"][:,1],yerr=Experimental["alpha"]["OH"]["Burns"][:,2],fmt="s",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax3.set_xlim(0.1,120)
    ax3.set_ylim(0,3)
    ax3.set_xscale("log")
    ax3.grid(True,dashes=[5,5])
    ax3.legend(loc=3,fontsize=20,title=r"$OH^{\bullet}$")

    # H
    ax4.tick_params(axis='both', which='major', labelsize=20)
    ax4.tick_params(axis='both', which='minor', labelsize=20)
    ax4.set_xlabel(r"$LET$ $keV/\mu m$",fontsize=24)
    ax4.set_ylabel(r"GValue $(Molecules/100eV)$",fontsize=24)
    ax4.errorbar(ElectronLET[:,0],ElectronSut["H^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronSut["OH^0"][:,2], linewidth=2, fmt="b--", label="TOPAS-Sut e-")
    ax4.errorbar(ProtonLET[:,0],  ProtonSut["H^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonSut["OH^0"][:,2],   linewidth=2, fmt="b-.", label="TOPAS-Sut Proton")
    ax4.errorbar(AlphaLET[:,0],   AlphaSut["H^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaSut["OH^0"][:,2],    linewidth=2, fmt="b:", label="TOPAS-Sut Alpha")
    ax4.errorbar(ElectronLET[:,0],ElectronRef["H^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronRef["OH^0"][:,2], linewidth=2, fmt="r--", label="TOPAS-Ref e-")
    ax4.errorbar(ProtonLET[:,0],  ProtonRef["H^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonRef["OH^0"][:,2],   linewidth=2, fmt="r-.", label="TOPAS-Ref Proton")
    ax4.errorbar(AlphaLET[:,0],   AlphaRef["H^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaRef["OH^0"][:,2],    linewidth=2, fmt="r:", label="TOPAS-Ref Alpha")    
    ax4.errorbar(Experimental["electron"]["H"]["Appleby"][:,0],Experimental["electron"]["H"]["Appleby"][:,1],yerr=Experimental["electron"]["H"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax4.errorbar(Experimental["electron"]["H"]["Elliot"][:,0],Experimental["electron"]["H"]["Elliot"][:,1],yerr=Experimental["electron"]["H"]["Elliot"][:,2],fmt="*",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax4.errorbar(Experimental["proton"]["H"]["Appleby"][:,0],Experimental["proton"]["H"]["Appleby"][:,1],yerr=Experimental["proton"]["H"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax4.errorbar(Experimental["alpha"]["H"]["Appleby"][:,0],Experimental["alpha"]["H"]["Appleby"][:,1],yerr=Experimental["alpha"]["H"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax4.set_xlim(0.1,120)
    ax4.set_ylim(0.1,0.8)
    ax4.set_xscale("log")
    ax4.grid(True,dashes=[5,5])
    ax4.legend(loc=3,fontsize=20,title=r"$H^{\bullet}$")

    # H2
    ax5.tick_params(axis='both', which='major', labelsize=20)
    ax5.tick_params(axis='both', which='minor', labelsize=20)
    ax5.set_xlabel(r"$LET$ $keV/\mu m$",fontsize=24)
    ax5.set_ylabel(r"GValue $(Molecules/100eV)$",fontsize=24)
    ax5.errorbar(ElectronLET[:,0],ElectronSut["H_2^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronSut["H_2^0"][:,2], linewidth=2, fmt="b--", label="TOPAS-Sut e-")
    ax5.errorbar(ProtonLET[:,0],  ProtonSut["H_2^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonSut["H_2^0"][:,2],   linewidth=2, fmt="b-.", label="TOPAS-Sut Proton")
    ax5.errorbar(AlphaLET[:,0],   AlphaSut["H_2^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaSut["H_2^0"][:,2],    linewidth=2, fmt="b:", label="TOPAS-Sut Alpha")
    ax5.errorbar(ElectronLET[:,0],ElectronRef["H_2^0"][:,1],xerr=ElectronLET[:,1],yerr=ElectronRef["H_2^0"][:,2], linewidth=2, fmt="r--", label="TOPAS-Ref e-")
    ax5.errorbar(ProtonLET[:,0],  ProtonRef["H_2^0"][:,1],  xerr=ProtonLET[:,1],  yerr=ProtonRef["H_2^0"][:,2],   linewidth=2, fmt="r-.", label="TOPAS-Ref Proton")
    ax5.errorbar(AlphaLET[:,0],   AlphaRef["H_2^0"][:,1],   xerr=AlphaLET[:,1],   yerr=AlphaRef["H_2^0"][:,2],    linewidth=2, fmt="r:", label="TOPAS-Ref Alpha")      
    ax5.errorbar(Experimental["electron"]["H2"]["Appleby"][:,0],Experimental["electron"]["H2"]["Appleby"][:,1],yerr=Experimental["electron"]["H2"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax5.errorbar(Experimental["electron"]["H2"]["Crumiere"][:,0],Experimental["electron"]["H2"]["Crumiere"][:,1],yerr=Experimental["electron"]["H2"]["Crumiere"][:,2],fmt="H",markersize=15,markerfacecolor="None",markeredgewidth=2,color="black")
    ax5.errorbar(Experimental["proton"]["H2"]["Appleby"][:,0],Experimental["proton"]["H2"]["Appleby"][:,1],yerr=Experimental["proton"]["H2"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="red")
    ax5.errorbar(Experimental["alpha"]["H2"]["Appleby"][:,0],Experimental["alpha"]["H2"]["Appleby"][:,1],yerr=Experimental["alpha"]["H2"]["Appleby"][:,2],fmt="o",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax5.errorbar(Experimental["alpha"]["H2"]["Crumiere"][:,0],Experimental["alpha"]["H2"]["Crumiere"][:,1],yerr=Experimental["alpha"]["H2"]["Crumiere"][:,2],fmt="H",markersize=15,markerfacecolor="None",markeredgewidth=2,color="blue")
    ax5.set_xlim(0.1,120)
    ax5.set_ylim(0.2,1.4)
    ax5.set_xscale("log")
    ax5.grid(True,dashes=[5,5])
    ax5.legend(loc=2,fontsize=20,title=r"$H_{2}$")

    # Ox-Red
    ElectronOx_Sut  = ElectronSut["OH^0"][:,1] + (2*ElectronSut["H2O2^0"][:,1])
    ElectronRed_Sut = ElectronSut["e_aq^-1"][:,1] + (2*ElectronSut["H_2^0"][:,1]) + ElectronSut["H^0"][:,1]
    ElectronOx_Sut_Err  = np.sqrt(ElectronSut["OH^0"][:,2]**2 + (2*ElectronSut["H2O2^0"][:,2])**2)
    ElectronRed_Sut_Err = np.sqrt(ElectronSut["e_aq^-1"][:,2]**2 + (2*ElectronSut["H_2^0"][:,2])**2 + ElectronSut["H^0"][:,2]**2)
    ElectronRatio_Sut   = ElectronOx_Sut/ElectronRed_Sut
    ElectronRatio_Sut_Err = ElectronRatio_Sut * np.sqrt((ElectronOx_Sut_Err/ElectronOx_Sut)**2 + (ElectronRed_Sut_Err/ElectronRed_Sut)**2)

    ElectronOx_Ref  = ElectronRef["OH^0"][:,1] + (2*ElectronRef["H2O2^0"][:,1])
    ElectronRed_Ref = ElectronRef["e_aq^-1"][:,1] + (2*ElectronRef["H_2^0"][:,1]) + ElectronRef["H^0"][:,1]
    ElectronOx_Ref_Err  = np.sqrt(ElectronRef["OH^0"][:,2]**2 + (2*ElectronRef["H2O2^0"][:,2])**2)
    ElectronRed_Ref_Err = np.sqrt(ElectronRef["e_aq^-1"][:,2]**2 + (2*ElectronRef["H_2^0"][:,2])**2 + ElectronRef["H^0"][:,2]**2)
    ElectronRatio_Ref   = ElectronOx_Ref/ElectronRed_Ref
    ElectronRatio_Ref_Err = ElectronRatio_Ref * np.sqrt((ElectronOx_Ref_Err/ElectronOx_Ref)**2 + (ElectronRed_Ref_Err/ElectronRed_Ref)**2)

    ProtonOx_Sut    = ProtonSut["OH^0"][:,1] + (2*ProtonSut["H2O2^0"][:,1])
    ProtonRed_Sut   = ProtonSut["e_aq^-1"][:,1] + (2*ProtonSut["H_2^0"][:,1]) + ProtonSut["H^0"][:,1]
    ProtonOx_Sut_Err  = np.sqrt(ProtonSut["OH^0"][:,1]**2 + (2*ProtonSut["H2O2^0"][:,1])**2)
    ProtonRed_Sut_Err = np.sqrt(ProtonSut["e_aq^-1"][:,1]**2 + (2*ProtonSut["H_2^0"][:,1])**2 + ProtonSut["H^0"][:,1]**2)
    ProtonRatio_Sut   = ProtonOx_Sut/ProtonRed_Sut
    ProtonRatio_Sut_Err = ProtonRatio_Sut * np.sqrt((ProtonOx_Sut_Err/ProtonOx_Sut)**2 + (ProtonRed_Sut_Err/ProtonRed_Sut)**2)

    ProtonOx_Ref    = ProtonRef["OH^0"][:,1] + (2*ProtonRef["H2O2^0"][:,1])
    ProtonRed_Ref   = ProtonRef["e_aq^-1"][:,1] + (2*ProtonRef["H_2^0"][:,1]) + ProtonRef["H^0"][:,1]
    ProtonOx_Ref_Err  = np.sqrt(ProtonRef["OH^0"][:,1]**2 + (2*ProtonRef["H2O2^0"][:,1])**2)
    ProtonRed_Ref_Err = np.sqrt(ProtonRef["e_aq^-1"][:,1]**2 + (2*ProtonRef["H_2^0"][:,1])**2 + ProtonRef["H^0"][:,1]**2)
    ProtonRatio_Ref   = ProtonOx_Ref/ProtonRed_Ref
    ProtonRatio_Ref_Err = ProtonRatio_Ref * np.sqrt((ProtonOx_Ref_Err/ProtonOx_Ref)**2 + (ProtonRed_Ref_Err/ProtonRed_Ref)**2)

    AlphaOx_Sut     = AlphaSut["OH^0"][:,1] + (2*AlphaSut["H2O2^0"][:,1])
    AlphaRed_Sut    = AlphaSut["e_aq^-1"][:,1] + (2*AlphaSut["H_2^0"][:,1]) + AlphaSut["H^0"][:,1]
    AlphaOx_Sut_Err  = np.sqrt(AlphaSut["OH^0"][:,1]**2 + (2*AlphaSut["H2O2^0"][:,1])**2)
    AlphaRed_Sut_Err = np.sqrt(AlphaSut["e_aq^-1"][:,1]**2 + (2*AlphaSut["H_2^0"][:,1])**2 + AlphaSut["H^0"][:,1]**2)
    AlphaRatio_Sut   = AlphaOx_Sut/AlphaRed_Sut
    AlphaRatio_Sut_Err = AlphaRatio_Sut * np.sqrt((AlphaOx_Sut_Err/AlphaOx_Sut)**2 + (AlphaRed_Sut_Err/AlphaRed_Sut)**2)

    AlphaOx_Ref     = AlphaRef["OH^0"][:,1] + (2*AlphaRef["H2O2^0"][:,1])
    AlphaRed_Ref    = AlphaRef["e_aq^-1"][:,1] + (2*AlphaRef["H_2^0"][:,1]) + AlphaRef["H^0"][:,1]
    AlphaOx_Ref_Err  = np.sqrt(AlphaRef["OH^0"][:,1]**2 + (2*AlphaRef["H2O2^0"][:,1])**2)
    AlphaRed_Ref_Err = np.sqrt(AlphaRef["e_aq^-1"][:,1]**2 + (2*AlphaRef["H_2^0"][:,1])**2 + AlphaRef["H^0"][:,1]**2)
    AlphaRatio_Ref   = AlphaOx_Ref/AlphaRed_Ref
    AlphaRatio_Ref_Err = AlphaRatio_Ref * np.sqrt((AlphaOx_Ref_Err/AlphaOx_Ref)**2 + (AlphaRed_Ref_Err/AlphaRed_Ref)**2)

    ax6.tick_params(axis='both', which='major', labelsize=20)
    ax6.tick_params(axis='both', which='minor', labelsize=20)
    ax6.set_xlabel(r"$LET$ $keV/\mu m$",fontsize=24)
    ax6.set_ylabel(r"Ox/Red Ratios",fontsize=24)
    ax6.errorbar(ElectronLET[:,0],ElectronRatio_Sut, yerr=ElectronRatio_Sut_Err*0, fmt="b--", linewidth=2, label="TOPAS-Sut e-")
    ax6.errorbar(ProtonLET[:,0],  ProtonRatio_Sut,   yerr=ProtonRatio_Sut_Err*0,   fmt="b-.", linewidth=2, label="TOPAS-Sut Proton")
    ax6.errorbar(AlphaLET[:,0],   AlphaRatio_Sut,    yerr=AlphaRatio_Sut_Err*0,    fmt="b:", linewidth=2, label="TOPAS-Sut Alpha")
    ax6.errorbar(ElectronLET[:,0],ElectronRatio_Ref, yerr=ElectronRatio_Ref_Err*0, fmt="r--", linewidth=2, label="TOPAS-Ref e-")
    ax6.errorbar(ProtonLET[:,0],  ProtonRatio_Ref,   yerr=ProtonRatio_Ref_Err*0,   fmt="r-.", linewidth=2, label="TOPAS-Ref Proton")
    ax6.errorbar(AlphaLET[:,0],   AlphaRatio_Ref,    yerr=AlphaRatio_Ref_Err*0,    fmt="r:", linewidth=2, label="TOPAS-Ref Alpha")
    ax6.set_xlim(0.1,120)
    ax6.set_ylim(0.95,1.05)
    ax6.set_xscale("log")
    ax6.grid(True,dashes=[5,5])
    ax6.legend(loc=2,fontsize=20)

    ax7.axis("off")
    Table = ax7.table(cellText=[['%1.3f +/- %1.3f'%(RefTimes[0],RefTimes[1]),'%1.3f +/- %1.3f'%(SutTimes[0],SutTimes[1])],\
                                ['%1.3f +/- %1.3f'%(RefTimes[2],RefTimes[3]),'%1.3f +/- %1.3f'%(SutTimes[2],SutTimes[3])],\
                                ['%1.3f +/- %1.3f'%(RefTimes[4],RefTimes[5]),'%1.3f +/- %1.3f'%(SutTimes[4],SutTimes[5])]],\
                      rowLabels=("Real","User","Sys"),\
                      colLabels=("TOPAS-Ref","TOPAS-Sut"), \
                      loc='center')

    Table.auto_set_font_size(False)
    Table.set_fontsize(24)
    Table.scale(1,3)

    plt.tight_layout()
    plt.savefig(join(args.outdir, 'Gvalue_LET-IRT.eps'),bbox_inches='tight')
    plt.savefig(join(args.outdir, 'Gvalue_LET-IRT.pdf'),bbox_inches='tight')


####################################################
### Define Main

def Main():
    parser = ArgumentParser()
    parser.add_argument('sut_dir', help='Result directory for MC under test')
    parser.add_argument('ref_dir', help='Result directory for benchmark')
    parser.add_argument('-o', '--outdir', default=None, help='Output directory')
    parser.add_argument('--sut_label', default='TOPAS under test')
    parser.add_argument('--ref_label', default='TOPAS benchmark')
    args = parser.parse_args()
    
    if not isdir(args.sut_dir):
        print('Could not find %s' % args.sut_dir)
        exit(1)

    if not isdir(args.ref_dir):
        print('Could not find %s' % args.ref_dir)
        exit(1)

    if not args.outdir:
        args.outdir = join('results', split(args.sut_dir.rstrip(os.sep))[-1])
    if not isdir(args.outdir):
        os.makedirs(args.outdir)

    plot_results(args.sut_dir, args.ref_dir, args)

####################################################
### Run Main

if __name__=="__main__":
    Main()

####################################################

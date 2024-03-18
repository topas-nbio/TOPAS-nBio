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
import matplotlib

####################################################
### Define Read Functions

def ReadGValues(File):
    GValues = {}

    f = open(File,"r")
    for z in f.readlines():
        A = z.split()
        if len(A) < 4:
            continue

        GVal = float(A[0])
        GErr = float(A[1])
        Time = float(A[2])
        Name = A[3]

        if Name in GValues:
            GValues[Name].append([Time,GVal,GErr])

        else:
            GValues[Name] = [[Time,GVal,GErr]]


    for i in GValues.keys():
        GValues[i] = np.array(GValues[i])

    return GValues


def Glob(match_path):
    os.system('ls ' + match_path + ' > list')
    listFiles = open('list','r')
    fnames = []
    for name in listFiles:
        name = name.split()[0] 
        fnames.append(name)
    os.system('rm list')
    return fnames

def average_results(match_path):
    fnames = Glob(match_path)
    if len(fnames) == 0:
        return None

    Result = {}
    n = 0.0
    for f in fnames:
        data = ReadGValues(f)
        if n == 0.0:
            Results = data
        else:
            for Mol in data.keys():
                Time1 = data[Mol][:,0]
                Data1 = data[Mol][:,1]
                Stdv1 = data[Mol][:,2]

                Results[Mol][:,1] += Data1
                Results[Mol][:,2] += Stdv1**2
        n += 1

    for i in Results.keys():
        Results[i] = np.array(Results[i])
        Results[i][:,1] /= n
        Results[i][:,2] = 1./n * np.sqrt(Results[i][:,2])

    return Results


def average_results_time(match_path):
    fnames = glob(match_path)
    if len(fnames) == 0:
        return None
    
    Real, RealE = 0.0, 0.0
    User, UserE = 0.0, 0.0
    Sys,  SysE  = 0.0, 0.0
    N = 0.0

    for f in fnames:
        Times = (subprocess.check_output("grep Execution " + f + " | awk '{print $2}' | sed 's/Real=//g' | sed 's/s//g' | awk '{sum+=$1}END{print sum/NR}'", shell=True))
        RealT = float(Times) #Times.split()[0])
        Real  += RealT
        RealE += RealT*RealT

        Times = (subprocess.check_output("grep Execution " + f + " | awk '{print $2}' | sed 's/User=//g' | sed 's/s//g' | awk '{sum+=$1}END{print sum/NR}'", shell=True))
        UserT = float(Times) #float(Times.split()[2])
        User  += UserT
        UserE += UserT*UserT

        Times = (subprocess.check_output("grep Execution " + f + " | awk '{print $2}' | sed 's/Sys=//g' | sed 's/s//g' | awk '{sum+=$1}END{print sum/NR}'", shell=True))
        SysT  = float(Times) #float(Times.split()[4])
        Sys  += SysT
        Sys  += SysT*SysT
        N += 1

    Real /= N
    User /= N
    Sys  /= N

    if N > 1:
        RealE = np.sqrt(N/(N-1) * (np.abs(RealE/N - Real**2)))
        UserE = np.sqrt(N/(N-1) * (np.abs(UserE/N - User**2)))
        SysE  = np.sqrt(N/(N-1) * (np.abs(SysE/N  - Sys**2)))

    else:
        RealE, UserE, SysE = 0, 0, 0

    return (Real,RealE,User,UserE,Sys,SysE)

####################################################
### Define Plot Function

def plot_results(sut_dir, ref_dir, args):
    concentrations = ["1e-3", "1e-2", "1e-1", "1e-0", "1e+1"]
    #concentrations = ["0.1031e-2", "0.1031e-1", "0.1031e-0", "0.1031e+1", "0.1031e+2"] 

    sut_T = average_results_time(sut_dir + '/*/log.out')
    ref_T = average_results_time(ref_dir + '/*/log.out')

    sut_G   = []
    ref_G   = []

    sut_SG   = {}
    ref_SG   = {}
    kobs = 9.7e8

    for T in concentrations: 
        sut_GT   = average_results(sut_dir + '/*/Methanol_%s_M_Nitrate_25e-3_M.phsp'%(T))
        ref_GT   = average_results(ref_dir + '/*/Methanol_%s_M_Nitrate_25e-3_M.phsp'%(T))

        sut_G.append(sut_GT)
        ref_G.append(ref_GT)

        for i in sut_GT.keys():
            V1 = np.copy(sut_GT[i][66])
            V2 = np.copy(ref_GT[i][66])
            V1[0] = 1e12 / (float(T)*kobs)
            V2[0] = 1e12 / (float(T)*kobs)
            
            if i in sut_SG:
                sut_SG[i].append(V1)
                ref_SG[i].append(V2)

            else:
                sut_SG[i]   = [V1]
                ref_SG[i]   = [V2]

    for i in sut_SG.keys():
        sut_SG[i]   = np.array(sut_SG[i])
        ref_SG[i]   = np.array(ref_SG[i])

    bench1 = np.genfromtxt('analysis/benchmark/Pastina1999.csv')
    
    bench2 = np.genfromtxt('analysis/benchmark/Hiroki2002.csv')

    bench3 = np.genfromtxt("analysis/benchmark/Ramos-Mendez2021.csv")

    fig = plt.figure(figsize=(8,8))
    ax1 = plt.subplot2grid((2,1),(0,0))
    ax2 = plt.subplot2grid((2,1),(1,0))

    matplotlib.rcParams['font.family']     = "sans-serif"
    matplotlib.rcParams['font.sans-serif'] = "Helvetica"
    matplotlib.rcParams['figure.dpi']      = 200

    ax1.set_title(r"$H_{2}O_{2}$", fontsize=24)

    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.tick_params(axis='both', which='minor', labelsize=20)

    ax1.errorbar(sut_SG["H2O2^0"][:,0],  sut_SG["H2O2^0"][:,1], yerr=sut_SG["H2O2^0"][:,2],  color="r", marker='o', linewidth=2, label=args.sut_label)
    ax1.errorbar(ref_SG["H2O2^0"][:,0],  ref_SG["H2O2^0"][:,1], yerr=ref_SG["H2O2^0"][:,2],  color="g", marker='x', dashes=[8,5], linewidth=2, label=args.ref_label)
    ax1.scatter(1e12/bench1[:,0], bench1[:,1], label='Pastina(1999)')
    ax1.scatter(1e12/(bench2[:,0]*kobs), bench2[:,1], label='Hiroki(2002)')
    ax1.scatter(bench3[:,0], bench3[:,1], label='Ramos-Mendez(2021)')

    ax1.legend(loc=0)

    ax1.set_xlim(5e1,2E6)
#    ax1.set_ylim([2,4.5])

    ax1.set_xscale('log')

    ax1.set_xlabel("Time (ps)",fontsize=20)

    ax1.set_ylabel("GValue",fontsize=20)

    ax2.set_axis_off()
    CellText = [["",""],["",""],["",""]]
    CellText[0][0] = str(round(ref_T[0],2))+" +- "+str(round(ref_T[1],2))
    CellText[0][1] = str(round(sut_T[0],2))+" +- "+str(round(sut_T[1],2))
    CellText[1][0] = str(round(ref_T[2],2))+" +- "+str(round(ref_T[3],2))
    CellText[1][1] = str(round(sut_T[2],2))+" +- "+str(round(sut_T[3],2))
    CellText[2][0] = str(round(ref_T[4],2))+" +- "+str(round(ref_T[5],2))
    CellText[2][1] = str(round(sut_T[4],2))+" +- "+str(round(sut_T[5],2))
    Table = ax2.table(cellText   = CellText,\
                      rowLabels  = ["Real (s)","User (s)", "Sys (s)"],\
                      colLabels  = ["Reference","Under Test"],\
                      colWidths  = [0.5,0.5],\
                      rowColours = ["lightskyblue"]*10,\
                      colColours = ["lightskyblue"]*10,\
                      cellLoc    = 'center',\
                      loc        = 'center',\
                      fontsize   = 30)
    Table.auto_set_font_size(False)
    Table.set_fontsize(15)
    Table.scale(1,1.5)

    fig.tight_layout()
    fig.savefig(join(args.outdir,'TimeEvolution_IRTManager_comp.pdf'))

    plt.clf()
    plt.cla()
    plt.close()

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

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

    Results = {}
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
                Results[Mol][:,2] += Stdv1
        n += 1

    for i in Results.keys():
        Results[i] = np.array(Results[i])
        Results[i][:,1] /= n
        Results[i][:,2] /= n

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
    Temps = np.arange(9)*10 + 10

    sut_T = average_results_time(sut_dir + '/*/log.out')
    ref_T = average_results_time(ref_dir + '/*/log.out')

    sut_G   = []
    ref_G   = []
    bench_G = []

    sut_SG   = {}
    ref_SG   = {}
    bench_SG = {}


    for T in Temps:
        sut_GT   = average_results(sut_dir + '/*/electron_Gvalue_Corrected_%s.phsp'%(T))
        ref_GT   = average_results(ref_dir + '/*/electron_Gvalue_Corrected_%s.phsp'%(T))
        bench_GT = ReadGValues("./analysis/benchmark/electron_Gvalue_ByHand_%s.phsp"%(T))

        sut_G.append(sut_GT)
        ref_G.append(ref_GT)
        bench_G.append(bench_GT)

        for i in sut_GT.keys():
            V1 = np.copy(sut_GT[i][-1])
            V2 = np.copy(ref_GT[i][-1])
            V3 = np.copy(bench_GT[i][-1])
            V1[0] = T
            V2[0] = T
            V3[0] = T
            if i in sut_SG:
                sut_SG[i].append(V1)
                ref_SG[i].append(V2)
                bench_SG[i].append(V3)

            else:
                sut_SG[i]   = [V1]
                ref_SG[i]   = [V2]
                bench_SG[i] = [V3]

    for i in sut_SG.keys():
        sut_SG[i]   = np.array(sut_SG[i])
        ref_SG[i]   = np.array(ref_SG[i])
        bench_SG[i] = np.array(bench_SG[i])

    fig = plt.figure(figsize=(20,18))
    ax1 = plt.subplot2grid((3,3),(0,0))
    ax2 = plt.subplot2grid((3,3),(0,1))
    ax3 = plt.subplot2grid((3,3),(0,2))
    ax4 = plt.subplot2grid((3,3),(1,0))
    ax5 = plt.subplot2grid((3,3),(1,1))
    ax6 = plt.subplot2grid((3,3),(1,2))
    ax7 = plt.subplot2grid((3,3),(2,0))
    ax8 = plt.subplot2grid((3,3),(2,1),colspan=2)
    #ax9 = plt.subplot2grid((3,3),(2,2))

    matplotlib.rcParams['font.family']     = "sans-serif"
    matplotlib.rcParams['font.sans-serif'] = "Helvetica"
    matplotlib.rcParams['figure.dpi']      = 200

    ax1.set_title(r"$e_{aq}^{-}$", fontsize=24)
    ax2.set_title(r"$OH$",         fontsize=24)
    ax3.set_title(r"$H$",          fontsize=24)
    ax4.set_title(r"$H_{2}$",      fontsize=24)
    ax5.set_title(r"$H_{2}O_{2}$", fontsize=24)
    ax6.set_title(r"$H^{+}$",      fontsize=24)
    ax7.set_title(r"$OH^{-}$",     fontsize=24)

    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.tick_params(axis='both', which='minor', labelsize=20)

    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='minor', labelsize=20)

    ax3.tick_params(axis='both', which='major', labelsize=20)
    ax3.tick_params(axis='both', which='minor', labelsize=20)

    ax4.tick_params(axis='both', which='major', labelsize=20)
    ax4.tick_params(axis='both', which='minor', labelsize=20)

    ax5.tick_params(axis='both', which='major', labelsize=20)
    ax5.tick_params(axis='both', which='minor', labelsize=20)

    ax6.tick_params(axis='both', which='major', labelsize=20)
    ax6.tick_params(axis='both', which='minor', labelsize=20)

    ax7.tick_params(axis='both', which='major', labelsize=20)
    ax7.tick_params(axis='both', which='minor', labelsize=20)

    ax1.errorbar(bench_G[0]["e_aq^-1"][:,0][::10],bench_G[0]["e_aq^-1"][:,1][::10],yerr=bench_G[0]["e_aq^-1"][:,2][::10], color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax2.errorbar(bench_G[0]["OH^0"][:,0][::10],   bench_G[0]["OH^0"][:,1][::10],   yerr=bench_G[0]["OH^0"][:,2][::10],    color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax3.errorbar(bench_G[0]["H^0"][:,0][::10],    bench_G[0]["H^0"][:,1][::10],    yerr=bench_G[0]["H^0"][:,2][::10],     color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax4.errorbar(bench_G[0]["H_2^0"][:,0][::10],  bench_G[0]["H_2^0"][:,1][::10],  yerr=bench_G[0]["H_2^0"][:,2][::10],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax5.errorbar(bench_G[0]["H2O2^0"][:,0][::10], bench_G[0]["H2O2^0"][:,1][::10], yerr=bench_G[0]["H2O2^0"][:,2][::10],  color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax6.errorbar(bench_G[0]["H3O^1"][:,0][::10],  bench_G[0]["H3O^1"][:,1][::10],  yerr=bench_G[0]["H3O^1"][:,2][::10],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax7.errorbar(bench_G[0]["OH^-1"][:,0][::10],  bench_G[0]["OH^-1"][:,1][::10],  yerr=bench_G[0]["OH^-1"][:,2][::10],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)

    ax1.errorbar(bench_G[-1]["e_aq^-1"][:,0][::10],bench_G[-1]["e_aq^-1"][:,1][::10],yerr=bench_G[-1]["e_aq^-1"][:,2][::10], color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax2.errorbar(bench_G[-1]["OH^0"][:,0][::10],   bench_G[-1]["OH^0"][:,1][::10],   yerr=bench_G[-1]["OH^0"][:,2][::10],    color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax3.errorbar(bench_G[-1]["H^0"][:,0][::10],    bench_G[-1]["H^0"][:,1][::10],    yerr=bench_G[-1]["H^0"][:,2][::10],     color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax4.errorbar(bench_G[-1]["H_2^0"][:,0][::10],  bench_G[-1]["H_2^0"][:,1][::10],  yerr=bench_G[-1]["H_2^0"][:,2][::10],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax5.errorbar(bench_G[-1]["H2O2^0"][:,0][::10], bench_G[-1]["H2O2^0"][:,1][::10], yerr=bench_G[-1]["H2O2^0"][:,2][::10],  color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax6.errorbar(bench_G[-1]["H3O^1"][:,0][::10],  bench_G[-1]["H3O^1"][:,1][::10],  yerr=bench_G[-1]["H3O^1"][:,2][::10],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax7.errorbar(bench_G[-1]["OH^-1"][:,0][::10],  bench_G[-1]["OH^-1"][:,1][::10],  yerr=bench_G[-1]["OH^-1"][:,2][::10],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2, label='Ramos-Mendez (2022)')

    ax1.errorbar(sut_G[0]["e_aq^-1"][:,0],sut_G[0]["e_aq^-1"][:,1],yerr=sut_G[0]["e_aq^-1"][:,2], color="r", linewidth=2)
    ax2.errorbar(sut_G[0]["OH^0"][:,0],   sut_G[0]["OH^0"][:,1],   yerr=sut_G[0]["OH^0"][:,2],    color="r", linewidth=2)
    ax3.errorbar(sut_G[0]["H^0"][:,0],    sut_G[0]["H^0"][:,1],    yerr=sut_G[0]["H^0"][:,2],     color="r", linewidth=2)
    ax4.errorbar(sut_G[0]["H_2^0"][:,0],  sut_G[0]["H_2^0"][:,1],  yerr=sut_G[0]["H_2^0"][:,2],   color="r", linewidth=2)
    ax5.errorbar(sut_G[0]["H2O2^0"][:,0], sut_G[0]["H2O2^0"][:,1], yerr=sut_G[0]["H2O2^0"][:,2],  color="r", linewidth=2)
    ax6.errorbar(sut_G[0]["H3O^1"][:,0],  sut_G[0]["H3O^1"][:,1],  yerr=sut_G[0]["H3O^1"][:,2],   color="r", linewidth=2)
    ax7.errorbar(sut_G[0]["OH^-1"][:,0],  sut_G[0]["OH^-1"][:,1],  yerr=sut_G[0]["OH^-1"][:,2],   color="r", linewidth=2, label='{} 10$^\circ$C'.format(args.sut_label))

    ax1.errorbar(sut_G[-1]["e_aq^-1"][:,0],sut_G[-1]["e_aq^-1"][:,1],yerr=sut_G[-1]["e_aq^-1"][:,2], color="r", linewidth=2, linestyle='--')
    ax2.errorbar(sut_G[-1]["OH^0"][:,0],   sut_G[-1]["OH^0"][:,1],   yerr=sut_G[-1]["OH^0"][:,2],    color="r", linewidth=2, linestyle='--')
    ax3.errorbar(sut_G[-1]["H^0"][:,0],    sut_G[-1]["H^0"][:,1],    yerr=sut_G[-1]["H^0"][:,2],     color="r", linewidth=2, linestyle='--')
    ax4.errorbar(sut_G[-1]["H_2^0"][:,0],  sut_G[-1]["H_2^0"][:,1],  yerr=sut_G[-1]["H_2^0"][:,2],   color="r", linewidth=2, linestyle='--')
    ax5.errorbar(sut_G[-1]["H2O2^0"][:,0], sut_G[-1]["H2O2^0"][:,1], yerr=sut_G[-1]["H2O2^0"][:,2],  color="r", linewidth=2, linestyle='--')
    ax6.errorbar(sut_G[-1]["H3O^1"][:,0],  sut_G[-1]["H3O^1"][:,1],  yerr=sut_G[-1]["H3O^1"][:,2],   color="r", linewidth=2, linestyle='--')
    ax7.errorbar(sut_G[-1]["OH^-1"][:,0],  sut_G[-1]["OH^-1"][:,1],  yerr=sut_G[-1]["OH^-1"][:,2],   color="r", linewidth=2, linestyle='--', label='{} 90$^\circ$C'.format(args.sut_label))

    ax1.errorbar(ref_G[0]["e_aq^-1"][:,0],ref_G[0]["e_aq^-1"][:,1],yerr=ref_G[0]["e_aq^-1"][:,2], color="g", linewidth=2)
    ax2.errorbar(ref_G[0]["OH^0"][:,0],   ref_G[0]["OH^0"][:,1],   yerr=ref_G[0]["OH^0"][:,2],    color="g", linewidth=2)
    ax3.errorbar(ref_G[0]["H^0"][:,0],    ref_G[0]["H^0"][:,1],    yerr=ref_G[0]["H^0"][:,2],     color="g", linewidth=2)
    ax4.errorbar(ref_G[0]["H_2^0"][:,0],  ref_G[0]["H_2^0"][:,1],  yerr=ref_G[0]["H_2^0"][:,2],   color="g", linewidth=2)
    ax5.errorbar(ref_G[0]["H2O2^0"][:,0], ref_G[0]["H2O2^0"][:,1], yerr=ref_G[0]["H2O2^0"][:,2],  color="g", linewidth=2)
    ax6.errorbar(ref_G[0]["H3O^1"][:,0],  ref_G[0]["H3O^1"][:,1],  yerr=ref_G[0]["H3O^1"][:,2],   color="g", linewidth=2)
    ax7.errorbar(ref_G[0]["OH^-1"][:,0],  ref_G[0]["OH^-1"][:,1],  yerr=ref_G[0]["OH^-1"][:,2],   color="g", linewidth=2, label='{} 10$^\circ$C'.format(args.ref_label))

    ax1.errorbar(ref_G[-1]["e_aq^-1"][:,0],ref_G[-1]["e_aq^-1"][:,1],yerr=ref_G[-1]["e_aq^-1"][:,2], color="g", linewidth=2, linestyle='--')
    ax2.errorbar(ref_G[-1]["OH^0"][:,0],   ref_G[-1]["OH^0"][:,1],   yerr=ref_G[-1]["OH^0"][:,2],    color="g", linewidth=2, linestyle='--')
    ax3.errorbar(ref_G[-1]["H^0"][:,0],    ref_G[-1]["H^0"][:,1],    yerr=ref_G[-1]["H^0"][:,2],     color="g", linewidth=2, linestyle='--')
    ax4.errorbar(ref_G[-1]["H_2^0"][:,0],  ref_G[-1]["H_2^0"][:,1],  yerr=ref_G[-1]["H_2^0"][:,2],   color="g", linewidth=2, linestyle='--')
    ax5.errorbar(ref_G[-1]["H2O2^0"][:,0], ref_G[-1]["H2O2^0"][:,1], yerr=ref_G[-1]["H2O2^0"][:,2],  color="g", linewidth=2, linestyle='--')
    ax6.errorbar(ref_G[-1]["H3O^1"][:,0],  ref_G[-1]["H3O^1"][:,1],  yerr=ref_G[-1]["H3O^1"][:,2],   color="g", linewidth=2, linestyle='--')
    ax7.errorbar(ref_G[-1]["OH^-1"][:,0],  ref_G[-1]["OH^-1"][:,1],  yerr=ref_G[-1]["OH^-1"][:,2],   color="g", linewidth=2, linestyle='--', label='{} 90$^\circ$C'.format(args.ref_label))

    ax1.set_xlim([1,1E7])
    ax2.set_xlim([1,1E7])
    ax3.set_xlim([1,1E7])
    ax4.set_xlim([1,1E7])
    ax5.set_xlim([1,1E7])
    ax6.set_xlim([1,1E7])
    ax7.set_xlim([1,1E7])

    ax1.set_ylim([2,4.5])
    ax2.set_ylim([2,5.5])
    ax3.set_ylim([0.45,0.7])
    ax4.set_ylim([0.2,0.6])
    ax5.set_ylim([0,0.8])
    ax6.set_ylim([2.5,4.5])
    ax7.set_ylim([0,0.7])

    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax3.set_xscale("log")
    ax4.set_xscale("log")
    ax5.set_xscale("log")
    ax6.set_xscale("log")
    ax7.set_xscale("log")

    ax1.set_xlabel("Time (ps)",fontsize=20)
    ax2.set_xlabel("Time (ps)",fontsize=20)
    ax3.set_xlabel("Time (ps)",fontsize=20)
    ax4.set_xlabel("Time (ps)",fontsize=20)
    ax5.set_xlabel("Time (ps)",fontsize=20)
    ax6.set_xlabel("Time (ps)",fontsize=20)
    ax7.set_xlabel("Time (ps)",fontsize=20)

    ax1.set_ylabel("GValue",fontsize=20)
    ax2.set_ylabel("GValue",fontsize=20)
    ax3.set_ylabel("GValue",fontsize=20)
    ax4.set_ylabel("GValue",fontsize=20)
    ax5.set_ylabel("GValue",fontsize=20)
    ax6.set_ylabel("GValue",fontsize=20)
    ax7.set_ylabel("GValue",fontsize=20)

    ax7.legend(loc=2)

    ax1.grid(True,dashes=[5,5])
    ax2.grid(True,dashes=[5,5])
    ax3.grid(True,dashes=[5,5])
    ax4.grid(True,dashes=[5,5])
    ax5.grid(True,dashes=[5,5])
    ax6.grid(True,dashes=[5,5])
    ax7.grid(True,dashes=[5,5])

    ax8.set_axis_off()
    CellText = [["",""],["",""],["",""]]
    CellText[0][0] = str(round(ref_T[0],2))+" +- "+str(round(ref_T[1],2))
    CellText[0][1] = str(round(sut_T[0],2))+" +- "+str(round(sut_T[1],2))
    CellText[1][0] = str(round(ref_T[2],2))+" +- "+str(round(ref_T[3],2))
    CellText[1][1] = str(round(sut_T[2],2))+" +- "+str(round(sut_T[3],2))
    CellText[2][0] = str(round(ref_T[4],2))+" +- "+str(round(ref_T[5],2))
    CellText[2][1] = str(round(sut_T[4],2))+" +- "+str(round(sut_T[5],2))
    Table = ax8.table(cellText   = CellText,\
                      rowLabels  = ["Real (s)","User (s)", "Sys (s)"],\
                      colLabels  = ["Reference","Under Test"],\
                      colWidths  = [0.5,0.5],\
                      rowColours = ["lightskyblue"]*10,\
                      colColours = ["lightskyblue"]*10,\
                      cellLoc    = 'center',\
                      loc        = 'center',\
                      fontsize   = 30)
    Table.auto_set_font_size(False)
    Table.set_fontsize(20)
    Table.scale(1,3)

    #ax9.set_axis_off()

    fig.tight_layout()
    fig.savefig(join(args.outdir,'TimeEvolution_IRTManager_comp.pdf'))

    plt.clf()
    plt.cla()
    plt.close()

    fig = plt.figure(figsize=(20,18))
    ax1 = plt.subplot2grid((3,3),(0,0))
    ax2 = plt.subplot2grid((3,3),(0,1))
    ax3 = plt.subplot2grid((3,3),(0,2))
    ax4 = plt.subplot2grid((3,3),(1,0))
    ax5 = plt.subplot2grid((3,3),(1,1))
    ax6 = plt.subplot2grid((3,3),(1,2))
    ax7 = plt.subplot2grid((3,3),(2,0))
    ax8 = plt.subplot2grid((3,3),(2,1),colspan=2)
    #ax9 = plt.subplot2grid((3,3),(2,2))

    matplotlib.rcParams['font.family']     = "sans-serif"
    matplotlib.rcParams['font.sans-serif'] = "Helvetica"
    matplotlib.rcParams['figure.dpi']      = 200

    ax1.set_title(r"$e_{aq}^{-}$", fontsize=24)
    ax2.set_title(r"$OH$",         fontsize=24)
    ax3.set_title(r"$H$",          fontsize=24)
    ax4.set_title(r"$H_{2}$",      fontsize=24)
    ax5.set_title(r"$H_{2}O_{2}$", fontsize=24)
    ax6.set_title(r"$H^{+}$",      fontsize=24)
    ax7.set_title(r"$OH^{-}$",     fontsize=24)

    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.tick_params(axis='both', which='minor', labelsize=20)

    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='minor', labelsize=20)

    ax3.tick_params(axis='both', which='major', labelsize=20)
    ax3.tick_params(axis='both', which='minor', labelsize=20)

    ax4.tick_params(axis='both', which='major', labelsize=20)
    ax4.tick_params(axis='both', which='minor', labelsize=20)

    ax5.tick_params(axis='both', which='major', labelsize=20)
    ax5.tick_params(axis='both', which='minor', labelsize=20)

    ax6.tick_params(axis='both', which='major', labelsize=20)
    ax6.tick_params(axis='both', which='minor', labelsize=20)

    ax7.tick_params(axis='both', which='major', labelsize=20)
    ax7.tick_params(axis='both', which='minor', labelsize=20)

    ax1.errorbar(bench_SG["e_aq^-1"][:,0], bench_SG["e_aq^-1"][:,1],yerr=bench_SG["e_aq^-1"][:,2], color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax2.errorbar(bench_SG["OH^0"][:,0],    bench_SG["OH^0"][:,1],   yerr=bench_SG["OH^0"][:,2],    color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax3.errorbar(bench_SG["H^0"][:,0],     bench_SG["H^0"][:,1],    yerr=bench_SG["H^0"][:,2],     color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax4.errorbar(bench_SG["H_2^0"][:,0],   bench_SG["H_2^0"][:,1],  yerr=bench_SG["H_2^0"][:,2],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax5.errorbar(bench_SG["H2O2^0"][:,0],  bench_SG["H2O2^0"][:,1], yerr=bench_SG["H2O2^0"][:,2],  color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax6.errorbar(bench_SG["H3O^1"][:,0],   bench_SG["H3O^1"][:,1],  yerr=bench_SG["H3O^1"][:,2],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2)
    ax7.errorbar(bench_SG["OH^-1"][:,0],   bench_SG["OH^-1"][:,1],  yerr=bench_SG["OH^-1"][:,2],   color="b", marker="o", dashes=[0,1], markerfacecolor="none", markersize=15, linewidth=2, label='Ramos-Mendez (2022)')

    ax1.errorbar(sut_SG["e_aq^-1"][:,0], sut_SG["e_aq^-1"][:,1],yerr=sut_SG["e_aq^-1"][:,2], color="r", linewidth=2)
    ax2.errorbar(sut_SG["OH^0"][:,0],    sut_SG["OH^0"][:,1],   yerr=sut_SG["OH^0"][:,2],    color="r", linewidth=2)
    ax3.errorbar(sut_SG["H^0"][:,0],     sut_SG["H^0"][:,1],    yerr=sut_SG["H^0"][:,2],     color="r", linewidth=2)
    ax4.errorbar(sut_SG["H_2^0"][:,0],   sut_SG["H_2^0"][:,1],  yerr=sut_SG["H_2^0"][:,2],   color="r", linewidth=2)
    ax5.errorbar(sut_SG["H2O2^0"][:,0],  sut_SG["H2O2^0"][:,1], yerr=sut_SG["H2O2^0"][:,2],  color="r", linewidth=2)
    ax6.errorbar(sut_SG["H3O^1"][:,0],   sut_SG["H3O^1"][:,1],  yerr=sut_SG["H3O^1"][:,2],   color="r", linewidth=2)
    ax7.errorbar(sut_SG["OH^-1"][:,0],   sut_SG["OH^-1"][:,1],  yerr=sut_SG["OH^-1"][:,2],   color="r", linewidth=2, label=args.sut_label)

    ax1.errorbar(ref_SG["e_aq^-1"][:,0], ref_SG["e_aq^-1"][:,1],yerr=ref_SG["e_aq^-1"][:,2], color="g", dashes=[8,5], linewidth=2)
    ax2.errorbar(ref_SG["OH^0"][:,0],    ref_SG["OH^0"][:,1],   yerr=ref_SG["OH^0"][:,2],    color="g", dashes=[8,5], linewidth=2)
    ax3.errorbar(ref_SG["H^0"][:,0],     ref_SG["H^0"][:,1],    yerr=ref_SG["H^0"][:,2],     color="g", dashes=[8,5], linewidth=2)
    ax4.errorbar(ref_SG["H_2^0"][:,0],   ref_SG["H_2^0"][:,1],  yerr=ref_SG["H_2^0"][:,2],   color="g", dashes=[8,5], linewidth=2)
    ax5.errorbar(ref_SG["H2O2^0"][:,0],  ref_SG["H2O2^0"][:,1], yerr=ref_SG["H2O2^0"][:,2],  color="g", dashes=[8,5], linewidth=2)
    ax6.errorbar(ref_SG["H3O^1"][:,0],   ref_SG["H3O^1"][:,1],  yerr=ref_SG["H3O^1"][:,2],   color="g", dashes=[8,5], linewidth=2)
    ax7.errorbar(ref_SG["OH^-1"][:,0],   ref_SG["OH^-1"][:,1],  yerr=ref_SG["OH^-1"][:,2],   color="g", dashes=[8,5], linewidth=2, label=args.ref_label)

    ax1.set_xlim([0,100])
    ax2.set_xlim([0,100])
    ax3.set_xlim([0,100])
    ax4.set_xlim([0,100])
    ax5.set_xlim([0,100])
    ax6.set_xlim([0,100])
    ax7.set_xlim([0,100])

    ax1.set_ylim([2,4.5])
    ax2.set_ylim([2,5.5])
    ax3.set_ylim([0.45,0.7])
    ax4.set_ylim([0.2,0.6])
    ax5.set_ylim([0,0.8])
    ax6.set_ylim([2.5,4.5])
    ax7.set_ylim([0,0.7])

    ax1.set_xlabel("Temperature (C)",fontsize=20)
    ax2.set_xlabel("Temperature (C)",fontsize=20)
    ax3.set_xlabel("Temperature (C)",fontsize=20)
    ax4.set_xlabel("Temperature (C)",fontsize=20)
    ax5.set_xlabel("Temperature (C)",fontsize=20)
    ax6.set_xlabel("Temperature (C)",fontsize=20)
    ax7.set_xlabel("Temperature (C)",fontsize=20)

    ax1.set_ylabel("GValue",fontsize=20)
    ax2.set_ylabel("GValue",fontsize=20)
    ax3.set_ylabel("GValue",fontsize=20)
    ax4.set_ylabel("GValue",fontsize=20)
    ax5.set_ylabel("GValue",fontsize=20)
    ax6.set_ylabel("GValue",fontsize=20)
    ax7.set_ylabel("GValue",fontsize=20)

    ax7.legend(loc=0)

    ax1.grid(True,dashes=[5,5])
    ax2.grid(True,dashes=[5,5])
    ax3.grid(True,dashes=[5,5])
    ax4.grid(True,dashes=[5,5])
    ax5.grid(True,dashes=[5,5])
    ax6.grid(True,dashes=[5,5])
    ax7.grid(True,dashes=[5,5])

    ax8.set_axis_off()
    CellText = [["",""],["",""],["",""]]
    CellText[0][0] = str(round(ref_T[0],2))+" +- "+str(round(ref_T[1],2))
    CellText[0][1] = str(round(sut_T[0],2))+" +- "+str(round(sut_T[1],2))
    CellText[1][0] = str(round(ref_T[2],2))+" +- "+str(round(ref_T[3],2))
    CellText[1][1] = str(round(sut_T[2],2))+" +- "+str(round(sut_T[3],2))
    CellText[2][0] = str(round(ref_T[4],2))+" +- "+str(round(ref_T[5],2))
    CellText[2][1] = str(round(sut_T[4],2))+" +- "+str(round(sut_T[5],2))
    Table = ax8.table(cellText   = CellText,\
                      rowLabels  = ["Real (s)","User (s)", "Sys (s)"],\
                      colLabels  = ["Reference","Under Test"],\
                      colWidths  = [0.5,0.5],\
                      rowColours = ["lightskyblue"]*10,\
                      colColours = ["lightskyblue"]*10,\
                      cellLoc    = 'center',\
                      loc        = 'center',\
                      fontsize   = 30)
    Table.auto_set_font_size(False)
    Table.set_fontsize(20)
    Table.scale(1,3)
    
    #ax9.set_axis_off()

    fig.tight_layout()
    fig.savefig(join(args.outdir,'TemperatureEvolution_IRTManager_comp.pdf'))

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

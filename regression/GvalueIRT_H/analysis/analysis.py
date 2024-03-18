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
import copy

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

# Function to calculate the x-coordinate for the first dataset (Huerta benchmark)
def huerta_1(x, y):
    xf = []
    yf = []
    yerr = []

    for i in range(len(x)):
        xf = np.append(xf, 1e12 / (x[i] * 2.1e8 + 1e-3 * 1.4e6))
        yf = np.append(yf, y[i] - 0.460526) # H2 at 0 scav
        yerr = np.append(yerr, y[i]*0.05) # 5% experimental uncertainty?

    return xf, yf, yerr

# Function to calculate the x-coordinate for the second dataset (Huerta benchmark)
def huerta_2(x, y):
    xf = []
    yf = []
    yerr = []

    for i in range(len(x)):
        xf = np.append(xf, 1e12/x[i])
        yf = np.append(yf, y[i] - 0.460526) # H2 at 0 scav
        yerr = np.append(yerr, y[i]*0.05) # 5% experimental uncertainty?

    return xf, yf, yerr

####################################################
### Define Plot Function

def plot_results(sut_dir, ref_dir, args):
    HCO2m_concentrations = ["1e-2", "1e-1", "300e-3", "1e-0"]
    NO3m_concentrations = ["1e-3", "24e-3"]

    sut_T = average_results_time(sut_dir + '/*/log.out')
    ref_T = average_results_time(ref_dir + '/*/log.out')

    # H2 G-values for baseline case of {H_2}^0 (without formate) as defined in Huerta Parajon 2008
    sut_SG_H20 = {}
    sut_SG_H20_err = {}
    ref_SG_H20   = {}
    ref_SG_H20_err = {}

    # H2 G-values for 1mM Nitrate
    sut_SG_H2_1mM = {}
    ref_SG_H2_1mM = {}

    # H2 G-values for 24mM Nitrate
    sut_SG_H2_24mM = {}
    ref_SG_H2_24mM = {}
    
    sut_GT = average_results(sut_dir + '/*/Br_1e-3_M_NO3m_1e-3_M.phsp')
    ref_GT = average_results(ref_dir + '/*/Br_1e-3_M_NO3m_1e-3_M.phsp')

    V1 = np.copy(sut_GT['H_2^0'][99]) # taking 1e+9 [99] instead of 1e+6 [66], since H2 not stabilised at 1us
    V2 = np.copy(ref_GT['H_2^0'][99])

    # Getting escape yield of H2, which represents {H_2}^0 in Huerta Parajon 2008. Multiplied to match number of concentrations
    sut_SG_H20['H_2^0'] = [V1[1]]*len(HCO2m_concentrations)
    sut_SG_H20_err['H_2^0'] = [V1[2]]*len(HCO2m_concentrations)

    ref_SG_H20['H_2^0'] = [V2[1]]*len(HCO2m_concentrations)
    ref_SG_H20_err['H_2^0'] = [V2[2]]*len(HCO2m_concentrations)

    kobsHCO2m = 2.1e8
    kobsNO3m = 1.4e6

    # Now getting H2 for different concentrations of Nitrate
    for N in NO3m_concentrations: 
        for T in HCO2m_concentrations:
            sut_GT   = average_results(sut_dir + '/*/HCO2m_%s_M_Br_1e-3_M_NO3m_%s_M.phsp'%(T,N))
            ref_GT   = average_results(ref_dir + '/*/HCO2m_%s_M_Br_1e-3_M_NO3m_%s_M.phsp'%(T,N))

            for i in sut_GT.keys():
                V1 = np.copy(sut_GT[i][99]) # taking 1e+9 [99] instead of 1e+6 [66], since H2 not stabilised at 1us
                V2 = np.copy(ref_GT[i][99])
                V1[0] = 1e12 / (float(T)*kobsHCO2m + 1e-3*kobsNO3m)
                V2[0] = 1e12 / (float(T)*kobsHCO2m + 1e-3*kobsNO3m)
                
                if N == '1e-3':
                    if i in sut_SG_H2_1mM:
                        sut_SG_H2_1mM[i].append(V1)
                        ref_SG_H2_1mM[i].append(V2)

                    else:
                        sut_SG_H2_1mM[i]   = [V1]
                        ref_SG_H2_1mM[i]   = [V2]
                else:
                    if i in sut_SG_H2_24mM:
                        sut_SG_H2_24mM[i].append(V1)
                        ref_SG_H2_24mM[i].append(V2)

                    else:
                        sut_SG_H2_24mM[i]   = [V1]
                        ref_SG_H2_24mM[i]   = [V2]

    ##### Subtraction to obtain H yields #####
              
    # Topas under test
    time_sut_1mM = {}
    value_sut_1mM = {}
    error_sut_1mM = {}

    time_sut_24mM = {}
    value_sut_24mM = {}
    error_sut_24mM = {}

    for key in sut_SG_H20: # only interested in H2
        for i in range(len(sut_SG_H2_1mM[key])):
            if i == 0:
                time_sut_1mM[key] = [sut_SG_H2_1mM[key][i][0]]
                value_sut_1mM[key] = [sut_SG_H2_1mM[key][i][1] - sut_SG_H20[key][i]]
                error_sut_1mM[key] = [sut_SG_H2_1mM[key][i][2] + sut_SG_H20_err[key][i]]

                time_sut_24mM[key] = [sut_SG_H2_24mM[key][i][0]]
                value_sut_24mM[key] = [sut_SG_H2_24mM[key][i][1] - sut_SG_H20[key][i]]
                error_sut_24mM[key] = [sut_SG_H2_24mM[key][i][2] + sut_SG_H20_err[key][i]]

            else:
                time_sut_1mM[key].append(sut_SG_H2_1mM[key][i][0])
                value_sut_1mM[key].append(sut_SG_H2_1mM[key][i][1] - sut_SG_H20[key][i])
                error_sut_1mM[key].append(sut_SG_H2_1mM[key][i][2] + sut_SG_H20_err[key][i])

                time_sut_24mM[key].append(sut_SG_H2_24mM[key][i][0])
                value_sut_24mM[key].append(sut_SG_H2_24mM[key][i][1] - sut_SG_H20[key][i])
                error_sut_24mM[key].append(sut_SG_H2_24mM[key][i][2] + sut_SG_H20_err[key][i])


    print('Time, 1mM: ', time_sut_1mM['H_2^0'])
    print('G-value, 1mM: ', value_sut_1mM['H_2^0'])
    print('Error, 1mM: ', error_sut_1mM['H_2^0'],'\n')

    print('Time, 24mM: ', time_sut_24mM['H_2^0'])
    print('G-value, 24mM: ', value_sut_24mM['H_2^0'])
    print('Error, 24mM: ', error_sut_24mM['H_2^0'],'\n')

    # Topas reference
    time_ref_1mM = {}
    value_ref_1mM = {}
    error_ref_1mM = {}

    time_ref_24mM = {}
    value_ref_24mM = {}
    error_ref_24mM = {}

    for key in ref_SG_H20: # only interested in H2
        for i in range(len(ref_SG_H2_1mM[key])):
            if i == 0:
                time_ref_1mM[key] = [ref_SG_H2_1mM[key][i][0]]
                value_ref_1mM[key] = [ref_SG_H2_1mM[key][i][1] - ref_SG_H20[key][i]]
                error_ref_1mM[key] = [ref_SG_H2_1mM[key][i][2] + ref_SG_H20_err[key][i]]

                time_ref_24mM[key] = [ref_SG_H2_24mM[key][i][0]]
                value_ref_24mM[key] = [ref_SG_H2_24mM[key][i][1] - ref_SG_H20[key][i]]
                error_ref_24mM[key] = [ref_SG_H2_24mM[key][i][2] + ref_SG_H20_err[key][i]]

            else:
                time_ref_1mM[key].append(ref_SG_H2_1mM[key][i][0])
                value_ref_1mM[key].append(ref_SG_H2_1mM[key][i][1] - ref_SG_H20[key][i])
                error_ref_1mM[key].append(ref_SG_H2_1mM[key][i][2] + ref_SG_H20_err[key][i])

                time_ref_24mM[key].append(ref_SG_H2_24mM[key][i][0])
                value_ref_24mM[key].append(ref_SG_H2_24mM[key][i][1] - ref_SG_H20[key][i])
                error_ref_24mM[key].append(ref_SG_H2_24mM[key][i][2] + ref_SG_H20_err[key][i])

    # Huerta data
    #bench1 = np.genfromtxt('analysis/benchmark/HuertaParajon.csv')
    bench1 = np.loadtxt("./analysis/benchmark/HuertaParajon_cleaned.csv")
    tr = list(zip(*bench1))
    x1 = np.array(tr[0][0:5])
    y1 = np.array(tr[1][0:5])
    xf1, yf1, yerr1 = huerta_1(x1,y1)

    x2 = np.array(tr[0][5:9])
    y2 = np.array(tr[1][5:9])
    xf2, yf2, yerr2 = huerta_2(x2,y2)

    fig = plt.figure(figsize=(8,8))
    ax1 = plt.subplot2grid((2,1),(0,0))
    ax2 = plt.subplot2grid((2,1),(1,0))

    matplotlib.rcParams['font.family']     = "sans-serif"
    matplotlib.rcParams['font.sans-serif'] = "Helvetica"
    matplotlib.rcParams['figure.dpi']      = 200

    ax1.set_title("H", fontsize=24)

    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.tick_params(axis='both', which='minor', labelsize=20)

    # print(final_sut["H_2^0"])
    ax1.errorbar(time_sut_1mM["H_2^0"],  value_sut_1mM["H_2^0"], yerr=error_sut_1mM["H_2^0"],  color="r", marker='o', linewidth=2, label='{}, 1mM'.format(args.sut_label))
    ax1.errorbar(time_ref_1mM["H_2^0"],  value_ref_1mM["H_2^0"], yerr=error_ref_1mM["H_2^0"],  color="k", marker='x', dashes=[8,5], linewidth=2, label='{}, 1mM'.format(args.ref_label))
    
    ax1.errorbar(time_sut_24mM["H_2^0"],  value_sut_24mM["H_2^0"], yerr=error_sut_24mM["H_2^0"],  color="b", marker='o', linewidth=2, label='{}, 24mM'.format(args.sut_label))
    ax1.errorbar(time_ref_24mM["H_2^0"],  value_ref_24mM["H_2^0"], yerr=error_ref_24mM["H_2^0"],  color="purple", marker='x', dashes=[8,5], linewidth=2, label='{}, 24mM'.format(args.ref_label))

    ax1.errorbar(xf1, yf1, yerr=yerr1, linestyle='', marker='o', markersize=5, capsize=2, color='green', label='Huerta Parajon 2008')
    ax1.errorbar(xf2, yf2, yerr=yerr2, linestyle='', marker='o', markersize=5, capsize=2, color='green')

    #ax1.scatter(1e12/bench1[:,0], bench1[:,1], label='Pastina(1999)')
    #ax1.scatter(bench2[:,0], bench3[:,1], label='Ramos-Mendez(2021)')

    ax1.legend(loc=0)

    ax1.set_xlim(5e1,1E7)
#   ax1.set_ylim([2,4.5])

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

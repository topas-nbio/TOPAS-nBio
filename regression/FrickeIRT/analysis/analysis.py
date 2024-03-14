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

####################################################
### Define Read Functions

def GetGValue(name, accepted):
    gvalue = {}
    gvalue2 = {}
    time = {}
    iFile = open(name, 'r')
    for line in iFile:
        line = line.split()
        mol = line[3]
        if mol in accepted:
            if not mol in gvalue:
                gvalue[mol] = []		
                gvalue2[mol] = []		
                time[mol] = []		
            else:
                gvalue2[mol].append(float(line[0])**2)
                gvalue[mol].append(float(line[0]))
                time[mol].append(float(line[2]))

    return (gvalue, gvalue2, time)


def Glob(match_path):
    os.system('ls ' + match_path + ' > list')
    listFiles = open('list','r')
    fnames = []
    for name in listFiles:
        name = name.split()[0] 
        fnames.append(name)
    os.system('rm list')
    return fnames

def GetFromTxt(name, bins):
    iFile = open(name,'r')
    clusterSizePerEvent = {}
    for line in iFile:
        line = line.split()
        if not 'Ionisation' in line[5]:
            continue
        event = int(line[4])
        if not event in clusterSizePerEvent:
            clusterSizePerEvent[event] = 1
        else:
            clusterSizePerEvent[event] += 1
  
    clusterSizes = np.asarray(clusterSizePerEvent.values())
    
    hist, bins = np.histogram(clusterSizes, np.linspace(0, bins, bins+1), density=True)
 
    return (bins[:-1], hist) 
 
def average_results(match_path, molecules):
    fnames = Glob(match_path) 
    if len(fnames) == 0:
        return None

    gvalue, gvalue2, time = GetGValue(fnames[0], molecules) 
    n = 1
    for i in range(1,len(fnames)):
        agvalue, agvalue2, atime = GetGValue(fnames[i], molecules)
        for j in agvalue.items():
            name = j[0]
            value = j[1]
            for j in range(len(value)):
                gvalue[name][j] += agvalue[name][j]
                gvalue2[name][j] += agvalue2[name][j]

        n += 1

    red = {}
    ox = {}
    for name in molecules:
        for i in range(len(gvalue[name])):
            g = gvalue[name][i]/n
            g2 = gvalue2[name][i]
            if n > 1:
                g2 = np.sqrt( (1.0/(n-1)) * (g2/n - g*g))
            else:
                g2 = 0.0
            gvalue[name][i] = g#100*g
            gvalue2[name][i] = g2#100*g2

    return (time, gvalue, gvalue2)

def average_results_time(match_path):
    fnames = Glob(match_path)
    if len(fnames) == 0:
        return None

    execution, executionE = 0.0, 0.0
    initialization, initializationE = 0.0, 0.0
    finalization, finalizationE = 0.0, 0.0
    n = 0.0
    for f in fnames:
        t = float(subprocess.check_output("grep Execution " + f + " | awk '{print $2}'", shell=True)[5:-2])
        execution += t
        executionE += t*t
        t = float(subprocess.check_output("grep Initialization " + f + " | awk '{print $2}'", shell=True)[5:-2])
        initialization += t
        initializationE += t*t
        t = float(subprocess.check_output("grep Finalization " + f + " | awk '{print $2}'", shell=True)[5:-2])
        finalization += t
        finalizationE += t*t
        n += 1

    execution /= n
    initialization /= n
    finalization /= n
    if n > 1:
        executionE = np.sqrt(n/(n-1) * (np.abs(executionE/n - execution**2)))
        initializationE = np.sqrt(n/(n-1) * (np.abs(initializationE/n - initialization**2)))
        finalizationE = np.sqrt(n/(n-1) * (np.abs(finalizationE/n - finalization**2)))
    else:
        executionE, initializationE, finalizationE = 0., 0., 0.

    return (initialization, initializationE, execution, executionE, finalization, finalizationE)

def GetRatio(Ref,Sut):
    Ratio = []
    for i in range(len(Ref)):
        Ref_Value = Ref[i]
        Sut_Value = Sut[i]

        if Ref_Value == Sut_Value:
            Ratio.append(1)

        elif Sut_Value == 0:
            Ratio.append(1)

        elif Ref_Value == 0:
            Ratio.append(1)

        else:
            Ratio.append(Sut_Value/Ref_Value)

    return np.array(Ratio)

####################################################
### Define Plot Function

def plot_results(sut_dir, ref_dir, args):
    tag = split(sut_dir.rstrip(os.sep))[-1]
                
    fig = plt.figure()

    molecules = ['Fe3']

    name = 'Fe3_Gvalue.phsp'
    namePrefix = sut_dir + '/*/' + name 
    sut = average_results(namePrefix,molecules)
    namePrefix = ref_dir + '/*/' + name 
    ref = average_results(namePrefix,molecules)

    bench = np.genfromtxt("./analysis/benchmark/Plante2011.txt")
    
    grid = plt.GridSpec(2,2)
    plt.subplot(grid[1])
    plt.errorbar(sut[0]['Fe3'], sut[1]['Fe3'], yerr=sut[2]['Fe3'], fmt='b--',label='TOPAS Test')
    plt.errorbar(ref[0]['Fe3'], ref[1]['Fe3'], yerr=ref[2]['Fe3'], fmt='r:',label='TOPAS Ref')
    plt.errorbar(5E13, 15.6, yerr=0.2, fmt='o', markerfacecolor="none", label='ICRU report 34')
    plt.plot(bench[:,0]*1E12,bench[:,1], label='Plante 2011')
    plt.xscale("log")
    plt.xlim(1,1E14)
    plt.yticks((0 ,5, 10, 15))
    plt.grid(True,dashes=[5,5])
    plt.ylim(0,16)
    plt.legend(loc=2, fontsize=8)

    plt.ylabel(r'G($Fe^{3+}$ / 100 eV)')
    plt.xlabel('Time (ps)')

    plt.tight_layout() 

    plt.subplot(grid[0])
    sut_time = average_results_time(sut_dir + '/*/' + 'log.out')
    ref_time = average_results_time(ref_dir + '/*/' + 'log.out')

    plt.axis('off')
    table = plt.table(cellText=[['%1.3f +/- %1.3f' % (sut_time[0],sut_time[1]), '%1.3f +/- %1.3f' % (ref_time[0],ref_time[1])],\
                                ['%1.3f +/- %1.3f' % (sut_time[2],sut_time[3]), '%1.3f +/- %1.3f' % (ref_time[2],ref_time[3])],\
                                ['%1.3f +/- %1.3f' % (sut_time[4],sut_time[5]), '%1.3f +/- %1.3f' % (ref_time[4],ref_time[5])],\
                                ['%1.3f +/- %1.3f' % (sut[1]["Fe3"][-1],sut[2]["Fe3"][-1]), '%1.3f +/- %1.3f' % (ref[1]["Fe3"][-1],ref[2]["Fe3"][-1])]],
                      rowLabels=('Init. (s)','Exec. (s)','Final. (s)', 'Value (/100eV)'),\
                      colLabels=(args.sut_label,args.ref_label),\
                      loc='center')

    plt.subplot(grid[3])
    plt.plot(sut[0]['Fe3'], GetRatio(ref[1]["Fe3"],sut[1]["Fe3"]))
    plt.xlim(1,1E14)
    plt.ylim(.9,1.1)
    plt.xlabel('Time (ps)')
    plt.ylabel('Test/Ref ratio')
    plt.xscale("log")

    plt.tight_layout() 
    plt.savefig(join(args.outdir, 'Gvalue.eps'), bbox_inches='tight')
    plt.savefig(join(args.outdir, 'Gvalue.pdf'), bbox_inches='tight')
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

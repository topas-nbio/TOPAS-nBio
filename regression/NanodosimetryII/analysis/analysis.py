#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os, subprocess
from os.path import isdir, join, split
from glob import glob
from functools import reduce

from pylab import rcParams
rcParams['legend.numpoints'] = 1

def Glob(match_path):
    os.system('ls ' + match_path + ' > list')
    listFiles = open('list','r')
    fnames = []
    for name in listFiles:
        name = name.split()[0] 
        fnames.append(name)
    os.system('rm list')
    return fnames

def GetFromBinary(name, bins, histories):
    clusterSizePerEvent = {}
    data = np.fromfile(name, 'i,f,f,f,i')
    for event in data['f4']:
        if not event in clusterSizePerEvent:
            clusterSizePerEvent[event] = 1
        else:
            clusterSizePerEvent[event] += 1
  
    clusterSizes = np.asarray(list(clusterSizePerEvent.values()))
    zeroClusterSize = histories - len(clusterSizes) 

    hist, bins = np.histogram(clusterSizes, np.linspace(0, bins, bins+1)) #, density=True)
    hist = np.asarray(hist,dtype='d')
    hist[0] = zeroClusterSize
    hist /= np.trapz(hist, bins[:-1])
    return (bins[:-1], hist) 

def GetM1andF2(bins, hist, yerror):
    M1, F2 = 0.0, 0.0
    eM1, eF2 = 0.0, 0.0
    for i in range(len(bins)):
        M1 += bins[i] * hist[i]
        eM1 += (bins[i] * yerror[i])**2
        if bins[i] >= 2:
            F2 += hist[i]
            eF2 += yerror[i]**2
    
    eM1, eF2 = np.sqrt(eM1), np.sqrt(eF2)
    return (M1, eM1, F2, eF2) 
 
def average_results(match_path,histories):
    fnames = Glob(match_path) 
    if len(fnames) == 0:
        return None

    maxClusterSize = 100 
    hClusterSize = np.zeros(maxClusterSize)
    hClusterSizeStdv = np.zeros(maxClusterSize)
    n = 0.0
    data = 0.0
    for f in fnames:
        data = GetFromBinary(f, maxClusterSize, histories)
        hClusterSize += data[1]
        hClusterSizeStdv += data[1]**2
        n += 1.0
    
    hClusterSize /= n

    if n > 1:
        hClusterSizeStdv = np.sqrt(n/(n-1.0) * (hClusterSizeStdv/n - hClusterSize**2))
    else:
        hClusterSizeStdv = np.zeros(maxClusterSize)

    bins = data[0]
    M1, eM2, F2, eF2 = GetM1andF2(bins, hClusterSize, hClusterSizeStdv)
    return (M1, eM2, F2, eF2) 

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

def plot_results(sut_dir, ref_dir, args):
    tag = split(sut_dir.rstrip(os.sep))[-1]
                
    fig = plt.figure(figsize=(5,6))
    grid = plt.GridSpec(3,1)

    sutM1F2 = np.zeros([40,4])
    refM1F2 = np.zeros([40,4])
    h = [100, 1000, 1000, 1000, 1000, 1000,\
         1000, 1000, 1000, 10000, 10000, 10000,\
         10000, 10000, 10000, 10000, 10000, 10000,\
         10000, 10000, 100, 1000, 1000, 1000, 1000, 1000,\
         1000, 1000, 1000, 10000, 10000, 10000,\
         10000, 10000, 10000, 10000, 10000, 10000,\
         10000, 10000]

    for i in range(40):
        if i < 10:
            name = 'ID_Run_000' + str(i) + '.phsp'
        else:
            name = 'ID_Run_00' + str(i) + '.phsp'

        namePrefix = sut_dir + '/*/' + name 
        sut = average_results(namePrefix, 10*h[i])

        namePrefix = ref_dir + '/*/' + name 
        ref = average_results(namePrefix, 10*h[i])

        sutM1F2[i,:] = np.array([sut[0], sut[1], sut[2], sut[3]])
        refM1F2[i,:] = np.array([ref[0], ref[1], ref[2], ref[3]])

    plt.subplot(grid[0:2])

    plt.errorbar(sutM1F2[:,0], sutM1F2[:,2], yerr=sutM1F2[:,3], ls='none',mfc='none',\
                         color='r',marker='o',label=args.sut_label,linewidth=2.0)

    plt.errorbar(refM1F2[:,0], refM1F2[:,2], yerr=refM1F2[:,3], ls='none',mfc='none',\
                          color='b',marker='s',label=args.ref_label,linewidth=2.0)

    electron = np.genfromtxt('analysis/benchmark/Conte_2018_F2_electron.csv',delimiter=',')
    plt.scatter(electron[:,0], electron[:,1], facecolors='none', edgecolors='y', marker='D',label='e-') 
    alpha = np.genfromtxt('analysis/benchmark/Conte_2018_F2_alpha.csv',delimiter=',')
    plt.scatter(alpha[:,0], alpha[:,1], facecolors='none', edgecolors='y', marker='v',label='alpha') 
    proton = np.genfromtxt('analysis/benchmark/Conte_2018_F2_proton.csv',delimiter=',')
    plt.scatter(proton[:,0], proton[:,1], facecolors='none', edgecolors='y', marker='o',label='proton')
    carbon = np.genfromtxt('analysis/benchmark/Conte_2018_F2_carbon.csv',delimiter=',')
    plt.scatter(carbon[:,0], carbon[:,1], facecolors='none', edgecolors='y', marker='s',label='carbon ion')

    plt.xlim(1e-2,15)
    plt.ylim(1e-3, 1.15)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('$F_{2}$')
    plt.xlabel('$M_{1}$')

    plt.legend(loc=0, borderaxespad=0., fontsize=10)
    plt.tight_layout() 

    #fig = plt.figure()
    plt.subplot(grid[-1])
    sut_time = average_results_time(sut_dir + '/*/' + 'log.out')
    ref_time = average_results_time(ref_dir + '/*/' + 'log.out')

    plt.axis('off')
    table = plt.table(cellText=[['%1.1f +/- %1.1f' % (sut_time[0],sut_time[1]), '%1.1f +/- %1.1f' % (ref_time[0],ref_time[1])],\
                                ['%1.1f +/- %1.1f' % (sut_time[2],sut_time[3]), '%1.1f +/- %1.1f' % (ref_time[2],ref_time[3])],\
                                ['%1.1f +/- %1.1f' % (sut_time[4],sut_time[5]), '%1.1f +/- %1.1f' % (ref_time[4],ref_time[5])]],\
                      rowLabels=('Init.','Exec.','Final.'),\
                      colLabels=(args.sut_label+' (s)',args.ref_label+' (s)'),\
                      loc='center'\
                      )

    plt.tight_layout() 
    plt.savefig(join(args.outdir, 'NanodosimetryII.eps'), bbox_inches='tight')
    plt.savefig(join(args.outdir, 'NanodosimetryII.pdf'), bbox_inches='tight')
    
if __name__ == '__main__':
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



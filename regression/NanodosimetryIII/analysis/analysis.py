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

def GetFromASCII(name, bins):
    clusterSize, ionClusterFreq, exClusterFreq, totalClusterFreq = np.loadtxt(name, unpack = True, dtype = int)

    ionClusterProb = np.zeros(bins)
    sum_freq = sum(ionClusterFreq)
    for i in range(len(ionClusterFreq)):
        ionClusterProb[i] = ionClusterFreq[i]/sum_freq

    append_zeros = np.zeros(bins - len(clusterSize))
    fclusterSize = np.append(clusterSize, append_zeros)

    return (fclusterSize, ionClusterProb)

# def GetFromBinary(name, bins):
#     clusterSizePerEvent = {}
#     data = np.fromfile(name, 'i,f,f,f,i,i')
#     data = data[np.where(data['f5'] > 0)]
#     for event, copy in zip(data['f4'],data['f5']):
#         if not event in clusterSizePerEvent:
#             clusterSizePerEvent[event] = {}
#         else:
#             if not copy in clusterSizePerEvent[event]:
#                 clusterSizePerEvent[event][copy] = 1
#             else:
#                 clusterSizePerEvent[event][copy] += 1

#     clusterSizes = []
#     for event in clusterSizePerEvent:
#         clusterSizes += clusterSizePerEvent[event].values()

#     hist, bins = np.histogram(clusterSizes, np.linspace(0, bins, bins+1), density=True)

#     return (bins[:-1], hist) 

def average_results(match_path):
    fnames = Glob(match_path) 
    if len(fnames) == 0:
        return None

    maxClusterSize = 350
    clusterFreq = np.zeros(maxClusterSize)
    clusterFreqStdv = np.zeros(maxClusterSize)
    n = 0.0
    data = 0.0
    for f in fnames:
        data = GetFromASCII(f, maxClusterSize)
        clusterFreq += data[1]
        clusterFreqStdv += data[1]**2
        n += 1.0
    
    clusterFreq /= n

    if n > 1:
        clusterFreqStdv = np.sqrt(n/(n-1.0) * (clusterFreqStdv/n - clusterFreq**2))
    else:
        clusterFreqStdv = np.zeros(maxClusterSize)

    bins = data[0]
    return (bins, clusterFreq, clusterFreqStdv)

# def average_results(match_path, normalize):
#     fnames = Glob(match_path) 
#     if len(fnames) == 0:
#         return None

#     maxClusterSize = 120
#     hClusterSize = np.zeros(maxClusterSize)
#     hClusterSizeStdv = np.zeros(maxClusterSize)
#     n = 0.0
#     data = 0.0
#     for f in fnames:
#         data = GetFromBinary(f, maxClusterSize)
#         hClusterSize += data[1]
#         hClusterSizeStdv += data[1]**2
#         n += 1.0
    
#     hClusterSize /= n

#     if n > 1:
#         hClusterSizeStdv = np.sqrt(n/(n-1.0) * (hClusterSizeStdv/n - hClusterSize**2))
#     else:
#         hClusterSizeStdv = np.zeros(maxClusterSize)

#     bins = data[0]
#     return (bins, hClusterSize, hClusterSizeStdv) 

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

def GetBenchmark(fileName):
    data = np.genfromtxt(fileName)
    return (data[:,0], data[:,1]) 

def plot_results(sut_dir, ref_dir, args):
    tag = split(sut_dir.rstrip(os.sep))[-1]
    p = ['p','a','li','be','b','c','n','o']
    pname = ['proton','alpha','lithium','beryllium','boron','carbon','nitrogen','oxygen']

    xlim = [14,22,32,45,55,70,80,110]
    fig = plt.figure()
    grid = plt.GridSpec(3,3) 

    for i in range(8):
        plt.subplot(grid[i])
        name = 'ID_Run_000' + str(i)+ '.phsp' 
        namePrefix = sut_dir + '/*/' + name 
        sut = average_results(namePrefix)
        namePrefix = ref_dir + '/*/' + name 
        ref = average_results(namePrefix)
 
        plt.step(sut[0][1:], sut[1][1:], where='mid',label=args.sut_label) 
        plt.step(ref[0][1:], ref[1][1:], where='mid',label=args.ref_label) 

        bench = GetBenchmark('analysis/benchmark/'+p[i]+'_ID_1_MeV.csv')

        plt.title(pname[i],fontsize=8)

        plt.step(bench[0][1:], bench[1][1:],where='mid',color='r', label='Ref.')
        plt.xlim((0,xlim[i]))

        plt.ylim((1e-4,1.0))
        plt.yscale('log')

        if i % 3. == 0:
            plt.ylabel('Probability')


        if i >= 5 and i <= 7:
           plt.xlabel('Ionization cluster size')

        if i == 5: 
           plt.legend(loc=0, borderaxespad=0., fontsize=7)

        plt.tight_layout() 
    
    plt.subplot(grid[-1]) 
    sut_time = average_results_time(sut_dir + '/*/' + 'log.out')
    ref_time = average_results_time(ref_dir + '/*/' + 'log.out')

    plt.axis('off')
    table = plt.table(cellText=[['%1.1f +/- %1.1f' % (sut_time[0],sut_time[1]), '%1.1f +/- %1.1f' % (ref_time[0],ref_time[1])],\
                                ['%1.1f +/- %1.1f' % (sut_time[2],sut_time[3]), '%1.1f +/- %1.1f' % (ref_time[2],ref_time[3])],\
                                ['%1.1f +/- %1.1f' % (sut_time[4],sut_time[5]), '%1.1f +/- %1.1f' % (ref_time[4],ref_time[5])]],\
                      rowLabels=('Init.','Exec.','Final.'),\
                      colLabels=(args.sut_label+' (s)',args.ref_label+' (s)'),\
                      loc='center',fontsize=10\
                      )

    plt.tight_layout() 
    plt.savefig(join(args.outdir, 'NanodosimetryIII.eps'), bbox_inches='tight')
    plt.savefig(join(args.outdir, 'NanodosimetryIII.pdf'), bbox_inches='tight')

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



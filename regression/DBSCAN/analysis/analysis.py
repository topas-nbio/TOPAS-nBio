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

def average_results(match_path, normalize):
    fnames = Glob(match_path) 
    if len(fnames) == 0:
        return None

    maxClusterSize = 60 #100
    hmeanSSB, hstdvSSB = np.zeros(maxClusterSize+1), np.zeros(maxClusterSize+1)
    hmeanDSB, hstdvDSB = np.zeros(maxClusterSize+1), np.zeros(maxClusterSize+1)
    meanSSB, meanDSB = 0., 0.
    stdvSSB, stdvDSB = 0., 0.
    n = 0.0
    for f in fnames:
        data = np.genfromtxt(f)
        meanDSB += data[:,2].mean()
        meanSSB += data[:,1].mean()
        stdvDSB += data[:,2].mean()**2
        stdvSSB += data[:,1].mean()**2

        maxSSB = int(data[:,1].max())
        ssb, bins = np.histogram(data[:,1], np.linspace(0,maxSSB,maxSSB+1), density=True)

        hmeanSSB[bins[:-1].astype('i')] += ssb 
        hstdvSSB[bins[:-1].astype('i')] += ssb**2

        maxDSB = int(data[:,2].max())
        dsb, bins = np.histogram(data[:,1], np.linspace(0,maxDSB,maxDSB+1), density=True)
        hmeanDSB[bins[:-1].astype('i')] += dsb 
        hstdvDSB[bins[:-1].astype('i')] += dsb**2 
        n += 1.0
    
    hmeanSSB /= n
    hmeanDSB /= n
    meanSSB /= n
    meanDSB /= n
    stdvSSB = np.sqrt(n/(n-1) * (stdvSSB/n - meanSSB**2))
    stdvDSB = np.sqrt(n/(n-1) * (stdvDSB/n - meanDSB**2))

    if n > 1:
        hstdvSSB = np.sqrt(1.0/(n-1.0) * (hstdvSSB/n - hmeanSSB**2))
        hstdvDSB = np.sqrt(1.0/(n-1.0) * (hstdvDSB/n - hmeanDSB**2))
    else:
        hstdvSSB = np.zeros(hmeanSSB.shape)
        hstdvDSB = np.zeros(hmeanDSB.shape)

    bins = np.linspace(0,maxClusterSize,maxClusterSize+1)
   
    ratio = meanDSB/meanSSB 
    ratioStd = ratio * np.sqrt( (stdvSSB/meanSSB)**2 + (stdvDSB/meanDSB)**2)
    return (bins, hmeanSSB, hstdvSSB, hmeanDSB, hstdvDSB, meanSSB, stdvSSB, meanDSB, stdvDSB, ratio, ratioStd)

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
    energies = {'0000': '0.5', '0001': '1.1', '0002': '2.3',\
                '0003': '5', '0004': '10.8', '0005': '23.2',\
                '0006': '50'}
                
    fig = plt.figure()
    i_plot = 1
    sut_dbscan = {}
    ref_dbscan = {}
    for idRun in sorted(energies): 
        energy = energies[idRun]
        name = 'DBSCAN_Run_'+idRun+'.phsp ' 

        namePrefix = sut_dir + '/*/' + name 
        sut_dbscan[energy] = average_results(namePrefix, True)

        namePrefix = ref_dir + '/*/' + name 
        ref_dbscan[energy] = average_results(namePrefix, True)

        plt.subplot(3,3,i_plot) 
        plt.title(energy + ' MeV')
        plt.yscale('log')
        plt.plot(sut_dbscan[energy][0], sut_dbscan[energy][1],\
                          'r-',label=args.sut_label,linewidth=2.0)

        plt.plot(ref_dbscan[energy][0], ref_dbscan[energy][1],\
                          'b-',label=args.ref_label,linewidth=2.0)

        if i_plot == 4:	
            plt.ylabel('Frequency')
        if i_plot > 4:
            plt.xlabel('Number of SSB')

        i_plot += 1

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.subplots_adjust(wspace=.5, hspace=.5) 
    plt.savefig(join(args.outdir, 'DBSCAN1.eps'), bbox_inches='tight')
    plt.savefig(join(args.outdir, 'DBSCAN1.pdf'), bbox_inches='tight')

    fig = plt.figure(figsize=(5,6))
    grid = plt.GridSpec(3,1)
    plt.subplot(grid[:2]) #2,1,1)
    plt.xscale('log')
    plt.xlim((0.4,60))
    plt.xlabel('Energy (MeV)')
    plt.ylabel('DSB/SSB')

    partrac = np.genfromtxt('analysis/benchmark/PARTRAC.csv') 
    plt.plot(partrac[:,0]*1e-6, partrac[:,1], '-', label='PARTRAC')

    Francis = np.genfromtxt('analysis/benchmark/Francis2011.csv') 
    plt.plot(Francis[:,0]*1e-6, Francis[:,1], 'o', label='Francis 2011')

    Leloup = np.genfromtxt('analysis/benchmark/Leloup.csv') 
    plt.errorbar(Leloup[:,0]*1e-6, Leloup[:,1], yerr=(Leloup[:,2]-Leloup[:,1]), marker='v', color='g', ls='none', label='Leloup exp')

    for energy, value in sut_dbscan.items():
        if '0.5' in energy:
            plt.errorbar(float(energy), value[9], yerr=value[10], color='r', marker='o', mfc='none', label=args.sut_label)
        else:
            plt.errorbar(float(energy), value[9], yerr=value[10], color='r', marker='o', mfc='none', )

    for energy, value in ref_dbscan.items():
        if '0.5' in energy:
            plt.errorbar(float(energy), value[9], yerr=value[10], color='b', marker='s', mfc='none', label=args.ref_label)
        else:
            plt.errorbar(float(energy), value[9], yerr=value[10], color='b', marker='s', mfc='none', )

    # From literature
    plt.legend(loc=0, prop={'size': 8})

    plt.subplot(grid[-1]) #2,1,2)
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
    plt.savefig(join(args.outdir, 'DBSCAN2.eps'), bbox_inches='tight')
    plt.savefig(join(args.outdir, 'DBSCAN2.pdf'), bbox_inches='tight')
    
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



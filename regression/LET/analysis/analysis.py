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

    n = 0.0
    mean, stdv = 0.0, 0.0
    for f in fnames:
        data = np.genfromtxt(f, delimiter=',')
        mean += data
        stdv += data**2 
        n += 1.0
    
    mean /= n
    if n > 1:
        stdv = np.sqrt(n/(n-1.0) * (stdv/n - mean**2))
    else:
        stdv = 0.0 

    return (mean, stdv)

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

    energies = ['0.5', '0.9', '1.6', '2.9', '5.3', '9.5', '17.1', '30.8', '55.5', '100']

    sut = np.zeros([len(energies),2])
    ref = np.zeros([len(energies),2])
    
    i = 0
    for energy in energies:
        name = 'energyDeposit_Run_000' + str(i) + '.csv'
        namePrefix = sut_dir + '/*/' + name 
        edep, edep_stdv = average_results(namePrefix, True)
        edep *= 1e3
        edep_stdv *= 1e3
        name = 'fluence_Run_000' + str(i) + '.csv'
        namePrefix = sut_dir + '/*/' + name 
        fluence, fluence_stdv = average_results(namePrefix, True)
        fluence *= 1e-6*421.16
        fluence_stdv *= 1e-6*421.16
        sut[i,0] = edep/fluence
        sut[i,1] = edep/fluence * np.sqrt((edep_stdv/edep)**2 + (fluence_stdv/fluence)**2)
         
        name = 'energyDeposit_Run_000' + str(i) + '.csv'
        namePrefix = ref_dir + '/*/' + name 
        edep, edep_stdv = average_results(namePrefix, True)
        edep *= 1e3
        edep_stdv *= 1e3
        name = 'fluence_Run_000' + str(i) + '.csv'
        namePrefix = ref_dir + '/*/' + name 
        fluence, fluence_stdv = average_results(namePrefix, True)
        fluence *= 1e-6*421.16
        fluence_stdv *= 1e-6*421.16
        ref[i,0] = edep/fluence
        ref[i,1] = edep/fluence * np.sqrt((edep_stdv/edep)**2 + (fluence_stdv/fluence)**2)

        i+=1

    grid = plt.GridSpec(3,1)
    plt.subplot(grid[0:2])
    energies = np.asarray(energies,'d')
    plt.errorbar(energies, sut[:,0], yerr=sut[:,1],\
                          color='r',marker='o',mfc='none',label=args.sut_label,linewidth=2.0)
    plt.errorbar(energies, ref[:,0], yerr=ref[:,1],\
                          color='b',marker='s',mfc='none',label=args.ref_label,linewidth=2.0)

    literature = np.genfromtxt('analysis/benchmark/PARTRAC.csv',delimiter=',')
    plt.plot(literature[:,0], literature[:,1], 'k-', label='PARTRAC')
 
    plt.xlim((0.4,200))
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('LET (keV/$\mu$m)')
    plt.xlabel('Proton energy (MeV)')

    plt.legend(loc=0, borderaxespad=0.)

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
                      loc='center'\
                      )

    plt.tight_layout() 
    plt.savefig(join(args.outdir, 'LET.eps'), bbox_inches='tight')
    plt.savefig(join(args.outdir, 'LET.pdf'), bbox_inches='tight')

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



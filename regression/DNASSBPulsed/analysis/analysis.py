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
        name = name.rsplit()[0] 
        if 'Strand' in name or 'header' in name:
            continue
        if os.path.getsize(name) <= 0:
            continue
        fnames.append(name)

    os.system('rm list')
    return fnames

def ReadGValues(File):
    f = open (File, 'r')

    Output = {}

    for z in f.readlines():
        A = z.split()

        if len(A) < 7:
            continue

        GVal = float(A[0])
        GErr = float(A[1])
        Mols = float(A[2])
        MolE = float(A[3])
        Time = float(A[5])
        Name = A[6]

        if Name in Output:
            Output[Name].append([Time,GVal,GErr,Mols])

        else:
            Output[Name] = [[Time, GVal, GErr, Mols]]

    for i in Output.keys():
        Output[i] = np.array(Output[i])

    return Output

def average_results(match_path):
    fnames = Glob(match_path) 
    if len(fnames) == 0:
        return None

    gvalue = 0
    g2value = {}
    i = 0
    for fileName in fnames: 
        avalue = ReadGValues(fileName)
        if i == 0:
            gvalue = avalue 
            for key in avalue:
                g2value[key] = avalue[key][:,1]**2
        else:
            for key in avalue:
                gvalue[key][:,3] += avalue[key][:,3]
                gvalue[key][:,1] += avalue[key][:,1]
                gvalue[key][:,2] += avalue[key][:,2]**2
                g2value[key] += avalue[key][:,1]**2
        i += 1

    for key in avalue:
        gvalue[key][:,3] /= i
        gvalue[key][:,1] /= i
        g2value[key] = np.sqrt(1.0/i * (g2value[key]/i - gvalue[key][:,1]**2))
        gvalue[key][:,2] = g2value[key] 

   
    ssb = 0.24 * gvalue['OHDeoxyriboseDamaged'][-1,1] +\
          0.08 * gvalue['HDeoxyriboseDamaged'][-1,1] 

    ssb_err = np.sqrt(\
          (0.24 * gvalue['OHDeoxyriboseDamaged'][-1,2])**2 +\
          (0.08 * gvalue['OHDeoxyriboseDamaged'][-1,2])**2)

    ssb *= 0.1036         # to umol/J
    ssb_err *= 0.1036     # to umol/J

    return (gvalue, ssb, ssb_err)

    
def average_results_time(match_path):
    fnames = Glob(match_path)
    if len(fnames) == 0:
        return None

    execution, executionE = 0.0, 0.0
    initialization, initializationE = 0.0, 0.0
    finalization, finalizationE = 0.0, 0.0
    n = 0.0
    for f in fnames:
        if os.path.getsize(f) <= 0:
            continue
        t = 0 
        try:
             subprocess.check_output("grep Execution " + f + " | awk '{print $2}' | sed 's/User=//g' | sed 's/s//g' | awk '{sum += $1} END { print sum/NR}'", shell=True)
             t = subprocess.check_output("grep Execution " + f + " | awk '{print $2}' | sed 's/User=//g' | sed 's/s//g' | awk '{sum += $1} END { print sum/NR}'", shell=True)
        except:
             continue

        if t:
            t = float(t) 
            execution += t
            executionE += t*t
            t = float(subprocess.check_output("grep Initialization " + f + " | awk '{print $2}' | sed 's/User=//g' | sed 's/s//g' | awk '{sum += $1} END { print sum/NR}'", shell=True))
            initialization += t
            initializationE += t*t
            t = float(subprocess.check_output("grep Finalization " + f + " | awk '{print $2}' | sed 's/User=//g' | sed 's/s//g' | awk '{sum += $1} END { print sum/NR}'", shell=True))
            finalization += t
            finalizationE += t*t
            n += 1

    execution /= n
    initialization /= n
    finalization /= n
    if n > 1:
        executionE = np.sqrt(1/n * (np.abs(executionE/n - execution**2)))
        initializationE = np.sqrt(1/n * (np.abs(initializationE/n - initialization**2)))
        finalizationE = np.sqrt(1/n * (np.abs(finalizationE/n - finalization**2)))
    else:
        executionE, initializationE, finalizationE = 0., 0., 0.

    return (initialization, initializationE, execution, executionE, finalization, finalizationE)

def plot_results(sut_dir, ref_dir, args):
    tag = split(sut_dir.rstrip(os.sep))[-1]
                
    fig = plt.figure(figsize=(5,6))

    concentrations = ['1e-4', '1e-3', '1e-2', '1e-1', '1e0']
    #concentrations = ['1e-3']

    sut = np.zeros([len(concentrations),2])
    ref = np.zeros([len(concentrations),2])
    
    i = 0
    for c in concentrations:
        name = 'IRT_GValue_' + c + '.phsp'
        namePrefix = sut_dir + '/*/' + name
        gvalue, ssb, ssb_err =  average_results(namePrefix)

        oFile = open(join(args.outdir, 'sut.csv'), 'w')
        print(gvalue.keys())
        for key in gvalue:
            for k in range(gvalue[key].shape[0]):
                oFile.write('%1.5e  %1.5e  %1.5e  %1.5e  %s\n' % (gvalue[key][k,0], gvalue[key][k,1], gvalue[key][k,2], gvalue[key][k,3], key))
            oFile.write('\n\n')
        oFile.close()

        sut[i,0] = ssb
        sut[i,1] = ssb_err

        namePrefix = ref_dir + '/*/' + name
        gvalue, ssb, ssb_err =  average_results(namePrefix)
        oFile = open(join(args.outdir, 'ref.csv'), 'w')
        for key in gvalue:
            for k in range(gvalue[key].shape[0]):
                oFile.write('%1.5e  %1.5e  %1.5e  %1.5e  %s\n' % (gvalue[key][k,0], gvalue[key][k,1], gvalue[key][k,2], gvalue[key][k,3],key))
            oFile.write('\n\n')
        oFile.close()

        ref[i,0] = ssb
        ref[i,1] = ssb_err
        i+=1

    grid = plt.GridSpec(3,1)
    plt.subplot(grid[0:2])
    concentrations = [float(c) for c in concentrations] 

    plt.errorbar(concentrations, sut[:,0], yerr=sut[:,1],\
                          color='r',marker='o',fmt='o',mfc='none',label=args.sut_label,linewidth=2.0)
    plt.errorbar(concentrations, ref[:,0], yerr=ref[:,1],\
                          color='b',marker='s',fmt='o',mfc='none',label=args.ref_label,linewidth=2.0)

    literatureA = np.genfromtxt('analysis/benchmark/pBR322_Milligan1993.dat')
    literatureB = np.genfromtxt('analysis/benchmark/pEC_Milligan1993.dat')
    literatureC = np.genfromtxt('analysis/benchmark/pUC18_Milligan1993.dat')

    plt.plot(literatureA[:,0], literatureA[:,1], 'o', label='pBR322, Milligan')
    plt.plot(literatureB[:,0], literatureB[:,1], 'o', label='pEC,   Milligan')
    plt.plot(literatureC[:,0], literatureC[:,1], 'o', label='pUC18, Milligan')
 
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('SSB (umol/J')
    plt.xlabel('DMSO [M]')

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
    plt.savefig(join(args.outdir, 'SSB_vs_DMSO.eps'), bbox_inches='tight')
    
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



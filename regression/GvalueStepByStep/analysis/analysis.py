#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os, subprocess
from os.path import isdir, join, split
from glob import glob
from functools import reduce
import matplotlib.lines as mlines

from pylab import rcParams
rcParams['legend.numpoints'] = 1

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
        for name, value in agvalue.items():
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
            gvalue[name][i] = g
            gvalue2[name][i] = g2

    OH = gvalue['OH^0']
    eOH = gvalue2['OH^0']
    eaq = gvalue['e_aq^-1']
    eeaq = gvalue2['e_aq^-1']
    H2O2 = gvalue['H2O2^0']
    eH2O2 = gvalue2['H2O2^0']
    H2 = gvalue['H_2^0']
    eH2 = gvalue2['H_2^0']
    H = gvalue['H^0']
    eH = gvalue2['H^0']
    red = np.zeros([len(OH),2])
    ox  = np.zeros([len(OH),2])
    ratio  = np.zeros([len(OH),2])
    for i in range(len(OH)):
         red[i,0] = eaq[i] + 2.0 * H2[i] + H[i]
         ox[i,0] = OH[i] + 2.0 * H2O2[i]
         red[i,1] = np.sqrt(eeaq[i]**2 + (2.0 * eH2[i])**2 + eH[i]**2)
         ox[i,1] = np.sqrt(eOH[i]**2 + (2.0 * eH2O2[i])**2)
    
    ratio[:,0] = ox[:,0]/red[:,0]
    ratio[:,1] = ratio[:,0] * np.sqrt( (red[:,1]/red[:,0])**2 + (ox[:,1]/ox[:,0])**2) 
    return (time, gvalue, gvalue2, red, ox, ratio)

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
                
    fig = plt.figure()

    molecules = ['OH^0','e_aq^-1','H3O^1','H2O2^0','H_2^0','H^0']
    moleculesName = ['$^{.}OH$','$e-_{aq}$','$H_{3}O^{+}$','$H_{2}O_{2}$','$H_{2}$','$H^{.}$']

    name = 'Gvalue.phsp'
    namePrefix = sut_dir + '/*/' + name 
    sut = average_results(namePrefix,molecules)
    namePrefix = ref_dir + '/*/' + name 
    ref = average_results(namePrefix,molecules)

    i = 1
    grid = plt.GridSpec(3,3)
    for molecule, name in zip(molecules, moleculesName): 
        plt.subplot(grid[i-1]) #3,3,i)
        plt.xlim((1,2e6))
        plt.xscale('log')
        plt.errorbar(sut[0][molecule], sut[1][molecule], yerr=sut[2][molecule],capsize=None,\
                          color='r',linestyle='-',label=args.sut_label,linewidth=0.8)
        #plt.step(sut[0][molecule], sut[1][molecule],\
                          #color='r')

        plt.errorbar(ref[0][molecule], ref[1][molecule], yerr=ref[2][molecule],capsize=None,\
                          color='b',linestyle='-',label=args.ref_label,linewidth=0.8)

        #plt.step(ref[0][molecule], ref[1][molecule],\
                          #color='b')


        # Benchmark 
        if i == 1:
            time = np.linspace(7e-3, 8.0, 1000)
            OHWang = 2.873 + 1.087 * np.exp(-time/0.1354) + 0.8786 * np.exp(-time/1.135)
            OHWang /= OHWang.max() / (5.0*1e-1/0.1036)
            plt.yticks(np.arange(2,6,1))
            plt.ylim((2.0,5.5))
            plt.plot(time*1e3, OHWang, linestyle='-.', linewidth=1.5)
            plt.errorbar(7, 5.0*1e-1/0.1035, yerr = 0.2*1e-1/0.1036, marker='o', color='orange', markersize=4)
            plt.errorbar(1e6, 2.5, yerr=0.12, marker='o', color='k', markersize=4)
            plt.errorbar(1e6, 2.65, yerr=0.30, marker='o', color='g', markersize=4)
            plt.errorbar(1e6, 2.87, yerr=0.14, marker='o', color='r', markersize=4)

            # TODO: Add the missing references for OH escape yield data points.
            refs = [mlines.Line2D([], [], color='orange', marker='o', linestyle='None', markersize=2, label='El Omar (2011)')]
            plt.legend(handles=refs, loc=0, fontsize=4)

        if i == 2:
            time = np.linspace(7e-3, 8.0, 1000)
            EAQWang = 2.885 + 1.153*np.exp(-time*0.2238) + 0.4572*np.exp(-time*3.369)
            EAQWang /= EAQWang.max() / (4.4*1e-1/0.1036)
            EAQBart = 2.5*(1 + 0.090*np.exp(-time/139) + 0.128*np.exp(-time/24.4) +\
                  0.255*np.exp(-time/3.51) + 0.118*np.exp(-time*0.480)) 

            plt.plot(time*1e3, EAQWang, linestyle='-.', linewidth=1.5)
            plt.plot(time*1e3, EAQBart, linestyle='-.', linewidth=1.5)
            plt.errorbar(7.0, 4.4*1e-1/0.1036, yerr=0.2*1e-1/0.1036, marker='o', color='g', markersize=4)
            plt.errorbar(1e6,2.67,yerr=0.15,marker='o', markersize=4)
            plt.errorbar(1e6,2.70,yerr=0.14,marker='o', markersize=4)
            plt.errorbar(1e6,2.65,yerr=0.05,marker='o', markersize=4)

            # TODO: Add the missing references for e-aq G-value data points.
            refs = [mlines.Line2D([], [], color='g', marker='o', linestyle='None', markersize=2, label='Sumiyoshi (2006)')]
            plt.legend(handles=refs, loc=0, fontsize=4)

 
        if i == 4:
            plt.yticks(np.arange(0,1,0.2))
            plt.errorbar(1e6,0.61,yerr=0.03,marker='o')
            plt.errorbar(1e6,0.69,yerr=0.03,marker='o')

        if i == 6:
            plt.yticks(np.arange(0.55,0.72,0.04))
  
        plt.ylabel('G(%s / 100 eV)' % name)
        plt.xlabel('Time (ps)')
        i += 1

    plt.subplot(grid[i-1]) 
    plt.yticks(np.arange(0.975,1.015,0.01))
    plt.errorbar(sut[0]['OH^0'], sut[5][:,0], yerr=sut[5][:,1], color='r', label=args.sut_label, linewidth=0.8)
    plt.errorbar(ref[0]['OH^0'], ref[5][:,0], yerr=ref[5][:,1], color='b', label=args.ref_label, linewidth=0.8)
    plt.ylabel('Ox/Red')
    plt.xlabel('Time (ps)')
    plt.xscale('log')
    plt.legend(loc=2, fontsize=6)

    plt.tight_layout() 

    plt.subplot(grid[7:])
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
    plt.savefig(join(args.outdir, 'Gvalue.pdf'), bbox_inches='tight')
    plt.savefig(join(args.outdir, 'Gvalue.eps'), bbox_inches='tight')
    
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



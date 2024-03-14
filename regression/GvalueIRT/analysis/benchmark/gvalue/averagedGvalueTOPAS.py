import sys, os
import numpy as np

try:
	prefix = sys.argv[1]
	outFileName = ''
	if len(sys.argv) == 3:
		outFileName = sys.argv[2]
		useUserName = True
	else:
		useUserName = False

except:
	exit(1)

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

ids = range(1,37,1)

molecules = ['OH^0', 'e_aq^-1', 'H3O^1', 'H2O2^0', 'H_2^0', 'H^0', 'OH^-1', 'HO2^0', 'O2^-1', 'O2^0']
#, 'OH^-1'] #, 'PRODUCT^0'] #, 'ThyOH^0'] #'NH3^0'] #'CH2CH2OH^0'] 
#'NH3^0', 'Cl^-1'] #'ThyOH^0'] #'NH3^0'] #'ThyOH^0'] #'NH3^0'] #'NO2^-2'] #, 'Cl^-1']

accepted = []
for i in range(len(ids)):
    if os.path.exists(prefix+'_'+str(ids[i])+'.phsp') and os.path.getsize(prefix+'_'+str(ids[i])+'.phsp') > 0:
        accepted.append(prefix+'_'+str(ids[i])+'.phsp')

gvalue, gvalue2, time = GetGValue(accepted[0], molecules) 
print 'Total accepted', len(accepted)

n = 1
for i in range(1,len(accepted)):
	agvalue, agvalue2, atime = GetGValue(accepted[i], molecules)
#	print accepted[i]
	for name, value in agvalue.iteritems():
		for j in range(len(value)):
			gvalue[name][j] += agvalue[name][j]
			gvalue2[name][j] += agvalue2[name][j]

	n += 1

if not useUserName:
	theName = prefix+'.gp'
else:
	theName = outFileName 

oFile = open(theName,'w')

for name in molecules:
	oFile.write('# %s\n' % name ) 
	print name
        for i in range(len(gvalue[name])):
		t = time[name][i]
		g = gvalue[name][i]/n
		g2 = gvalue2[name][i]
		g2 = np.sqrt( (1.0/(n-1)) * (g2/n - g*g))
		oFile.write('%1.5f  %1.5f  %1.5f  %1.5f\n' % (t, g, g2*np.sqrt(n), g2 * np.sqrt(n)))
#		if t == 1.01552: #i == 0:
#                        print '  at 1.0 ps. %1.3f  +/-  StatUnc: %1.3f%%  Std/Mean: %1.3f%%' % (g, 100*g2/g, 100*g2*np.sqrt(n)/g)
                if i == len(gvalue[name])-1:
                        print '  at 1.0 us. %1.3f  +/-  StatUnc: %1.3f or %1.3f%%  Std/Mean: %1.3f%%' % (g, g2, 100*g2/g, 100*g2*np.sqrt(n)/g)

	oFile.write('\n\n')
	
oFile.close()


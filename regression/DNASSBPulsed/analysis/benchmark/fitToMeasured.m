cd ~/Projects/IRT_DNADamage_Rezaee/analysis/benchmark/
pBR322=load('pBR322_Milligan1993.dat');
pEC=load('pEC_Milligan1993.dat');
pUC18=load('pUC18_Milligan1993.dat');
concentration = [pBR322(:,1); pEC(:,1); pUC18(:,1)];
Gvalue_SSB = [pBR322(:,2); pEC(:,2); pUC18(:,2)];
plot(concentration, Gvalue_SSB);
[ordered, index] = sort(concentration)
oc = concentration(index);
og = Gvalue_SSB(index);
concentration = oc;
Gvalue_SSB = og;
logConc = log10(concentration);
logGvalue = log10(Gvalue_SSB);

f = fit(logConc, logGvalue, 'poly1');
loglog(concentration, 10.^(f.p1 * logConc + f.p2), concentration, Gvalue_SSB, 'o')
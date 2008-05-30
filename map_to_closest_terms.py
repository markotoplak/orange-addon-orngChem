import orngCA
import sys, math, getopt
import orange
import numpy

def dist(a,b):
    return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

k=5
j=3

optlist, args=getopt.getopt(sys.argv[1:],"i:m:k:j:o1:o2")
opt=dict(optlist)
print opt
data=orange.ExampleTable(opt["-i"])
mid="-m" in opt and int(opt["-m"]) or 0
k="-k" in opt and int(opt["-k"]) or k
j="-j" in opt and int(opt["-j"]) or j

names1=[v.name for v in data.domain.variables]
mid=data.domain.getmetas().keys()[mid]
names2=[str(d[mid]) for d in data]

data=[map(lambda a:float(a) or 1e-6, e) for e in data]
ca=orngCA.CA(data, names2, names1)
row=ca.getPrincipalRowProfilesCoordinates()
col=ca.getPrincipalColProfilesCoordinates()
print row
print col
ofile="-o1" in opt and opt["-o1"] or "term-chem.tab"
file=open(ofile, "w")
file.write("term\t")
file.write("\t".join(["chemical"+str(i)+"\t"+"dist"+str(i) for i in range(k)])+"\n")
file.write("d\t"+"\t".join(["d\tc" for i in range(k)])+"\n\n")
for r, name in zip(row, names2):
    file.write(name)
    d=[(n, dist(r,c)) for c, n in zip(col, names1)]
    d.sort(lambda a,b:cmp(a[1],b[1]))
    #print d
    d=d[0:min(k,len(d))]
    for f in d:
        file.write("\t"+str(f[0])+"\t"+str(f[1]))
    file.write("\n")
file.close()


ofile="-o2" in opt and opt["-o2"] or "chem-term.tab"
file=open(ofile, "w")
file.write("chemical\t")
file.write("\t".join(["term"+str(i)+"\t"+"dist"+str(i) for i in range(j)])+"\n")
file.write("d\t"+"\t".join(["d\tc" for i in range(j)])+"\n\n")
for r, name in zip(col, names1):
    file.write(name)
    d=[(n, dist(r,c)) for c, n in zip(row, names2)]
    d.sort(lambda a,b:cmp(a[1],b[1]))
    d=d[0:min(k,len(d))]
    for f in d:
        file.write("\t"+str(f[0])+"\t"+str(f[1]))
    file.write("\n")
file.close()

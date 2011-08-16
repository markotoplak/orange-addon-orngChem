import orngCA
import sys
import orange
if len(sys.argv)>2:
    mid=int(sys.argv[1].strip("-"))
    data=orange.ExampleTable(sys.argv[2])
else:
    mid=0
    data=orange.ExampleTable(sys.argv[1])
names1=[v.name for v in data.domain.variables]
mid=data.domain.getmetas().keys()[mid]
names2=[str(d[mid]) for d in data]
#import go
#go.setDataDir("./")
#go.loadGO()
#map go term id's to term names
#if mapTermNames:
#    names2=[go.loadedGO.termDict[t.strip()].name for t in names2]
data=[map(lambda a:float(a) or 1e-6, e) for e in data]
ca=orngCA.CA(data, names2, names1)
ca.Biplot()
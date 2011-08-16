import orange
import Numeric
import sys, getopt

inputFilename="sens_profile.tab"
outputFilename="sens_profile_t.tab"
metaVarName="name"
metaindex=0

opt=dict(getopt.getopt(sys.argv[1:], "i:n:m:o:")[0])

inputFilename=opt.get("-i", None) or inputFilename
outputFilename=opt.get("-o", None) or outputFilename
metaVarName=opt.get("-n", None) or metaVarName
metaindex=opt.get("-m") or metaindex
try:
    metaindex=int(metaindex)
except:
    pass

data=orange.ExampleTable(inputFilename)
try:
    meta=data.domain.getmetas().keys()[metaindex]
except:
    meta=data.domain[metaindex]
vars=[str(var.name) for var in data.domain.variables if str(var.name)!=meta]
metavals=map(str,[e[meta] for e in data])
domain=orange.Domain([orange.FloatVariable(m) for m in metavals],0)
mid=orange.newmetaid()
domain.addmeta(mid, orange.StringVariable(metaVarName))
table=orange.ExampleTable(domain)
print vars
#print metavals
for var in vars:
    ex=orange.Example(domain)
    for e, m in zip(data, metavals):
        ex[m]=e[var]
    ex[mid]=var
    table.append(ex)
table.save(outputFilename)



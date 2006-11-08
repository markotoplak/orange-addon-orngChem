import orange
import orngChem
import sys, getopt

smilesFilename="smiles.tab"
fragmentsFilename="fragments.txt"
outputFilename="fragmentmap.tab"
binary=False

opt=dict(getopt.getopt(sys.argv[1:], "s:f:o:b")[0])

smilesFilename=opt.get("-s",None) or smilesFilename
fragmetnsFilename=opt.get("-f", None) or fragmentsFilename
outputFilename=opt.get("-o", None) or outputFilename
binary=opt.has_key("-b")

smilesData=orange.ExampleTable(smilesFilename)
smilesData=smilesData.filter(orange.Filter(lambda e:not e[1].isSpecial()))
smilesDict=dict([(str(e[0]), str(e[1])) for e in smilesData])
revSmilesDict=dict([(val, key) for key, val in smilesDict.items()])

fragments=map(lambda s:s.strip(), open(fragmetnsFilename).read().split("\n"))
fragmentMap=orngChem.map_fragments(fragments, smilesDict.values(), binary)

vars=[orange.FloatVariable(frag) for frag in fragments]
mid=orange.newmetaid()
cvar=orange.StringVariable("chemical")
domain=orange.Domain(vars,0)
domain.addmeta(mid,cvar)
table=orange.ExampleTable(domain)

for chem, fragmap in fragmentMap.items():
    ex=orange.Example(domain)
    for fragment, val in fragmap.items():
        ex[fragment]=val
    ex[mid]=chem
    table.append(ex)

table.save(outputFilename)
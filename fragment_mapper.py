import orange
import orngChem_Old as orngChem
import sys, getopt

smilesFilename="smiles.tab"
fragmentsFilename="fragments.txt"
outputFilename="fragmentmap.tab"
binary=False

opt=dict(getopt.getopt(sys.argv[1:], "s:f:o:a:b")[0])

smilesFilename=opt.get("-s",None) or smilesFilename
fragmetnsFilename=opt.get("-f", None) or fragmentsFilename
outputFilename=opt.get("-o", None) or outputFilename
attrName=opt.get("-a", 1)
binary=opt.has_key("-b")

smilesData=orange.ExampleTable(smilesFilename)
smilesData=smilesData.filter(orange.Filter(lambda e:not e[attrName].isSpecial()))
#smilesDict=dict([(str(e[0]), str(e[1])) for e in smilesData])
#revSmilesDict=dict([(val, key) for key, val in smilesDict.items()])
smilesCodes=[str(e[attrName]) for e in smilesData]

fragments=map(lambda s:s.strip(), open(fragmetnsFilename).read().split("\n"))
fragmentMap=orngChem.map_fragments(fragments, smilesCodes, binary)

vars=[orange.FloatVariable(frag) for frag in fragments]
mid=orange.newmetaid()
cvar=orange.StringVariable("chemical")
vars=smilesData.domain.attributes+vars+(smilesData.domain.classVar and [smilesData.domain.classVar] or [])

domain=orange.Domain(vars,1)
domain.addmetas(smilesData.domain.getmetas())
#domain.addmeta(mid,cvar)
table=orange.ExampleTable(domain)
for e in smilesData:
    ex=orange.Example(domain)
    for v in smilesData.domain.variables+smilesData.domain.getmetas().values():
        ex[v]=e[v]
    fragmap=fragmentMap[str(e[attrName])]
    for frag, val in fragmap.items():
        ex[frag]=val
    table.append(ex)
"""
for chem, fragmap in fragmentMap.items():
    ex=orange.Example(domain)
    for fragment, val in fragmap.items():
        ex[fragment]=val
    ex[mid]=chem
    table.append(ex)
"""
table.save(outputFilename)
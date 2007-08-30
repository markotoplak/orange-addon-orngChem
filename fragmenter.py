import orange
import orngChem_Old as orngChem
import getopt
import sys

smilesFilename="smiles.tab"
fragmentFilename="fragments.txt"
freq=0.4
opt=dict(getopt.getopt(sys.argv[1:], "s:f:o:a:")[0])
smilesFilename=opt.get("-s", None) or smilesFilename
fragmentFilename=opt.get("-o", None) or fragmentFilename
freq=float(opt.get("-f", freq))
attrName=opt.get("-a",None) or 1

if smilesFilename.endswith(".tab"):
    smilesData=orange.ExampleTable(smilesFilename)
    smiles=map(str, [e[attrName] for e in smilesData if not e[attrName].isSpecial()])
else:
    smiles=map(lambda s:s.strip(), open(smilesFilename).read().split("\n"))
    
fragments=orngChem.find_fragments(smiles,freq)
file=open(fragmentFilename, "w")
file.write("\n".join(fragments))
file.close()




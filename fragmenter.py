import orange
import orngChem
import getopt
import sys

smilesFilename="smiles.tab"
fragmentFilename="fragments.txt"
freq=0.4
opt=dict(getopt.getopt(sys.argv[1:], "s:f:o:")[0])
smilesFilename=opt.get("-s", None) or smilesFilename
fragmentFilename=opt.get("-o", None) or fragmentFilename
freq=float(opt.get("-f", None)) or freq

if smilesFilename.endswith(".tab"):
    smilesData=orange.ExampleTable(smilesFilename)
    smiles=map(str, [e[1] for e in smilesData if not e[1].isSpecial()])
else:
    smiles=map(lambda s:s.strip(), open(smilesFilename).read().split("\n"))
    
fragments=orngChem.find_fragments(smiles,freq)
file=open(fragmentFilename, "w")
file.write("\n".join(fragments))
file.close()




import orange
import go
import orngChem_Old as orngChem
import sys, getopt

smilesFilename="smiles.tab"
sensFilename="sens.tab"
fragmentsFilename="fragments.txt"
outputFilename="sens_profile.tab"
subsetFilename=None
goDataDir=None
slimsSubset=True
aspect="P"
fragmentBased=True

opt=dict(getopt.getopt(sys.argv[1:], "s:S:b:f:a:g:o:lm")[0])

smilesFilename=opt.get("-s",None) or smilesFilename
fragmetnsFilename=opt.get("-f", None) or fragmentsFilename
sensFilename=opt.get("-S", None) or sensFilename
outputFilename=opt.get("-o", None) or outputFilename
subsetFilename=opt.get("-b", None)
goDataDir=opt.get("-g", None)
slimsSubset=opt.has_key("-l")
fragmentBased=not opt.has_key("-m")
aspect=opt.get("-a", None) or aspect

smilesData=orange.ExampleTable(smilesFilename)
smilesData=smilesData.filter(orange.Filter(lambda e:not e[1].isSpecial()))
smilesDict=dict([(str(e[0]), str(e[1])) for e in smilesData])
revSmilesDict=dict([(val, key) for key, val in smilesDict.items()])

sensData=orange.ExampleTable(sensFilename)
sensDict=dict([(str(e[0]), e) for e in sensData])

genes=map(str, [e[0] for e in sensData])

if goDataDir:
    go.setDataDir(goDataDir)
go.loadGO()
go.loadAnnotation("sgd")
go.setSlims("goslim_yeast")

terms=go.findTerms(genes, slimsOnly=slimsSubset, aspect=aspect, reportEvidence=False)

if subsetFilename:
    file=open(subsetFilename)
    subset=map(lambda s:s.strip(), file.read().split("\n"))
    l=filter(lambda t:t[0] in subset or go.loadedGO.termDict[t[0]].name in subset, terms.items())
    terms=dict(l)

if fragmentBased:
    #print terms
    fragments=map(lambda s:s.strip(), open(fragmetnsFilename).read().split("\n"))
    fragmentMap=orngChem.map_fragments(fragments, smilesDict.values())
    domain=orange.Domain([orange.FloatVariable(frag) for frag in fragments],0)
    mid1=orange.newmetaid()
    mid2=orange.newmetaid()
    domain.addmeta(mid1, orange.StringVariable("GOTerm id"))
    domain.addmeta(mid2, orange.StringVariable("GOTerm name"))
    table=orange.ExampleTable(domain)
    matrix=[]
    
    for term, genes in terms.items():
        ex=orange.Example(domain)
        for frag in fragments:
            chemicals=filter(lambda a: fragmentMap[a][frag], fragmentMap.keys())
            avgSens=0.0
            for g in genes:
                for c in chemicals:
                    avgSens+=float(sensDict[g][revSmilesDict[c]])
            avgSens/=len(genes)*len(chemicals)
            ex[frag]=avgSens
        ex[mid1]=term
        ex[mid2]=go.loadedGO.termDict[term].name
        table.append(ex)
    table.save(outputFilename)
else:
    domain=orange.Domain(sensData.domain.variables[1:],0)
    mid1=orange.newmetaid()
    mid2=orange.newmetaid()
    domain.addmeta(mid1, orange.StringVariable("GOTerm id"))
    domain.addmeta(mid2, orange.StringVariable("GOTerm name"))
    table=orange.ExampleTable(domain)
    for term, genes in terms.items():
        ex=orange.Example(domain)
        for chem in sensData.domain.variables[1:]:
            avgSens=0.0
            for g in genes:
                avgSens+=float(sensDict[g][chem])
            avgSens/=len(genes)
            ex[chem]=avgSens
        ex[mid1]=term
        ex[mid2]=go.loadedGO.termDict[term].name
        table.append(ex)
    table.save(outputFilename)

            
                
    
    

        





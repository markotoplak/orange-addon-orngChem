import go
go.loadAnnotation()
go.loadGO()
from sets import Set
terms=["GO:0030036",
       "GO:0043037",
       "GO:0030528",
       "GO:0016044",
       "GO:0006873",
       "GO:0006281",
       "GO:0042026",
       "GO:0007047",
       "GO:0000754",
       "GO:0008202",
       "GO:0007017"]

booneGenes=open("genes.txt").read().split("\n")

file=open("terms_with_genes.txt","w")
#for each term find all genes that map to this term and those that are from boone's paper
for term in terms:
    genes=go.findGenes([term]).keys()
    #test
    t=go.GOTermFinder(genes, aspect="P")
    #t=go.findTerms(genes)
    #if term not in t.keys():
    #    print "!!!!!ERROR!!!!!", term, t
    #print genes
    genesBoone=filter(lambda a: go.mapGeneName(a) in genes, booneGenes)
    #t1=go.GOTermFinder(genesBoone, aspect="P")
    #print genesBoone
    t=go.findTerms(genesBoone, aspect=["P"])
    #print t
    #print genesBoone
    #print Set(t).difference(Set(t1)), Set(t1).difference(Set(t))
    #if term not in t.keys():
    #    print "!!!!!ERROR BOONE!!!!!", term, go.loadedGO.termDict[term].name
    #print genesBoone
    file.write(term+"\t"+",".join(genes)+"\t"+",".join(genesBoone)+"\n")
file.close()
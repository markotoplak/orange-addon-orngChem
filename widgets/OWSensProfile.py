"""
<name>Sensitivity profile analysis</name>
"""

import orange
import obiChem
from OWWidget import *
import OWGUI
import obiGO


class OWSensProfile(OWWidget):
    settingsList = ["aspectNum", "fragFreq", "smilesFileName"]
    
    def __init__(self, parent=None, signalManager=None, name="Sens profile"):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea=False)
        self.inputs = [("Sensitivity", ExampleTable, self.sensData)]
        self.outputs = [("Molecule fragmetns", ExampleTable), ("Slim-based molecular profiles", ExampleTable), ("Molecule-based slim profiles", ExampleTable), ("Slim-based fragment profiles", ExampleTable), ("Fragment-based slim profiles", ExampleTable)]

        self.aspectNum = 0
        self.fragFreq = 0.1
        self.smilesFileName = ""
        self.slimsOnly = 1
        self.maxTerms = 100
        
        self.loadSettings()
        
        OWGUI.radioButtonsInBox(self.controlArea, self, "aspectNum",
                                ["Process", "Function", "Component"], "Aspect")
        
        OWGUI.checkBox(self.controlArea, self, "slimsOnly", "Slims only",
                       tooltip="Search for slim GO terms only")
        
        OWGUI.spin(self.controlArea, self, "maxTerms", 10, 1000, box="Number of terms",
                   tooltip="Maximum number of GO terms")
        
        self.findTermsButton = OWGUI.button(self.controlArea, self, "Find terms",
                                            tooltip="Run term search",
                                            callback=self.findTerms)
        
        OWGUI.button(self.controlArea, self, "Load selected terms",
                     tooltip="Load selected terms from a file",
                     callback=self.loadSelectedTerms)
        
        OWGUI.doubleSpin(self.controlArea, self, "fragFreq", 0.05, 0.95, 0.05,
                         box="Fragment frequency",
                         tooltip="Minimum fragment frequency")
        
        self.findFragmentsButton = OWGUI.button(self.controlArea, self, "Find fragments",
                                                tooltip="Run fragment search", 
                                                callback=self.findFragments)
        
        OWGUI.button(self.controlArea, self, "Load smiles",
                     tooltip="Load fragment SMILES descriptors from a file",
                     callback=self.loadSmiles)
        
        self.commitButton = OWGUI.button(self.controlArea, self, "Commit",
                                         callback=self.commitData)
        
        self.findTermsButton.setEnabled(0)
        self.findFragmentsButton.setEnabled(0)
        self.commitButton.setEnabled(0)
        self.resize(100,100)
        
        self.fragments = []
        self.sensDict = {}
        self.smilesDict = {}
        self.reverseSmilesDict = {}
        self.data = None
        self.fragmentMap = None
        self.terms = None
        self.ontology = None
        self.annotations = None
        self.selectedTerms = []
        self.lastDir = os.path.expanduser("~/")
        
        
        self.loadSmiles(self.smilesFileName)
        
        self.loadSettings()
        
    def findFragments(self):
        if self.smilesDict:
            smiles = [self.smilesDict[chem] for chem in self.chemicals if chem in self.smilesDict]
            fragments = []
            miner = obiChem.FragmentMiner(smiles, [], minSupport=self.fragFreq)
            for fragment in miner.SearchIterator():
                fragments.append(fragment.ToSmiles())
                
            self.fragments = fragments
#            self.fragmentMap = obiChem.map_fragments(self.fragments, smiles)
            self.fragmentMap = map_fragments(self.fragments, smiles)

            if self.terms:
                self.commitButton.setEnabled(1)

    def findTerms(self):
        if not self.data:
            return 
        genes = [str(e[0]) for e in self.data]
        aspect = ["P", "F", "C"][self.aspectNum]
        if self.ontology is None:
            self.ontology = obiGO.Ontology()
            self.ontology.SetSlimsSubset("goslim_yeast")
        
        if self.annotations is None:
            self.annotations = obiGO.Annotations("yeast", ontology=self.ontology)
        
        self.progressBarInit()
        
        self.terms = self.annotations.GetEnrichedTerms(genes, slimsOnly=self.slimsOnly, aspect=aspect, progressCallback=self.progressBarSet)
#        self.terms=go.findTerms(genes, slimsOnly=self.slimsOnly, aspect=[aspect],reportEvidence=False, progressCallback=self.progressBarSet)
        
        self.progressBarFinished()
        self.filterTerms()
        if self.fragments:
            self.commitButton.setEnabled(1)

    def filterTerms(self):
        ll = []
        for key, val in self.terms.items():
            if self.selectedTerms and key in self.selectedTerms:
                ll.append((key, val))
            elif not self.selectedTerms:
                ll.append((key, val))
        #ll.sort(lambda a,b:-cmp(a[1][1], b[1][1]))
        #ll=ll[:max(len(ll),self.maxTerms)]
        self.terms=dict(ll)

    def sendFragments(self):
        vars = [orange.FloatVariable(frag) for frag in self.fragments]
        mid = orange.newmetaid()
        chemVar = orange.StringVariable("chemical name")
        domain = orange.Domain(vars,0)
        domain.addmeta(mid, chemVar)
        table = orange.ExampleTable(domain)
        for chem, map in self.fragmentMap.items():
            val = [map[frag] for frag in self.fragments]
            e = orange.Example(domain, val)
            e[mid] = chem
            table.append(e)
        self.send("Molecule fragmetns", table)

    def sendFragmentTerms(self):
        matrix=[]
        for frag in self.fragments:
            vec=[]
            for term in self.terms.keys():
                genes=self.terms[term][0]
                chem=filter(lambda a:self.fragmentMap[a][frag], self.fragmentMap.keys())
                avgSens=0.0
                for g in genes:
                    for c in chem:
                        avgSens+=self.sensDict[g][self.reverseSmilesDict[c]]
                avgSens/=len(genes)*len(chem)
                vec.append(avgSens)
            matrix.append(vec)
            
        vars=[orange.FloatVariable(term) for term in self.terms.keys()]
        mid=orange.newmetaid()
        domain=orange.Domain(vars,0)
        domain.addmeta(mid, orange.StringVariable("fragment"))
        table=orange.ExampleTable(domain)
        for frag, vec in zip(self.fragments, matrix):
            e=orange.Example(domain, vec)
            e[mid]=frag
            table.append(e)
        self.send("Slim-based fragment profiles", table)
        
        import Numeric
        matrix=Numeric.transpose(matrix)
        vars=[orange.FloatVariable(frag) for frag in self.fragments]
        mid=orange.newmetaid()
        mid1=orange.newmetaid()
        domain=orange.Domain(vars,0)
        domain.addmeta(mid, orange.StringVariable("term id"))
        domain.addmeta(mid1, orange.StringVariable("term name"))
        table=orange.ExampleTable(domain)
        for term_id, vec in zip(self.terms.keys(),matrix):
            e=orange.Example(domain, list(vec))
            term = self.ontology[term_id]
            e[mid]=term_id
            e[mid1]=term.name
            table.append(e)
        self.send("Fragment-based slim profiles", table)

    def sendMoleculeTerms(self):
        matrix=[]
        for chem in self.data.domain.variables[1:]:
            vec=[]
            for term in self.terms.keys():
                genes=self.terms[term][0]
                avgSens=0.0
                for g in genes:
                    avgSens+=self.sensDict[g][chem]
                avgSens/=len(genes)
                vec.append(avgSens)
            matrix.append(vec)
            
        vars=[orange.FloatVariable(term) for term in self.terms.keys()]
        mid=orange.newmetaid()
        domain=orange.Domain(vars,0)
        domain.addmeta(mid, orange.StringVariable("chemical name"))
        table=orange.ExampleTable(domain)
        for chem, vec in zip(self.chemicals,matrix):
            e=orange.Example(domain, vec)
            e[mid]=chem
            table.append(e)
        self.send("Slim-based molecular profiles", table)
        
        import Numeric
        matrix=Numeric.transpose(matrix)
        vars=[orange.FloatVariable(chem) for chem in self.chemicals]
        mid=orange.newmetaid()
        mid1=orange.newmetaid()
        domain=orange.Domain(vars,0)
        domain.addmeta(mid, orange.StringVariable("term id"))
        domain.addmeta(mid1, orange.StringVariable("term name"))
        table=orange.ExampleTable(domain)
        for term_id, vec in zip(self.terms.keys(),matrix):
            e=orange.Example(domain, list(vec))
            term = self.ontology[term_id]
            e[mid]=term_id
            e[mid1]=term.name
            table.append(e)
        self.send("Molecule-based slim profiles", table)

    def sensData(self, data=None):
        self.data=data
        self.send("Molecule fragmetns", None)
        self.send("Slim-based molecular profiles", None)
        self.send("Molecule-based slim profiles", None)
        self.send("Slim-based fragment profiles", None)
        self.send("Fragment-based slim profiles", None)
        if self.data:
            self.findTermsButton.setEnabled(1)
            self.sensDict={}
            geneAttr = 0
            for e in data:
                self.sensDict[str(e[geneAttr])] = e
#                self.sensDict[go.mapGeneName(str(e[0]))]=e
            self.chemicals=[var.name for var in data.domain.variables[1:]]
        else:
            self.findTermsButton.setEnabled(0)
            self.commitButton.setEnabled(0)
            self.sensDict={}
            self.chemicals=[]
            
    def commitData(self):
        if self.fragments:
            self.sendFragments()
            self.sendFragmentTerms()
        self.sendMoleculeTerms()
        
    def loadSmiles(self, filename=""):
        if not filename:
            filename = str(QFileDialog.getOpenFileName(self, "Choose a smiles file", self.lastDir, "Text files (*.txt *.tab)"))
        try:
            data = orange.ExampleTable(filename)
            self.smilesDict = [(str(e[0]), str(e[1])) for e in data if not e[1].isSpecial()]
        except:
            file = open(filename)
            data = file.read()
            lines = data.split("\n")
            self.smilesDict = [l.split("\t") for l in lines]
        self.smilesDict=dict([(a[0].strip(), a[1].strip()) for a in self.smilesDict if a[1].strip() or a[1].strip()!="?"])
        self.reverseSmilesDict = dict([(val, key) for key, val in self.smilesDict.items()])
        self.smilesFileName = filename
        self.lastDir = os.path.dirname(filename)
        self.findFragmentsButton.setEnabled(1)

    def loadSelectedTerms(self):
        filename=str(QFileDialog.getOpenFileName(self, "Choose file containing terms", self.lastDir, "Text file (*.txt)"))
        data=open(filename)
        data=data.read().split("\n")
        self.selectedTerms = data
        self.lastDir = os.path.dirname(filename)
        
        
def map_fragments(fragments, smiles):
    from collections import defaultdict
    map = defaultdict(dict)
    for mol in smiles:
        for frag in fragments:
            map[mol][frag] = obiChem.ContaindIn(mol, frag)
    return dict(map)


if __name__=="__main__":
    import sys
    app=QApplication(sys.argv)
    w=OWSensProfile()
    d=orange.ExampleTable("../sens.tab")
    w.sensData(d)
    w.show()
    app.exec_()
    w.saveSettings()
    
            
        
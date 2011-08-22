"""<name>Molecule Match</name>
<description>Selection of molecules based on SMARTS patterns</description>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact> 
<priority>2020</priority>"""

import orange
import obiChem
from OWWidget import *
import OWGUI
import sys
import os

try:
    import pybel
    HAVE_PYBEL = True
except ImportError:
    HAVE_PYBEL = False
    

class OWMoleculeMatch(OWWidget):
    settingsList = ["recentFragments"]
    def __init__(self, parent=None, signalManager=None, name="Molecule Match"):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea=False)
        
        self.inputs = [("Molecules", ExampleTable, self.SetMoleculeTable)]
        self.outputs = [("Selected molecules", ExampleTable), ("All molecules", ExampleTable)]
        
        self.fragment = ""
        self.recentFragments = []
        self.comboSelection = 0 
        self.smilesAttr = 0
        self.smilesAttrList = []
        self.data = None
        self.loadSettings()

        ##GUI        
        self.lineEdit = OWGUI.lineEdit(self.controlArea, self, "fragment",
                                       box="Pattern",
                                       tooltip="A SMARTS pattern to filter the examples with.",
                                       callback=self.LineEditSelect)
        
        OWGUI.separator(self.controlArea)
        
        self.fragmentComboBox = OWGUI.comboBox(self.controlArea, self, 
                                   "comboSelection", "Recent Patterns",
                                   tooltip="Select a recently used pattern.",
                                   items=self.recentFragments,
                                   callback=self.RecentSelect)
        
        OWGUI.separator(self.controlArea)
        
        self.smilesAttrCombo = OWGUI.comboBox(self.controlArea, self, "smilesAttr",
                                              "Examples SMILES attribute",
                                              tooltip="Select the attribute in the input example \
                                              table that stores the SMILES molecule description", 
                                              callback=self.Process)
        
        OWGUI.rubber(self.controlArea)
        self.resize(100, 100)
        
        if not HAVE_PYBEL:
            self.error(10, "This widget requires pybel module (part of openbabel python extension).")

    def SetMoleculeTable(self, data=None):
        self.data = data
        if data:
            self.SetSmilesAttrCombo()
            self.Process()

    def SetSmilesAttrCombo(self):
        fragmenter = obiChem.Fragmenter()
        attr = [fragmenter.FindSmilesAttr(self.data)]
        self.smilesAttrList = attr
        if self.smilesAttr > len(attr)-1:
            self.smilesAttr = 0
        self.smilesAttrCombo.clear()
        icons = OWGUI.getAttributeIcons()
        for a in attr:
            self.smilesAttrCombo.addItem(QIcon(icons[a.varType]), a.name)
            
    def LineEditSelect(self):
        self.Process()
        self.recentFragments.insert(0, self.fragment)
        self.fragmentComboBox.insertItem(0, self.fragment)
        if len(self.recentFragments) > 10:
            self.recentFragments = self.recentFragments[:10]

    def RecentSelect(self):
        self.fragment = self.recentFragments[self.comboSelection]
        self.Process()

    def Process(self):
        if not HAVE_PYBEL:
            return
        from functools import partial
        newVar = orange.Variable.make(str(self.fragment), orange.VarTypes.Continuous)
#        newVar = orange.FloatVariable(str(self.fragment), numberOfDecimals=0)
        
        def getVal(var, smilesAttr, smarts ,example, returnWhat):
            mol = pybel.readstring("smi", str(example[smilesAttr]))
            if smarts.findall(mol):
                return var(1)
            else:
                return var(0)
            
        newVar.getValueFrom  = partial(getVal, newVar, self.smilesAttrList[self.smilesAttr], pybel.Smarts(str(self.fragment)))
        domain = orange.Domain(self.data.domain.attributes + [newVar], self.data.domain.classVar)
        domain.addmetas(self.data.domain.getmetas())
        
        selected = []
        all = []
        for example in self.data:
            ex = orange.Example(domain, example)
            if ex[newVar]==newVar(1):
                selected.append(example)
            all.append(ex)
        all = orange.ExampleTable(domain, all)
        if selected:
            selected = orange.ExampleTable(selected)
        else:
            selected = orange.ExampleTable(self.data.domain)
        self.send("Selected molecules", selected)
        self.send("All molecules", all)

if __name__=="__main__":
    app = QApplication(sys.argv)
    w = OWMoleculeMatch()
    w.show()
    data = orange.ExampleTable("E:/orangecvs/Chemo/smiles.tab")
    w.SetMoleculeTable(data)
    app.exec_()
    w.saveSettings()
    
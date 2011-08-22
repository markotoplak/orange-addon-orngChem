"""<name>Fragmenter</name>
<description>Find frequent fragments in molecules</description>
"""

from OWWidget import *
import OWGUI

import pybel
import obiChem 

class OWFragmenter(OWWidget):
    settingsList = ["minSupport", "maxSupport"]
    contextHandlers = {"active": DomainContextHandler("active", ["activeSmilesAttr"]),
                       "inactive": DomainContextHandler("inactive", ["inactiveSmilesAttr"])
                       }
    
    def __init__(self, parent=None, signalManager=None, name="Fragmenter"):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea=False)
        
        self.inputs = [("Active chemicals", ExampleTable, self.setActive), ("Inactive chemicals", ExampleTable, self.setInactive)]
        self.outputs = [("Fragments", ExampleTable)]
        
        self.minSupport = 0.2
        self.maxSupport = 0.2
        self.activeSmilesAttr = 0
        self.inactiveSmilesAttr = 0
        
        self.loadSettings()
        #  GUI
        
        box = OWGUI.widgetBox(self.controlArea, "Active Chemicals Set") 
        OWGUI.doubleSpin(box, self, "minSupport", 0.05, 0.95, step=0.05, 
                         label="Min. active frequency", 
                         tooltip="Minimal fragment frequency in the active chemicals set.",
                         callback=self.updateFreq)
        
        self.activeSmilesAttrCB = OWGUI.comboBox(box, self, "activeSmilesAttr",
                                                 label="SMILES attribute",
                                                 callback=self.updateAttrs)
        
        box = OWGUI.widgetBox(self.controlArea, "Inactive Chemicals Set")
        OWGUI.doubleSpin(box, self, "maxSupport", 0.05, 0.95, step=0.05,
                         label="Max. inactive frequency",
                         tooltip="Maximal fragment frequency in the inactive chemicals set.",
                         callback=self.updateFreq)
        
        self.inactiveSmilesAttrCB = OWGUI.comboBox(box, self, "inactiveSmilesAttr",
                                                   label="SMILES attribute",
                                                   callback=self.updateAttrs)
        
        OWGUI.button(self.controlArea, self, "Run",
                     callback=self.run)
        
        self.activeData = None
        self.inactiveData = None
        self.activeDataAttrs = []
        self.inactiveDataAttrs = []
        
        self.resize(100, 100)
        
    def updateFreq(self):
        pass
    
    def updateAttrs(self):
        #TODO: update info box
        pass
    
    def setActive(self, data=None):
        self.activeData = data
        self.activeSmilesAttrCB.clear()
        if data is not None:
            fragmenter = obiChem.Fragmenter()
            attrs = [fragmenter.FindSmilesAttr(data)]
            self.activeDataAttrs = attrs
            self.activeSmilesAttrCB.addItems([attr.name for attr in attrs])
            
    def setInactive(self, data=None):
        self.inactiveData = data
        self.inactiveSmilesAttrCB.clear()
        if data is not None:
            fragmenter = obiChem.Fragmenter()
            attrs = [fragmenter.FindSmilesAttr(data)]
            self.inactiveDataAttrs = attrs
            self.inactiveSmilesAttrCB.addItems([attr.name for attr in attrs])
            
    def run(self):
        active = self.getActiveSmiles()
        inactive = self.getInactiveSmiles()
        if not active:
            return 
        miner = obiChem.FragmentMiner(active, inactive, minSupport=self.minSupport,
                                      maxSupport=self.maxSupport,
                                      addWholeRings=True,
                                      canonicalPruning=True,
                                      findClosed=True)
        fragments = []
        import zlib, math
        estimation = len(zlib.compress(" ".join(active), 9)) # estimation of fragments
        estimation = math.log(len(set("".join(active)).difference(set("()[]"))) ** estimation) * (1.0 - self.minSupport) / 100
        
        self.progressBarInit()
        for i, fragment in enumerate(miner.SearchIterator()):
            fragments.append((fragment.ToSmiles(), fragment.Support()))
            self.progressBarSet(min(100.0 * i / estimation, 100.0))
        self.progressBarFinished()
        
        def getVal(var, smilesAttr, smarts ,example, returnWhat):
            mol = pybel.readstring("smi", str(example[smilesAttr]))
            if smarts.findall(mol):
                return var(1)
            else:
                return var(0)
            
        fragmentVars = []
        from functools import partial
        for fragment, support in fragments:
            var = orange.FloatVariable(fragment, numberOfDecimals=0)
            var.getValueFrom = partial(getVal, var, self.activeDataAttrs[self.activeSmilesAttr], pybel.Smarts(fragment))
            fragmentVars.append(var)
        newdomain = orange.Domain(self.activeData.domain.attributes + fragmentVars, self.activeData.domain.classVar)
        newdomain.addmetas(self.activeData.domain.getmetas())
        
        data = orange.ExampleTable(newdomain, self.activeData)
        self.send("Fragments", data)
        
    def getActiveSmiles(self):
        if self.activeData and self.activeDataAttrs:
            attr = self.activeDataAttrs[self.activeSmilesAttr]
            return [str(ex[attr]) for ex in self.activeData if not ex[attr].isSpecial()]
        else:
            return []
        
    def getInactiveSmiles(self):
        if self.inactiveData and self.inactiveDataAttrs:
            attr = self.inactiveDataAttrs[self.inactiveSmilesAttr]
            return [str(ex[attr]) for ex in self.inactiveData if not ex[attr].isSpecial()]
        else:
            return []
        
        
        
            
            
        
        
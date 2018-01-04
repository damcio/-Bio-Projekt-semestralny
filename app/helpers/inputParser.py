import subprocess
from defs import ROOT_DIR

import os
from Bio.PDB import PDBList, PDBParser, MMCIFParser

nonoWords = ['MG', 'MN', 'HOH']

class structureInfo():
    def __init__(self):
        self.basePairs = list()

def onlineInput(structureName, fileFormat=None):
    pdbl = PDBList()

    return pdbl.retrieve_pdb_file(pdb_code=structureName, file_format=fileFormat)

def readBasePairs(annotatedOutput):
    out = iter(str(annotatedOutput).split('\n'))
    pairs = list(list())
    for line in out:
        tmpList = list()
        if 'Base-pairs' in line.strip():
            tmpLine = next(out).strip()
            while tmpLine and 'Residue conformations' not in tmpLine:
                tmpList.append(tmpLine.split()[0])
                tmpLine = next(out).strip()
            pairs.append(tmpList)
    return pairs

def readStrand(annotatedOutput):
   pdbp = PDBParser()
   struct = pdbp.get_structure('file', annotatedOutput)
   strands = []

   for model in struct.get_models():
       strand = []
       for res in model.get_residues():
           tempName = res.resname.lstrip()
           if tempName not in nonoWords:
               if len(tempName) > 1:
                   tempName = tempName[-1]
               strand.append(tempName)
       strands.append(strand)

   return strands


def annotate(filename=None):
    stdout = subprocess.run([ROOT_DIR+'/ext/MC_Annotate/MC-Annotate', filename], stdout=subprocess.PIPE, universal_newlines=True).stdout
    with open(ROOT_DIR+'/ext/MC_Annotate/output/'+filename.split('/')[-1]+'.out', 'w+') as f:
        f.write(stdout)
    return readBasePairs(stdout)


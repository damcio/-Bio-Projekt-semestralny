import subprocess
from defs import ROOT_DIR

import os
from Bio.PDB import PDBList, PDBParser, MMCIFParser

def onlineInput(structureName, fileFormat=None):
    pdbl = PDBList()

    return pdbl.retrieve_pdb_file(pdb_code=structureName, file_format=fileFormat)

def readModels(annotatedOutput):
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

def annotate(filename=None):
    stdout = subprocess.run([ROOT_DIR+'/ext/MC_Annotate/MC-Annotate', filename], stdout=subprocess.PIPE, universal_newlines=True).stdout
    return readModels(stdout)


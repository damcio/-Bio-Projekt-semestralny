import subprocess
from defs import ROOT_DIR

import os
from Bio.PDB import PDBList, PDBParser, MMCIFParser

nonoWords = ['MG', 'MN', 'HOH']


class structureInfo():
    def __init__(self):
        self.basePairs = list()


def onlineInput(structureName, fileFormat=None):
    if not fileFormat:
        fileFormat = 'pdb'
    pdbl = PDBList()

    return pdbl.retrieve_pdb_file(pdb_code=structureName, file_format=fileFormat,
                                  pdir=ROOT_DIR + '/downloadedStructures/')


def readBasePairs(annotatedOutput):
    out = iter(str(annotatedOutput).split('\n'))
    pairs = list(list())
    for line in out:
        tmpList = list()
        if 'Base-pairs' in line.strip():
            tmpLine = next(out).strip().split()
            while tmpLine and 'Residue' not in tmpLine[0]:
                try:
                    if 'Ww/' in tmpLine[3] and tmpLine[6] == 'cis':
                        tmpPair = tmpLine[0].split('-')
                        tmpList.append((int(tmpPair[0][1:]), int(tmpPair[1][1:])))
                except IndexError:
                    pass
                tmpLine = next(out).strip().split()
            pairs.append(tmpList)
    return pairs


def readStrand(filePath):
    pdbp = PDBParser()
    struct = pdbp.get_structure('file', filePath)
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


def traceDepth(pair, stacks, depth):
    for stackPair in stacks[depth]:
        if pair[1] > stackPair[1]:
            depth += 1
            return depth
    return depth


def makeDotNotation(strand, basepairs):
    brackets = [('(', ')'), ('[', ']')]
    output = ['.' for _ in strand]
    stacks = [[], []]
    depth = 0
    for pair in basepairs:
        tempDepth = depth
        depth = traceDepth(pair, stacks, depth)
        while tempDepth > depth:
            tempDepth = depth
            depth = traceDepth(pair, stacks, depth)

        stacks[depth].append(pair)
        depth = 0

    for i, stack in enumerate(stacks):
        for pair in stack:
            output[pair[0]-1], output[pair[1]-1] = brackets[i][0], brackets[i][1]

    return "".join(output)

def annotate(filename=None):
    stdout = subprocess.run([ROOT_DIR + '/ext/MC_Annotate/MC-Annotate', filename], stdout=subprocess.PIPE,
                            universal_newlines=True).stdout
    with open(ROOT_DIR + '/ext/MC_Annotate/output/' + filename.split('/')[-1] + '.out', 'w+') as f:
        f.write(stdout)
    return readBasePairs(stdout)

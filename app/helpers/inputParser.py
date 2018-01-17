import subprocess
from collections import defaultdict, namedtuple

from defs import ROOT_DIR

import os
from Bio.PDB import PDBList, PDBParser, MMCIFParser

nonoWords = ['MG', 'MN', 'HOH', 'K']


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


def readModels(filePath):
    pdbp = PDBParser()
    struct = pdbp.get_structure('file', filePath)
    models = list()

    for model in struct.get_models():
        for chain in model.get_chains():
            strand = defaultdict(list)
            chainName = chain.id
            for res in chain.get_residues():
                tempName = res.resname.lstrip()
                if tempName not in nonoWords:
                    if len(tempName) > 1:
                        tempName = tempName[-1]
                    strand[chainName].append(tempName)
            models.append(strand)

    with open(filePath, 'r+') as file:
        for line in file:
            if 'MISSING RESIDUES' in line:
                missing_residues = []
                missing_res_tuple = namedtuple('missing_res_tuple', ['chain', 'residue', 'index'])
                for _ in range(5):
                    next(file)
                temp_line = file.readline().rstrip().split()
                while '470' not in temp_line[1]:
                    missing_residues.append(missing_res_tuple(temp_line[3], temp_line[2], int(temp_line[4])))
                    temp_line = file.readline().rstrip().split()

                for model in models:
                    for res in missing_residues:
                        if res.chain in model.keys():
                            model[res.chain].insert(res.index - 1,
                                                    res.residue)
                break

    return models


def traceDepth(pair, stacks, depth):
    for stackPair in stacks[depth]:
        if pair[1] > stackPair[1]:
            if pair[0] < stackPair[1]:
                depth += 1
                return depth
    return depth


def makeDotNotation(strand, basepairs):
    brackets = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]
    output = ['.' if i != '-' else '-' for i in strand]
    stacks = [[], [], [], [], []]
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
            output[pair[0] - 1], output[pair[1] - 1] = brackets[i][0], brackets[i][1]

    return "".join(output)


def annotate(filename=None):
    stdout = subprocess.run([ROOT_DIR + '/ext/MC_Annotate/MC-Annotate', filename], stdout=subprocess.PIPE,
                            universal_newlines=True).stdout
    with open(ROOT_DIR + '/ext/MC_Annotate/output/' + filename.split('/')[-1] + '.out', 'w+') as f:
        f.write(stdout)
    return readBasePairs(stdout)

#!/usr/bin/python3
import os

import pytest
from Bio.PDB.PDBParser import PDBParser

from app.helpers import inputParser
from defs import ROOT_DIR


def test_pdbListDownload():
    structureName = '1FAT'
    file = inputParser.onlineInput(structureName)
    assert structureName.lower() in file
    with open(file) as f:
        assert structureName in f.readline()


def test_pdbListDownloadWithFormat():
    structureName = '2LBK'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)
    assert structureName.lower() in file
    assert fileFormat in file
    with open(file) as f:
        assert structureName in f.readline()


def test_annotateDownloadedfile():
    structureName = '2lbk'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    models = inputParser.annotate(file)

    assert len(models) == 8


def test_strandLensMultipleModels():
    structureName = '2lbk'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    models = inputParser.readModels(file)

    assert len(models) == 8
    assert models[0]['A'][0] == 'G'
    assert models[0]['A'][-1] == 'C'
    for model in models:
        assert len(model['A']) == 17


def test_strandTwoChains():
    structureName = '2z74'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    models = inputParser.readModels(file)

    assert len(models[0].keys()) == 2
    assert models[0]['A'][0] == 'A'
    assert models[0]['A'][-1] == 'U'
    assert models[0]['B'][0] == 'C'
    assert models[0]['B'][-1] == 'A'

def test_strandLensToughModel():
    structureName = '1ehz'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    models = inputParser.readModels(file)

    assert 1 == len(models[0].keys())
    assert 'G' == models[0]['A'][0]
    assert 'G' == models[0]['A'][9]
    for model in models:
        assert len(model['A']) == 76


def test_annotate():
    fileName = ROOT_DIR + '/downloadedStructures/pdb2lbk.ent'

    models = inputParser.annotate(fileName)

    assert len(models) == 8


def test_annotateCanonicalOnly():
    fileName = ROOT_DIR + '/downloadedStructures/pdb1ehz.ent'

    models = inputParser.annotate(fileName)

    assert 1 == len(models)
    assert 18 == len(models[0])


def test_properStrandRead():
    structureName = '1ehz'
    fileFormat = 'pdb'
    desiredOutput = list('gCGGAUUUAgCUCAGuuGGGAGAGCgCCAGAcUgAAgAucUGGAGgUCcUGUGuuCGaUCCACAGAAUUCGCACCA'.upper())

    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    strand = inputParser.readModels(file)[0]
    # txt = ''.join(strand['A'])
    # with open('tmp.txt', 'w') as file:
    #     file.write(txt)

    assert desiredOutput == strand['A']


def test_properLongStrandRead():
    structureName = '3g78'
    fileFormat = 'pdb'
    desiredOutput = list("GuGUGCCCGGCAUGGGUGCAGUCUAUAGGGUGAGAGUCCCGAACUGUGAAGGCAGAAGUA\
ACAGUUAGCCUAACGCAAGGGUGUCCGUGGCGACAUGGAAUCUGAAGGAAGCGGACGGCA\
AACCUUCGGUCUGAGGAACACGAACUUCAUAUGAGGCUAGGUAUCAAUGGAUGAGUUUGC\
AUAACAAAACAAAGUCCUUUCUGCCAAAGUUGGUACAGAGUAAAUGAAGCAGAUUGAUGA\
AGGGAAAGACUGCAUUCUUACCCGGGGAGGUCUGGAAACAGAAGUCAGCAGAAGUCAUAG\
UACCCUGUUCGCAGGGGAAGGACGGAACAAGUAUGGCGUUCGCGCCUAAGCUUGAACCGC\
CGUAUACCGAACGGUACGUACGGUGGUGUGAGAGGAGUUCGCUCUACUCUAU".upper())

    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    models = inputParser.readModels(file)

    txt = ''.join(models[0]['A'])
    base_pairs = inputParser.annotate(file)

    assert desiredOutput == models[0]['A']


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalOnlyTwoLayers():
    structureName = '1ehz'
    fileFormat = 'pdb'
    desiredOutput = "(((((((...(((.....[..)))..(((...........)))......((((..]....)))).)))))))...."

    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    models = inputParser.readModels(file)

    base_pairs = inputParser.annotate(file)

    dotNotation = inputParser.make_dot_notation(models[0], base_pairs[0])

    assert desiredOutput == dotNotation


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalOnlyOneLayer():
    structureName = '2lbk'
    fileFormat = 'pdb'
    desiredOutput = ".(((((.....)))))."

    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    strands = inputParser.readModels(file)

    models = inputParser.annotate(file)

    dotNotation = inputParser.make_dot_notation(strands[0], models[0])

    assert desiredOutput == dotNotation


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalHugeStructure():
    structureName = '3g78'
    fileFormat = 'pdb'
    desiredOutput = ".(((((.....)))))."
    desiredCount = '-.{[.(.(((<..(((((((((((...(((.......)))..(((((...{{{.{{{...\
                    )))))..(((...(((..((((.((((((....))))))))))...]>..)))...))).\
                    ..(((((((((((.(.....)...((((.......(......(...((((((..((((..\
                    [[[[[.))))...)))).}}}.}}}...))...)......)....)))))))))...)))\
                    )))...)))))))))))...}))).)(...((((....))))...).......(((.(..\
                    ..((..........))...))))....(.(((..(((......)))...))).).(((.(\
                    ((((((((....)))..)))))).)))...----------------------'.count('.')

    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    strands = inputParser.readModels(file)

    models = inputParser.annotate(file)
    strandmod = strands[0]
    strandmod.insert(0, '-')
    strandmod.extend(['-'] * 22)
    dotNotation = inputParser.make_dot_notation(strandmod, models[0])

    assert desiredCount == dotNotation.count('.')

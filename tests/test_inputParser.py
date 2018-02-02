#!/usr/bin/python3
import os

import pytest
from Bio.PDB.PDBParser import PDBParser

from app.helpers import inputParser
from defs import ROOT_DIR


def test_pdbListDownload():
    structureName = '1FAT'
    file = inputParser.online_input(structureName)
    assert structureName.lower() in file
    with open(file) as f:
        assert structureName in f.readline()


def test_pdbListDownloadWithFormat():
    structureName = '2LBK'
    fileFormat = 'pdb'
    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)
    assert structureName.lower() in file
    assert fileFormat in file
    with open(file) as f:
        assert structureName in f.readline()


def test_annotateDownloadedfile():
    structureName = '2lbk'
    fileFormat = 'pdb'
    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    models = inputParser.annotate(file)

    assert len(models) == 8


def test_strandLensMultipleModels():
    structureName = '2lbk'
    fileFormat = 'pdb'
    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    models = inputParser.read_models_from_pdb_file(file)

    assert len(models) == 8
    assert models[0]['A'][0] == 'G'
    assert models[0]['A'][-1] == 'C'
    for model in models:
        assert len(model['A']) == 17


def test_strandLensToughModel():
    structureName = '1ehz'
    fileFormat = 'pdb'
    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    models = inputParser.read_models_from_pdb_file(file)

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

    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    strand = inputParser.read_models_from_pdb_file(file)[0]
    # txt = ''.join(strand['A'])
    # with open('tmp.txt', 'w') as file:
    #     file.write(txt)

    assert desiredOutput == strand['A']


def test_properLongStrandRead():
    structureName = '3g78'
    fileFormat = 'pdb'
    desiredOutput = list("uGUGCCCGGCAUGGGUGCAGUCUAUAGGGUGAGAGUCCCGAACUGUGAAGGCAGAAGUA\
ACAGUUAGCCUAACGCAAGGGUGUCCGUGGCGACAUGGAAUCUGAAGGAAGCGGACGGCA\
AACCUUCGGUCUGAGGAACACGAACUUCAUAUGAGGCUAGGUAUCAAUGGAUGAGUUUGC\
AUAACAAAACAAAGUCCUUUCUGCCAAAGUUGGUACAGAGUAAAUGAAGCAGAUUGAUGA\
AGGGAAAGACUGCAUUCUUACCCGGGGAGGUCUGGAAACAGAAGUCAGCAGAAGUCAUAG\
UACCCUGUUCGCAGGGGAAGGACGGAACAAGUAUGGCGUUCGCGCCUAAGCUUGAACCGC\
CGUAUACCGAACGGUACGUACGGUGGUGUG".upper())

    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    models = inputParser.read_models_from_pdb_file(file)

    assert desiredOutput == models[0]['A']

def test_get_missing_residues_from_strandA_3G78():
    structure_name = '3g78'
    file_format = 'pdb'
    desiredCount = '-.{[.(.(((<..(((((((((((...(((.......)))..(((((...{{{.{{{...\
)))))..(((...(((..((((.((((((....))))))))))...]>..)))...))).\
..(((((((((((.(.....)...((((.......(......(...((((((..((((..\
[[[[[.))))...)))).}}}.}}}...))...)......)....)))))))))...)))\
)))...)))))))))))...}))).)(...((((....))))...).......(((.(..\
..((..........))...))))....(.(((..(((......)))...))).).(((.(\
((((((((....)))..)))))).)))...----------------------'.count('-')

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    missing_residues = inputParser.get_missing_residues_from_pdb(file)

    assert desiredCount == len(missing_residues)


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalOnlyTwoLayers():
    structureName = '1ehz'
    fileFormat = 'pdb'
    desiredOutput = "(((((((...(((.....[..)))..(((...........)))......((((..]....)))).)))))))...."

    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    strands = inputParser.read_models_from_pdb_file(file)

    models = inputParser.annotate(file)

    dotNotation = inputParser.make_dot_notation(strands[0], models[0])

    assert desiredOutput == dotNotation


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalOnlyOneLayer():
    structureName = '2lbk'
    fileFormat = 'pdb'
    desiredOutput = ".(((((.....)))))."

    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    strands = inputParser.read_models_from_pdb_file(file)

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

    file = inputParser.online_input(structure_name=structureName, file_format=fileFormat)

    strands = inputParser.read_models_from_pdb_file(file)

    models = inputParser.annotate(file)
    dotNotation = inputParser.make_dot_notation(strands[0]['A'], models[0])

    assert desiredCount == dotNotation.count('.')

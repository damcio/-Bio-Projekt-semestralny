#!/usr/bin/python3
import os

import pytest
from Bio.PDB.PDBParser import PDBParser

from app.helpers import inputParser
from defs import ROOT_DIR


def test_pdbListDownload():
    structure_name = '1FAT'
    file = inputParser.online_input(structure_name)
    assert structure_name.lower() in file
    with open(file) as f:
        assert structure_name in f.readline()


def test_pdbListDownloadWithFormat():
    structure_name = '2LBK'
    file_format = 'pdb'
    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)
    assert structure_name.lower() in file
    assert file_format in file
    with open(file) as f:
        assert structure_name in f.readline()


def test_annotateDownloadedfile():
    structure_name = '2lbk'
    file_format = 'pdb'
    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    models = inputParser.annotate(file)

    assert len(models) == 8


def test_strandLensMultipleModels():
    structure_name = '2lbk'
    file_format = 'pdb'
    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    models = inputParser.read_models_from_pdb_file(file)

    assert len(models) == 8
    assert models[0]['A'][0] == 'G'
    assert models[0]['A'][-1] == 'C'
    for model in models:
        assert len(model['A']) == 17


def test_strandTwoChains():
    structure_name = '2z74'
    file_format = 'pdb'
    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    models = inputParser.read_models_from_pdb_file(file)

    assert len(models[0].keys()) == 2
    assert models[0]['A'][0] == 'A'
    assert models[0]['A'][-1] == 'U'
    assert models[0]['B'][0] == 'G'
    assert models[0]['B'][-1] == 'A'


def test_strandLensToughModel():
    structure_name = '1ehz'
    file_format = 'pdb'
    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

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
    structure_name = '1ehz'
    file_format = 'pdb'
    desired_output = list('gCGGAUUUAgCUCAGuuGGGAGAGCgCCAGAcUgAAgAucUGGAGgUCcUGUGuuCGaUCCACAGAAUUCGCACCA'.upper())

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    strand = inputParser.read_models_from_pdb_file(file)[0]

    assert desired_output == strand['A']


def test_properLongStrandRead():
    structure_name = '3g78'
    file_format = 'pdb'
    desired_output = list("uGUGCCCGGCAUGGGUGCAGUCUAUAGGGUGAGAGUCCCGAACUGUGAAGGCAGAAGUA\
ACAGUUAGCCUAACGCAAGGGUGUCCGUGGCGACAUGGAAUCUGAAGGAAGCGGACGGCA\
AACCUUCGGUCUGAGGAACACGAACUUCAUAUGAGGCUAGGUAUCAAUGGAUGAGUUUGC\
AUAACAAAACAAAGUCCUUUCUGCCAAAGUUGGUACAGAGUAAAUGAAGCAGAUUGAUGA\
AGGGAAAGACUGCAUUCUUACCCGGGGAGGUCUGGAAACAGAAGUCAGCAGAAGUCAUAG\
UACCCUGUUCGCAGGGGAAGGACGGAACAAGUAUGGCGUUCGCGCCUAAGCUUGAACCGC\
CGUAUACCGAACGGUACGUACGGUGGUGUG".upper())

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    models = inputParser.read_models_from_pdb_file(file)

    assert desired_output == models[0]['A']


def test_get_missing_residues_from_strandA_3G78():
    structure_name = '3g78'
    file_format = 'pdb'
    desired_count = '-.{[.(.(((<..(((((((((((...(((.......)))..(((((...{{{.{{{...\
)))))..(((...(((..((((.((((((....))))))))))...]>..)))...))).\
..(((((((((((.(.....)...((((.......(......(...((((((..((((..\
[[[[[.))))...)))).}}}.}}}...))...)......)....)))))))))...)))\
)))...)))))))))))...}))).)(...((((....))))...).......(((.(..\
..((..........))...))))....(.(((..(((......)))...))).).(((.(\
((((((((....)))..)))))).)))...----------------------'.count('-')

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    missing_residues = inputParser.get_missing_residues_from_pdb(file)

    assert desired_count == len(missing_residues)


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalOnlyTwoLayers():
    structure_name = '1ehz'
    file_format = 'pdb'
    desired_output = "(((((((...(((.....[..)))..(((...........)))......((((..]....)))).)))))))...."

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    strands = inputParser.read_models_from_pdb_file(file)

    models = inputParser.annotate(file)

    dot_notation = inputParser.make_dot_notation(strands[0], models[0])

    assert desired_output == dot_notation


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalOnlyOneLayer():
    structure_name = '2lbk'
    file_format = 'pdb'
    desired_output = ".(((((.....)))))."

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    strands = inputParser.read_models_from_pdb_file(file)

    models = inputParser.annotate(file)

    dot_notation = inputParser.make_dot_notation(strands[0], models[0])

    assert desired_output == dot_notation


@pytest.mark.skip("might get dropped")
def test_buildDotNotationCanonicalHugeStructure():
    structure_name = '3g78'
    file_format = 'pdb'
    desired_count = '-.{[.(.(((<..(((((((((((...(((.......)))..(((((...{{{.{{{...\
)))))..(((...(((..((((.((((((....))))))))))...]>..)))...))).\
..(((((((((((.(.....)...((((.......(......(...((((((..((((..\
[[[[[.))))...)))).}}}.}}}...))...)......)....)))))))))...)))\
)))...)))))))))))...}))).)(...((((....))))...).......(((.(..\
..((..........))...))))....(.(((..(((......)))...))).).(((.(\
((((((((....)))..)))))).)))...----------------------'.count('.')

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    strands = inputParser.read_models_from_pdb_file(file)

    models = inputParser.annotate(file)
    dot_notation = inputParser.make_dot_notation(strands[0]['A'], models[0])

    assert desired_count == dot_notation.count('.')

#!/usr/bin/python3
import os

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

    strands = inputParser.readStrand(file)

    assert len(strands) == 8
    assert strands[0][0] == 'G'
    assert strands[0][-1] == 'C'
    for strand in strands:
        assert len(strand) == 17

def test_strandLensToughModel():
    structureName = '1ehz'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    strands = inputParser.readStrand(file)

    assert len(strands) == 1
    assert strands[0][0] == 'G'
    assert strands[0][9] == 'G'
    for strand in strands:
        assert len(strand) == 76

def test_annotate():
    fileName = ROOT_DIR+'/tests/testFiles/pdb2lbk.ent'

    models = inputParser.annotate(fileName)

    assert len(models) == 8
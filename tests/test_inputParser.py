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

def test_annotate():
    fileName = ROOT_DIR+'/tests/testFiles/pdb2lbk.ent'

    models = inputParser.annotate(fileName)

    assert len(models) == 8
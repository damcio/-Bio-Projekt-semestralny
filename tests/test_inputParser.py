#!/usr/bin/python3
from app.helpers import inputParser


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

def test_parseDownloadedfile():
    structureName = '2lbk'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    inputParser.localInput(file)


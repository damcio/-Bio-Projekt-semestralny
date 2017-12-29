#!/usr/bin/python3
from app.helpers import inputParser


def test_pdbListDownload():
    structureName = '1FAT'
    file = inputParser.onlineInput(structureName)
    assert structureName.lower() in file

def test_pdbListDownloadWithFormat():
    structureName = '2lbk'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)
    assert structureName in file
    assert fileFormat in file

def test_parseDownloadedfile():
    structureName = '2lbk'
    fileFormat = 'pdb'
    file = inputParser.onlineInput(structureName=structureName, fileFormat=fileFormat)

    inputParser.localInput(file)


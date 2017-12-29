#!/usr/bin/python3
import inputParser

def test_pdbListDownload():
    structureName = '1FAT'
    file = inputParser.onlineInput(structureName)
    assert 'cif' in file
    assert structureName.lower() in file

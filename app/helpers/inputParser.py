from Bio.PDB import PDBList, PDBParser, MMCIFParser

def onlineInput(structureName, fileFormat=None):
    pdbl = PDBList()

    return pdbl.retrieve_pdb_file(pdb_code=structureName, file_format=fileFormat)

def localInput(fileName, fileFormat=None):
    if 'mmcif' in fileName:
        parser = MMCIFParser()
    elif 'pdb' in fileName:
        parser = PDBParser()

    pass

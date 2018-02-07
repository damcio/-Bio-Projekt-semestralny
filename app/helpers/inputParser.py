import subprocess
from collections import OrderedDict, namedtuple

from defs import ROOT_DIR

import os
from Bio.PDB import PDBList, PDBParser

nono_words = ['MG', 'MN', 'HOH', 'K', 'G6P']
base_pair_tuple = namedtuple('base_pair', ['strand', 'position'])


def build_txt_strand_from_chains(model):
    """
    Concatenates chains in given model and builds a single strand from them

    :param model: Model aggregate of chains
    :return: RNA strand text string
    """
    txt = ''
    for _, chain in model.items():
        txt += ''.join(chain)

    return txt


def online_input(structure_name, file_format=None):
    """
    Uses BioPython's PDBList to download tertiary structure from PDB

    :param structure_name: PDB correct name of RNA structure
    :param file_format: File format to pull from database, currently only 'pdb' is supported
    :return:
    """
    if not file_format:
        file_format = 'pdb'
    pdbl = PDBList()

    return pdbl.retrieve_pdb_file(pdb_code=structure_name, file_format=file_format,
                                  pdir=ROOT_DIR + '/downloadedStructures/')


def read_base_pairs_from_mcannotate(annotated_output):
    """
    Parses MC-Annotate output for listed base pairs

    :param annotated_output: Output from MC-Annotate's analysis of tertiary structure
    :return: List of secondary structure's base pairs
    """
    out = iter(str(annotated_output).split('\n'))
    pairs = list(list())
    for line in out:
        tmp_list = list()
        if 'Base-pairs' in line.strip():
            tmp_line = next(out).strip().split()
            while tmp_line and 'Residue' not in tmp_line[0]:
                try:
                    if 'Ww/' in tmp_line[3] and tmp_line[6] == 'cis':
                        tmp_pair = tmp_line[0].split('-')
                        tmp_list.append((base_pair_tuple(tmp_pair[0][0], int(tmp_pair[0][1:])),
                                         base_pair_tuple(tmp_pair[1][0], int(tmp_pair[1][1:]))))
                except IndexError:
                    pass
                tmp_line = next(out).strip().split()
            pairs.append(tmp_list)
    return pairs


def make_offset_dict(strand):
    """
    Helper function for calculating offsets at which chains fall into in a model

    :param strand: Model to calculate offsets from, OrderedDict with chain names as keys
    :return: Dictionary of offsets
    """
    offset_dict = OrderedDict.fromkeys(strand.keys())
    temp_length = 0
    for chain, items in strand.items():
        offset_dict[chain] = temp_length
        temp_length += len(items)

    return offset_dict


def fix_base_pairs(strand, basepairs):
    """
    MC-Annotate provides base pairs between chains, with residue positions relative to chain,
    this function offsets positions so that pair positions are aligned to one concatenated chain

    :param strand: Model with chains to align
    :param basepairs: MC-Annotate's base pairs list
    :return: Offset base pairs list
    """
    offset_dict = make_offset_dict(strand)

    fixed_basepairs = []

    for pair in basepairs:
        fixed_basepairs.append(
            (
                base_pair_tuple(pair[0].strand, pair[0].position + offset_dict[pair[0].strand]),
                base_pair_tuple(pair[1].strand, pair[1].position + offset_dict[pair[1].strand])
            )
        )

    return fixed_basepairs


def get_missing_residues_from_pdb(file_path):
    """
    PDBParser does not include "missing" residues that are mentioned in PDB file,
    this function parses PDB file and returns mentioned missing residues

    :param file_path: Path to PDB file
    :return: List of missing residues, with information about which chain they come from
    """
    with open(file_path, 'r+') as file:
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

                return missing_residues
    return []


def read_models_from_pdb_file(file_path):
    """
    Parses PDB file for strands

    :param file_path: Path to PDB file
    :return: Returns list of models and chains they contain
    """
    pdbp = PDBParser()
    struct = pdbp.get_structure('file', file_path)
    models = list()

    for model in struct.get_models():
        strand = OrderedDict()
        for chain in model.get_chains():
            chain_name = chain.id
            strand.setdefault(chain_name, [])
            for res in chain.get_residues():
                temp_name = res.resname.lstrip()
                if temp_name not in nono_words:
                    if len(temp_name) > 1:
                        temp_name = temp_name[-1]
                    strand[chain_name].append(temp_name)
        models.append(strand)

    return models


def read_complete_models(file_path):
    """
    Appends residues missing from PDBParser to its output

    :param file_path: Path to PDB file
    :return: List of complete models
    """
    models = read_models_from_pdb_file(file_path)
    missing_residues = get_missing_residues_from_pdb(file_path)
    for model in models:
        for res in missing_residues:
            if res.chain in model.keys():
                model[res.chain].insert(res.index - 1,
                                        res.residue)
    return models

    return models


def trace_depth(pair, stacks, depth):
    """
    Helper function for tracing base pairs "depth" while building dot bracket notation string

    :param pair: Pair to analyze
    :param stacks: Current pair "stacks"
    :param depth: Current depth of main function
    :return: Depth at which make_dot_notation function should now operate
    """
    for stackPair in stacks[depth]:
        if pair[1].position > stackPair[1].position:
            if pair[0].position < stackPair[1].position:
                depth += 1
                return depth
    return depth


def make_dot_notation(strand, basepairs):
    """
    Makes Vienna-format text string secondary structure representation from
    MC-Annotate and PDBParser parsed input

    :param strand: Strand to build secondary structure of
    :param basepairs: List of base pairs found in the strand
    :return: Dot bracket notation text string
    """
    brackets = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]
    output = ['.' if i != '-' else '-' for i in strand]
    stacks = [[], [], [], [], []]
    depth = 0
    for pair in basepairs:
        temp_depth = depth
        depth = trace_depth(pair, stacks, depth)
        while temp_depth > depth:
            temp_depth = depth
            depth = trace_depth(pair, stacks, depth)

        stacks[depth].append(pair)
        depth = 0

    for i, stack in enumerate(stacks):
        for pair in stack:
            output[pair[0].position - 1], output[pair[1].position - 1] = brackets[i][0], brackets[i][1]

    return "".join(output)


def annotate_basepairs(file_path):
    """
    Runs external MC-Annotate program and returns base pairs found by it

    :param file_path: Path to PDB file
    :return: List of base pairs found by MC-Annotate
    """
    try:
        stdout = subprocess.run([ROOT_DIR + '/ext/MC_Annotate/MC-Annotate', file_path], stdout=subprocess.PIPE,
                                universal_newlines=True).stdout
    except FileNotFoundError:
        print("MC_Annotate not found! Please put the binary file in ./ext/MC_Annotate/")

    with open(ROOT_DIR + '/ext/MC_Annotate/output/' + file_path.split('/')[-1] + '.out', 'w+') as f:
        f.write(stdout)
    return read_base_pairs_from_mcannotate(stdout)

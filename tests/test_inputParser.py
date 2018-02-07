#!/usr/bin/python3
from collections import OrderedDict

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

    models = inputParser.annotate_basepairs(file)

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

    models = inputParser.annotate_basepairs(fileName)

    assert len(models) == 8


def test_annotateCanonicalOnly():
    fileName = ROOT_DIR + '/downloadedStructures/pdb1ehz.ent'

    models = inputParser.annotate_basepairs(fileName)

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


def test_make_offset_dict_three_chains():
    strand = OrderedDict.fromkeys('AZB')
    strand['A'] = ['C', 'G', 'A', 'G']
    strand['Z'] = ['U', 'A']
    strand['B'] = ['C', 'A', 'C']

    offset_dict = inputParser.make_offset_dict(strand)

    assert offset_dict['A'] == 0
    assert offset_dict['Z'] == 4
    assert offset_dict['B'] == 6


def test_fix_base_pairs_two_strand_model():
    structure_name = '3g78'
    file_format = 'pdb'

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)
    models = inputParser.read_complete_models(file)

    base_pairs = inputParser.annotate_basepairs(file)[0]

    fixed_base_pairs = inputParser.fix_base_pairs(models[0], base_pairs)

    pair_with_z_strand = [pair for pair in fixed_base_pairs if pair[1].strand == 'Z']

    assert pair_with_z_strand[0][1].position == 418


def test_build_dot_notation_canonical_only_two_layers():
    structure_name = '1ehz'
    file_format = 'pdb'
    desired_output = "(((((((...(((.....[..)))..(((...........)))......((((..]....)))).)))))))...."

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    models = inputParser.read_models_from_pdb_file(file)

    base_pairs = inputParser.annotate_basepairs(file)

    text_strand = inputParser.build_txt_strand_from_chains(models[0])

    dot_notation = inputParser.make_dot_notation(text_strand, base_pairs[0])

    assert desired_output == dot_notation


def test_build_dot_notation_canonical_only_one_layer():
    structure_name = '2lbk'
    file_format = 'pdb'
    desired_output = ".(((((.....)))))."

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    models = inputParser.read_models_from_pdb_file(file)
    text_strand = inputParser.build_txt_strand_from_chains(models[0])

    base_pairs = inputParser.annotate_basepairs(file)

    dot_notation = inputParser.make_dot_notation(text_strand, base_pairs[0])

    assert desired_output == dot_notation


def test_build_dot_notation_canonical_huge_structure():
    structure_name = '3g78'
    file_format = 'pdb'
    desired_count = '-.{[.(.(((<..(((((((((((...(((.......)))..(((((...{{{.{{{...\
)))))..(((...(((..((((.((((((....))))))))))...]>..)))...))).\
..(((((((((((.(.....)...((((.......(......(...((((((..((((..\
[[[[[.))))...)))).}}}.}}}...))...)......)....)))))))))...)))\
)))...)))))))))))...}))).)(...((((....))))...).......(((.(..\
..((..........))...))))....(.(((..(((......)))...))).).(((.(\
((((((((....)))..)))))).)))...----------------------.]]]]]...'.count('.')

    file = inputParser.online_input(structure_name=structure_name, file_format=file_format)

    models = inputParser.read_complete_models(file)
    text_strand = inputParser.build_txt_strand_from_chains(models[0])

    base_pairs = inputParser.annotate_basepairs(file)
    fixed_base_pairs = inputParser.fix_base_pairs(models[0], base_pairs[0])

    dot_notation = inputParser.make_dot_notation(text_strand, fixed_base_pairs)

    assert desired_count < dot_notation.count('.')

import argparse
import app.helpers.inputParser as inputParser
from app import nussinov

parser = argparse.ArgumentParser(description='Provide secondary-structure name')
parser.add_argument('structure_name', metavar='N', type=str, help='structure name from PDB')
parser.add_argument('prediction_mode', metavar='M', type=str,
                    help='mode for prediction: nuss - nussinov algorithm, mcan - mcannotate prediction')


def pretty_print_model(i, model, dot_notation):
    """
    Pretty prints secondary structure

    :rtype: void
    :param i: number of model to print
    :param model: text representation of RNA model to print
    :param dot_notation: dot_notation representation of RNA model to print
    """
    print("Model %d" % i)
    print(">strand_joined")
    for chunk in [model[i:i + 60] for i in range(0, len(model), 60)]:
        print(chunk)

    for chunk in [dot_notation[i:i + 60] for i in range(0, len(dot_notation), 60)]:
        print(chunk)


def nussinov_part(file_path):
    """
    Nussinov-algorithm branch of program

    :rtype: void
    :param file_path: path to pdb file
    """
    models = inputParser.read_complete_models(file_path)

    for i, model in enumerate(models):
        joined_chain_model = inputParser.build_txt_strand_from_chains(model)
        base_pair = nussinov.nussinov_algorithm(joined_chain_model)
        dot_notation = inputParser.make_dot_notation(joined_chain_model, base_pair)
        pretty_print_model(i, joined_chain_model, dot_notation)


def mcannotate_part(file_path):
    """

    :rtype: void
    :param file_path: path to pdb file
    """
    base_pairs = inputParser.annotate_basepairs(file_path)
    models = inputParser.read_complete_models(file_path)

    for i, (model, base_pair) in enumerate(zip(models, base_pairs)):
        fixed_base_pairs = inputParser.fix_base_pairs(model, base_pair)
        joined_chain_model = inputParser.build_txt_strand_from_chains(model)
        dot_notation = inputParser.make_dot_notation(joined_chain_model, fixed_base_pairs)
        pretty_print_model(i, joined_chain_model, dot_notation)


if __name__ == '__main__':
    args = parser.parse_args()
    prediction_mode = args.prediction_mode
    filepath = inputParser.online_input(args.structure_name, 'pdb')
    if prediction_mode == 'nuss':
        nussinov_part(filepath)
    elif prediction_mode == 'mcan':
        mcannotate_part(filepath)
    else:
        print("Unrecognized prediction mode!")

import argparse
import app.helpers.inputParser as inputParser

parser = argparse.ArgumentParser(description='Provide secondary-structure name')
parser.add_argument('structure_name', metavar='N', type=str, help='structure name from PDB')
parser.add_argument('prediction_mode', metavar='M', type=str, help='mode for prediction: nuss - nussinov algorithm, mcan - mcannotate prediction')


def pretty_print_model(model, dot_notation):
    pass

def nussinov_part(filename):
    pass


def mcannotate_part(filename):
    base_pairs = inputParser.annotate_basepairs(filename)
    models = inputParser.read_complete_models(filename)

    for i, (model, base_pair) in enumerate(zip(models, base_pairs)):
        fixed_base_pairs = inputParser.fix_base_pairs(model, base_pair)
        dot_notation = inputParser.make_dot_notation(inputParser.build_txt_strand_from_chains(model), fixed_base_pairs)
        print(dot_notation)


def main(args):
    prediction_mode = args.prediction_mode
    filename = inputParser.online_input(args.structure_name, 'pdb')
    if prediction_mode == 'nuss':
        nussinov_part(filename)
    elif prediction_mode == 'mcan':
        mcannotate_part(filename)
    else:
        print("Unrecognized prediction mode!")


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

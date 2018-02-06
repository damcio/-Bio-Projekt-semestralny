import argparse
import app.helpers.inputParser as inputParser

parser = argparse.ArgumentParser(description='Provide secondary-structure name')
parser.add_argument('structure_name', metavar='N', type=str, help='structure name from PDB')
parser.add_argument('prediction-mode', metavar='M', type=str, help='mode for prediction: nuss - nussinov algorithm, mcan - mcannotate prediction')


def main(args):
    filename = inputParser.online_input(args.structure_name, 'pdb')
    annotated = inputParser.annotate_basepairs(filename)
    inputParser.read_models_from_pdb_file(filename)
    print(annotated)


if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    main(args)

import argparse
import app.helpers.inputParser as inputParser

parser = argparse.ArgumentParser(description='Provide secondary-structure name')
parser.add_argument('structure_name', metavar='N', type=str, help='structure name from PDB')


def main(args):
    filename = inputParser.onlineInput(args.structure_name, 'pdb')
    annotated = inputParser.annotate(filename)
    inputParser.readModels(filename)
    print(annotated)


if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    main(args)

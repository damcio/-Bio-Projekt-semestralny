import argparse
import app.helpers.inputParser as inputParser

parser = argparse.ArgumentParser(description='Provide secondary-structure name')
parser.add_argument('structure_name', metavar='N', type=str, help='structure name from PDB')


def main(args):
    filename = inputParser.onlineInput(args.structure_name)
    a = 1
    pass

if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    main(args)

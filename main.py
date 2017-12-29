import argparse
import inputParser

parser= argparse.ArgumentParser(description='Provide secondary-structure name')
parser.add_argument('structure name', metavar='N', type=str, help='structure name from PDB')


def main():
    pass


if __name__ == '__main__':
    args = parser.parse_args()
    main()

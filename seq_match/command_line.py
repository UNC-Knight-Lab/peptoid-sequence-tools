from .seq_match import seq_match
import os
import argparse


def main():
    # Use current working directory as default input/output folder
    cwd = os.getcwd()
    input_folder = cwd
    output_folder = cwd

    # Parse arguments for input/output folder
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='Full path to input folder where mass spec text data is. ')
    parser.add_argument('-o', type=str, help='Full path to output folder where analysis file will be created.')
    args = parser.parse_args()

    if args.i:
        input_folder = args.i
    if args.o:
        output_folder = args.o

    # Call sequence matching function
    seq_match(input_folder, output_folder)

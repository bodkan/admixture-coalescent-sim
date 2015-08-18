"""Split the output file from a coalescent simulator containing genealogies of
several independent loci into multiple files (one file per locus).

Usage:
    reformat_tree_output.py --input=<FILENAME> --output_prefix=<STRING>
    reformat_tree_output.py -h | --help

Options:
    --input=<FILENAME>        input file with gene trees from a coalescent simulator
    --output_prefix=<STRING>  prefix of output files
    -h --help
"""

from docopt import docopt

def main(input_filename, output_prefix):
    try:
        with open(input_filename, "r") as tree_file:
            ith_locus = 0
            for line in tree_file:
                # skip the header lines
                if not (line.startswith("//") or line.startswith("[")): continue

                # // marks the beginning of a new locus -- when it's reached
                # open a new output file and close the old one
                if line.startswith("//"):
                    if ith_locus > 0: output_file.close()

                    ith_locus += 1
                    output_filename = "output/" + output_prefix + "_" + str(ith_locus) + ".txt"
                    output_file = open(output_filename, "w")

                    continue

                # dump the line into a separate output file removing [ and ]
                # characters from the beginning of the line
                # (the number in brackets specifies the length of a given
                # segment which is read by downstream R scripts and which gets
                # ignored unless the brackets are filtered out)
                print(line.replace("[", "").replace("]", ""), end = "", file = output_file)

    except FileNotFoundError:
        print("Input file {} not found!".format(input_filename))

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args["--input"], args["--output_prefix"])

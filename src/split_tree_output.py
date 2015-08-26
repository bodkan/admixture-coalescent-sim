"""Parse the file created by a coalescent simulator containing genealogies of
several independent loci and place the genealogies of a specified locus
into a separate output file.

Usage:
    reformat_tree_output.py --input=<FILENAME> --locus_number=<NUMBER> --output_prefix=<STRING>
    reformat_tree_output.py -h | --help

Options:
    --input=<FILENAME>        input file with gene trees from a coalescent simulator
    --locus_number=<NUMBER>   index of locus to pull out of the input file
    --output_prefix=<STRING>  prefix of he output file
    -h --help
"""

from docopt import docopt

def main(input_filename, locus_number, output_prefix):
    try:
        with open(input_filename, "r") as tree_file:
            ith_locus = 0
            for line in tree_file:
                # skip the header lines
                if not (line.startswith("//") or line.startswith("[")) or line == "\n": continue

                # // marks the beginning of a new locus -- when it's reached
                # open a new output file and close the old one
                elif line.startswith("//"):
                    ith_locus += 1
                    if ith_locus == locus_number:
                        output_filename = "output/" + output_prefix + str(locus_number) + ".newick"
                        output_file = open(output_filename, "w")
                    if ith_locus > locus_number:
                        output_file.close()
                        break

                # dump the line into a separate output file
                elif ith_locus == locus_number:
                    print(line, end = "", file = output_file)

    except FileNotFoundError:
        print("Input file {} not found!".format(input_filename))

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args["--input"], int(args["--locus_number"]), args["--output_prefix"])

import sys
import time
from optparse import OptionParser

formats = ['fastq', 'fasta', 'fa', 'fas']
fasta_files = []
separator = "###"


def check_format(name):
    for f in formats:
        if str(name).endswith(f):
            return True
    return False


def main(args=[]):
    usage = "usage: %prog [options] arg \nProgram 1. rename sequences"
    parser = OptionParser(usage, version='%prog version 1.0')

    parser.add_option("-f", "--fasta", dest="FASTA", help="fasta file with seqs", action="store", type="string")
    parser.add_option("-o", "--output", dest="OUTPUT", help="fasta file with seqs", action="store", type="string", default="re_named.fasta")

    (options, arg) = parser.parse_args(args)

    t_st = time.time()

    # read-write fasta file
    linie_counter = 0

    hl_out = open(options.OUTPUT, 'w')
    hl_in = open(options.FASTA, 'r')
    for line in hl_in.readlines():
        if str(line).startswith(">"):
            linie_counter += 1
            line = ">Contig_" + str(linie_counter) + "\n"

        hl_out.write(line)
    hl_out.close()


if __name__ == "__main__":
    main(sys.argv[1:])

import os
import sys
import time
from multiprocessing import Pool
from optparse import OptionParser

formats = ['fastq', 'fasta', 'fa', 'fas']


def alignmet_function(command_string):
    os.system(command_string)
    return command_string


def check_format(name):
    for f in formats:
        if str(name).endswith(f):
            return True
    return False


def main(args=[]):
    usage = "usage: %prog [options] arg \nProgram make alignment for all fasta files in directory"
    parser = OptionParser(usage, version='%prog version 1.0')

    parser.add_option("-f", "--folder", dest="FOLDER", help="folder", action="store", type="string")
    parser.add_option("-p", "--processor", dest="PROCESSOR", help="processors number", action="store", type="int", default=1)
    parser.add_option("-s", "--start_with", dest="STARTWITH", help="fasta file must start with this", action="store", type="string")
    parser.add_option("-o", "--out", dest="OUT", help="out directory", action="store", type="string", default="_alignment_out")

    (options, arg) = parser.parse_args(args)

    # Entering program    ##########################################

    t_st = time.time()

    if not os.path.isdir(options.FOLDER):
        sys.stdout.write('\nWrong input directory!')
        return

    # input/output things  ##########################################

    dir_for_alignment = options.OUT
    if not os.path.exists(dir_for_alignment):
        os.mkdir(dir_for_alignment)

    # alignment commands  ##########################################

    dirList = os.listdir(options.FOLDER)
    all_cmds_list = []
    fasta_list = []
    for fasta_file in dirList:
        if not check_format(fasta_file):
            continue

        if options.STARTWITH != "":
            if not str(fasta_file).startswith(options.STARTWITH):
                continue

        fasta_list.append(fasta_file)

    # dirList = os.listdir(options.FOLDER)
    for fasta_file in fasta_list:
        name_arr = fasta_file.split(".")

        # com_str = "einsi " + os.curdir + os.sep + options.FOLDER + os.sep + str(
        #     fasta_file) + " > " + os.curdir + os.sep + dir_for_alignment + os.sep + name_arr[0] + "_align.fasta"
        com_str = "mafft --auto " + os.curdir + os.sep + options.FOLDER + os.sep + str(
            fasta_file) + " > " + os.curdir + os.sep + dir_for_alignment + os.sep + name_arr[0] + "_align.fasta"

        all_cmds_list.append(com_str)

    # alignment   ##########################################

    pool = Pool(processes=int(options.PROCESSOR))  # start 4 worker processes
    pool.map(alignmet_function, all_cmds_list)  # prints "[0, 1, 4,..., 81]"

    # closing program    ##########################################

    t_end = time.time()
    sys.stdout.write("\nTime [s]: " + str(t_end - t_st))


if __name__ == "__main__":
    main(sys.argv[1:])

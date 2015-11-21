import os
import sys
import time
from Bio import SeqIO
from optparse import OptionParser

formats = ['fastq', 'fasta', 'fa', 'fas']


class AlignParseObject:
    def __init__(self, position):
        self.position = position
        self.nucl_dic = {"A": 0, "C": 0, "T": 0, "G": 0}

    def add_nucl(self, nucl):
        nucl = str(nucl).upper()
        if nucl in self.nucl_dic:
            self.nucl_dic[nucl] += 1

    def get_consensus(self):

        cons_list = []
        for k, v in self.nucl_dic.iteritems():
            cons_list.append((k, v))
        cons_list = sorted(cons_list, key=lambda length: length[1], reverse=True)

        if int(cons_list[0][1]) == 0:
            return ""
        return cons_list[0][0]


def alignmet_function(command_string):
    os.system(command_string)
    return command_string


def check_format(name):
    for f in formats:
        if str(name).endswith(f):
            return True
    return False


def main(args=[]):
    usage = "usage: %prog [options] arg \nProgram make consensus sequence from alignment and clean consensus from gaps"
    parser = OptionParser(usage, version='%prog version 1.0')

    parser.add_option("-f", "--folder", dest="FOLDER", help="folder", action="store", type="string")
    parser.add_option("-s", "--start_with", dest="STARTWITH", help="fasta file must start with this", action="store", type="string")
    parser.add_option("-o", "--out", dest="OUT", help="out directory", action="store", type="string", default="_alignment_consensus")

    (options, arg) = parser.parse_args(args)

    # Entering program    ##########################################

    t_st = time.time()

    if not os.path.isdir(options.FOLDER):
        sys.stdout.write('\nWrong input directory!')
        return

    log_info_hlr = open(options.FOLDER + os.sep + "outinfo.log", "w")

    # Input/output things  ##########################################

    dir_for_alignment = options.OUT
    if not os.path.exists(dir_for_alignment):
        os.mkdir(dir_for_alignment)

    # Alignment commands  ##########################################

    dirList = os.listdir(options.FOLDER)
    fasta_list = []

    for fasta_file in dirList:
        if not check_format(fasta_file):
            continue

        if options.STARTWITH != "":
            if not str(fasta_file).startswith(options.STARTWITH):
                continue

        fasta_list.append(fasta_file)

    counter = 0
    for fasta_file in fasta_list:
        counter += 1
        log_info_hlr.write("\n" + fasta_file)
        sys.stdout.write("\n" + fasta_file + " " + str(counter) + " of " + str(len(fasta_list)))

        seq_new_name = ""
        dic_positionInSeq_AlignParseObjects = {}
        for seq in SeqIO.parse(open(os.curdir + os.sep + options.FOLDER + os.sep + fasta_file, 'r'), "fasta"):  # SeqIO, AlignIO
            s = str(seq.seq)
            seq_new_name = seq.name

            for i in xrange(len(s)):
                if i not in dic_positionInSeq_AlignParseObjects:
                    ob = AlignParseObject(i)
                    dic_positionInSeq_AlignParseObjects[i] = ob
                ob = dic_positionInSeq_AlignParseObjects[i]
                ob.add_nucl(s[i])

        seq = ""
        for k, v in dic_positionInSeq_AlignParseObjects.iteritems():
            seq += v.get_consensus();

        # generating results
        cluster_file_name = str(fasta_file)[0:str(fasta_file).index(".fasta")] + "_" + "cons" + ".fasta"
        file_hlr_cluster_seqs = open(dir_for_alignment + os.sep + cluster_file_name, "w")
        file_hlr_cluster_seqs.write(">" + seq_new_name + "\n")
        file_hlr_cluster_seqs.write(seq + "\n")

        file_hlr_cluster_seqs.close()

    # Closing program    ##########################################

    t_end = time.time()
    sys.stdout.write("\nTime [s]: " + str(t_end - t_st))

    log_info_hlr.close()


if __name__ == "__main__":
    main(sys.argv[1:])

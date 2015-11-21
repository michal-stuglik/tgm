import os
import sys
import time
from Bio import SeqIO
from optparse import OptionParser

formats = ['fastq', 'fasta', 'fa', 'fas']
seqname_cluster_number = {}


def lookup_cluster(clusters_num_seqList_dic={}):
    dic = {}  # return {sequencename:list_of_culsters}
    for cluster_num, sequence_list in clusters_num_seqList_dic.iteritems():
        for seqname in sequence_list:
            if seqname not in dic:
                dic[seqname] = []
            l = dic[seqname]
            l.append(cluster_num)
            dic[seqname] = l
    return dic


def make_first_cluster(sequence_culstering_dic={}):
    cluster_num_seqList_dic = {}
    c = 0
    for k, v in sequence_culstering_dic.iteritems():
        c += 1
        l = []
        l.append(k)
        for s in v:
            if s not in l:
                l.append(s)
        cluster_num_seqList_dic[c] = l

    return cluster_num_seqList_dic


def merge_clusters(lookup_dic={}, clusters_num_seqList_dic={}):
    for seqname, list_of_clusters in lookup_dic.items():
        if len(list_of_clusters) < 2: continue;

        list_to_aggreg_segs_in_clusters = []

        for c in list_of_clusters:
            if c not in clusters_num_seqList_dic: continue;
            l = clusters_num_seqList_dic[c]
            for s in l:
                if s not in list_to_aggreg_segs_in_clusters: list_to_aggreg_segs_in_clusters.append(s)
            del clusters_num_seqList_dic[c]
        if len(list_to_aggreg_segs_in_clusters) > 0: clusters_num_seqList_dic[c] = list_to_aggreg_segs_in_clusters;

    return clusters_num_seqList_dic


def are_there_more_clusters_for_one_seq(lookup_dic={}):
    are_there = False
    for seqname, cluster_list in lookup_dic.iteritems():
        if len(cluster_list) > 1:
            are_there = True
            break
    return are_there


def do_all_clustering(sequence_culstering_dic={}):
    clusters_num_seqList_dic = make_first_cluster(sequence_culstering_dic)
    lookup_dic = lookup_cluster(clusters_num_seqList_dic)

    bool_do_more_mergeing = are_there_more_clusters_for_one_seq(lookup_dic)

    while bool_do_more_mergeing:
        clusters_num_seqList_dic = merge_clusters(lookup_dic, clusters_num_seqList_dic)

        lookup_dic = lookup_cluster(clusters_num_seqList_dic)
        bool_do_more_mergeing = are_there_more_clusters_for_one_seq(lookup_dic)

    return clusters_num_seqList_dic


def generate_clustering_dictionatry(IDENTITY, dic_of_seqs_in_file={}):
    sequence_culstering_dic = {}

    used_sequenceNames_list = []
    for s1 in dic_of_seqs_in_file.values():
        for s2 in dic_of_seqs_in_file.values():
            if s1.name == s2.name: continue;

            seq_list_1 = list(str(s1.seq))
            seq_list_2 = list(str(s2.seq))

            ident = 0
            non_ident = 0

            for i in xrange(len(seq_list_1)):

                try:
                    if seq_list_1[i] == "-": continue;
                    if seq_list_2[i] == "-": continue;

                    if seq_list_1[i] == seq_list_2[i]:
                        ident += 1
                    elif seq_list_1[i] != seq_list_2[i]:
                        non_ident += 1
                except:
                    print "Unexpected error:", sys.exc_info()[0]
                    print "\n" + str(len(seq_list_1))
                    print "\n" + s1.name
                    print "\n" + s2.name
                    break

            if ident + non_ident == 0: continue  # non-overlaping seqs

            overall_ident = float(ident) / float(ident + non_ident) * 100

            if overall_ident >= float(IDENTITY):

                # s1
                if s1.name not in sequence_culstering_dic:
                    l = []
                else:
                    l = sequence_culstering_dic[s1.name]
                l.append(s2.name)
                sequence_culstering_dic[s1.name] = l
                used_sequenceNames_list.append(s1.name)

                # s2
                if s2.name not in sequence_culstering_dic:
                    l = []
                else:
                    l = sequence_culstering_dic[s2.name]
                l.append(s1.name)
                sequence_culstering_dic[s2.name] = l
                used_sequenceNames_list.append(s2.name)

        # new clusters for diverged sequnece
        for s1 in dic_of_seqs_in_file.values():
            if s1.name not in used_sequenceNames_list:
                l = [s1.name]
                sequence_culstering_dic[s1.name] = l

    return sequence_culstering_dic


def generate_cluster_files(fasta_file, dir_for_alignment, assesed_clusters_dic={}):
    for c_num, c_list in assesed_clusters_dic.iteritems():
        cluster_file_name = str(fasta_file)[0:str(fasta_file).index(".fasta")] + "_" + str(c_num) + ".fasta"
        file_hlr_cluster_seqs = open(dir_for_alignment + os.sep + cluster_file_name, "w")
        SeqIO.write(c_list, file_hlr_cluster_seqs, "fasta")
        file_hlr_cluster_seqs.close()


def alignmet_function(command_string):
    os.system(command_string)
    return command_string


def check_format(name):
    for f in formats:
        if str(name).endswith(f):
            return True
    return False


def main(args=[]):
    usage = "usage: %prog [options] arg \nProgram make assess alignment for all fasta files in directory, if this is nessesery it split alignment on more than one file"
    parser = OptionParser(usage, version='%prog version 1.0')

    parser.add_option("-f", "--folder", dest="FOLDER", help="folder", action="store", type="string")
    parser.add_option("-p", "--processor", dest="PROCESSOR", help="processors number", action="store", type="int", default=1)
    parser.add_option("-s", "--start_with", dest="STARTWITH", help="fasta file must start with this", action="store", type="string")
    parser.add_option("-i", "--identity", dest="IDENTITY", help="level of identity", action="store", type="int", default=99)
    parser.add_option("-o", "--out", dest="OUT", help="out directory", action="store", type="string", default="_alignment_assesment_out")

    (options, arg) = parser.parse_args(args)

    # Entering program    ##########################################

    t_st = time.time()

    if not os.path.isdir(options.FOLDER):
        sys.stdout.write('\nWrong input directory!')
        return

    log_info_hlr = open(options.FOLDER + os.sep + "outoutinfo.log", "w")

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

    for fasta_file in fasta_list:

        log_info_hlr.write(str(fasta_file))
        sys.stdout.write("\n" + str(fasta_file))

        dic_of_seqs_in_file = {}
        for seq in SeqIO.parse(open(os.curdir + os.sep + options.FOLDER + os.sep + fasta_file, 'r'), "fasta"):  # SeqIO, AlignIO
            dic_of_seqs_in_file[seq.name] = seq

        log_info_hlr.write("\ngenerating clustering dictionary: " + str(fasta_file))
        sequences_culstering_dic = generate_clustering_dictionatry(options.IDENTITY, dic_of_seqs_in_file)

        clusters_num_seqList_dic = do_all_clustering(sequences_culstering_dic)

        cluster_numerator = 0
        for cluster_number, sequenceNameList in clusters_num_seqList_dic.iteritems():

            sequence_to_file_list = []

            print "\n" + fasta_file + " sequenceNameList count: " + str(len(sequenceNameList))
            for s in sequenceNameList:
                seq = dic_of_seqs_in_file[s]
                sequence_to_file_list.append(seq)

            # results
            cluster_numerator += 1
            generate_cluster_files(fasta_file, str(os.curdir + os.sep + dir_for_alignment), {cluster_numerator: sequence_to_file_list})

    # Closing program    ##########################################

    t_end = time.time()
    sys.stdout.write("\nTime [s]: " + str(t_end - t_st))
    log_info_hlr.close()


if __name__ == "__main__":
    main(sys.argv[1:])

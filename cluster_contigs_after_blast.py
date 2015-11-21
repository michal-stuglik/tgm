import os
import sys
import time
from Bio import SeqIO
from optparse import OptionParser

req_length = 1
merged_db = 'merged_db.fasta'
formats = ['fastq', 'fasta', 'fa', 'fas']
fasta_files = []
separator = "###"

relations_dic = {}
seq_name_order_dic = {}
ordered_list = []


class AlignmentPairObject:
    def __init__(self, name, align_length, shorter_length, reverse):
        self.name = name
        self.align_length = align_length
        self.shorter_length = shorter_length
        self.reverse = reverse


def get_true_stranded_seq_list(non_ordered_seq_list=[]):
    global ordered_list
    global seq_name_order_dic

    ordered_list = []
    seq_name_order_dic = {}

    to_sort_in_future_list = set_seqs_in_order(non_ordered_seq_list, seq_name_order_dic)

    while len(to_sort_in_future_list) > 0:
        to_sort_in_future_list = set_seqs_in_order(to_sort_in_future_list, seq_name_order_dic)

    return ordered_list


def set_seqs_in_order(non_ordered_seq_list=[], seq_name_order_dic={}):
    to_sort_in_future_list = []

    for seq in non_ordered_seq_list:

        if len(seq_name_order_dic) == 0:
            ordered_list.append(seq)
            seq_name_order_dic[seq.name] = 0
        else:

            added_to_list = False
            for seq_name_in_list, rev_counter in seq_name_order_dic.iteritems():

                if seq_name_in_list in relations_dic[seq.name]:
                    do_rc = relations_dic[seq.name][seq_name_in_list]

                    # relations analysis
                    if do_rc:
                        if seq_name_order_dic[seq_name_in_list] == 0:
                            rc_seq = seq.reverse_complement(id=seq.id, name=True, description=True)
                            ordered_list.append(rc_seq)
                            seq_name_order_dic[rc_seq.name] = 1
                        elif seq_name_order_dic[seq_name_in_list] == 1:
                            ordered_list.append(seq)
                            seq_name_order_dic[seq.name] = 0
                    else:
                        if seq_name_order_dic[seq_name_in_list] == 0:
                            ordered_list.append(seq)
                            seq_name_order_dic[seq.name] = 0
                        elif seq_name_order_dic[seq_name_in_list] == 1:
                            rc_seq = seq.reverse_complement(id=seq.id, name=True, description=True)
                            ordered_list.append(rc_seq)
                            seq_name_order_dic[rc_seq.name] = 1

                    added_to_list = True
                    break

            if not added_to_list:
                to_sort_in_future_list.append(seq)

    return to_sort_in_future_list


def set_relations_dic(ap, c1, c2):
    global relations_dic
    if c1 not in relations_dic:
        d = {c2: ap.reverse}
        relations_dic[c1] = d
    else:
        d = relations_dic[c1]
        if c2 not in d:
            d[c2] = ap.reverse
            relations_dic[c1] = d


def sort_seq_in_list(sequences_list=[]):
    new_list = []
    for seq in sequences_list:
        new_list.append((seq.name, len(str(seq.seq))))

    new_list = sorted(new_list, key=lambda length: length[1], reverse=True)

    sorted_list = []
    for t in new_list:
        sorted_list.append(t[0])

    return__seq_list = []
    for seq_name in sorted_list:
        for seq in sequences_list:
            if seq_name == seq.name:
                return__seq_list.append(seq)
                break

    return return__seq_list


def merge_contigs(contigs_lookup={}, contigs_clusters={}):
    # {seq_name:dic}
    for contig_name, dic_of_cluster_num in contigs_lookup.items():

        # I am interested in contigs in more than 1 clusters
        if len(dic_of_cluster_num) < 2:
            continue

        n_c_dics = {}

        for cluster_number in dic_of_cluster_num.keys():
            for contig_name_2 in contigs_clusters[cluster_number]:
                if contig_name_2 not in n_c_dics:
                    n_c_dics[contig_name_2] = 0
            del contigs_clusters[cluster_number]

        contigs_clusters[cluster_number] = n_c_dics
        break

    return contigs_clusters


def get_contigs_list_as_dic(contigs_clusters={}):
    contigs_dics = {}
    for n, contigs in contigs_clusters.iteritems():
        for c in contigs:
            if c not in contigs_dics:
                contigs_dics[c] = 0

    return contigs_dics


def contigs_lookup(contigs_dics={}, contigs_clusters={}):
    lookup_dic = {}
    lookup_dic = lookup_dic.fromkeys(contigs_dics.keys(), None)

    for n, contigs in contigs_clusters.iteritems():
        for c in contigs:
            l = lookup_dic[c]
            if l == None:
                l = {}
            l[n] = None
            lookup_dic[c] = l

    seq_contigs_id = {}
    for seq, dic in lookup_dic.iteritems():
        if len(dic) > 1:
            seq_contigs_id[seq] = dic

    return seq_contigs_id


def check_format(name):
    for f in formats:
        if str(name).endswith(f):
            return True
    return False


def check_reqs(min_identity, identity, min_evalue, evalue):
    if identity < min_identity:
        return False
    if evalue > min_evalue:
        return False
    return True


def main(args=[]):
    usage = "usage: %prog [options] arg \nProgram 1. select best ortolog seq for  blast query, 2. "
    parser = OptionParser(usage, version='%prog version 1.0')

    parser.add_option("-b", "--blast", dest="BLASTRESULT_FILENAME", help="Txt file from blast", metavar="FILE")
    parser.add_option("-f", "--query_fasta", dest="QUERY_FASTA", help="fasta file with seqs for query in bast", action="store", type="string")
    parser.add_option("-d", "--database_fasta", dest="DATABASE_FASTA", help="fasta file with seqs for subject in bast", action="store", type="string")
    parser.add_option("-o", "--output_folder", dest="OUTPUT_FOLDER", help="output folder")
    parser.add_option("-v", "--identity_value", dest="IDENTITY_VALUE", help="identity value", type="float", action='store', default=0.00)
    parser.add_option("-k", "--align2shortest_value", dest="ALIGN2SHORTEST_VALUE", help="identity value", type="float", action='store', default=0.00)

    parser.add_option("-q", "--query", dest="QUERY_COLUMN", help="query column", type="int", action='store', default=1)
    parser.add_option("-s", "--subject", dest="SUBJECT_COLUMN", help="subject column", type='int', action='store', default=2)
    parser.add_option("-e", "--evalue", dest="EVALUE_COLUMN", help="e-value column", type="int", action='store', default=11)
    parser.add_option("-i", "--identity", dest="IDENTITY_COLUMN", help="identity column", type="int", action='store', default=3)
    parser.add_option("-a", "--alignment", dest="ALIGNMENT_COLUMN", help="alignment column", type="int", action='store', default=4)

    parser.add_option("-g", "--query_start", dest="QUERY_START_COLUMN", help="query start column", type="int", action='store', default=7)
    parser.add_option("-z", "--query_end", dest="QUERY_END_COLUMN", help="query end column", type='int', action='store', default=8)

    parser.add_option("-j", "--subject_start", dest="SUBJECT_START_COLUMN", help="subject start column", type="int", action='store', default=9)
    parser.add_option("-l", "--subject_end", dest="SUBJECT_END_COLUMN", help="subject end column", type='int', action='store', default=10)

    parser.add_option("-p", "--progress", dest="BOOL_PROGRESS", help="show progress", action='store', default='True')

    (options, arg) = parser.parse_args(args)

    # Entering program    ##########################################

    t_st = time.time()

    if not os.path.isdir(options.OUTPUT_FOLDER):
        sys.stdout.write('\nWrong output directory!')
        return

    log_info_hlr = open(options.OUTPUT_FOLDER + os.sep + "outinfo.log", "w")

    # Sequence length dictionary   ##########################################

    db_seq_length_dic = {}
    for seq in SeqIO.parse(open(options.QUERY_FASTA, 'r'), "fasta"):
        db_seq_length_dic[seq.name] = len(str(seq.seq))

    # Parsing blast file   ##########################################

    linie_counter = 0
    pair_counter = 0
    pairs_dic = {}

    hl_in = open(options.BLASTRESULT_FILENAME, 'r')
    for line in hl_in.readlines():
        linie_counter += 1
        if str(line).startswith("#"):
            continue
        try:
            query = str(line).split('\t')[options.QUERY_COLUMN - 1]
            if len(query.split(' ')) > 1:
                query = query.split(' ')[0]
            subject = str(line).split('\t')[options.SUBJECT_COLUMN - 1]
            identity = float((line).split('\t')[options.IDENTITY_COLUMN - 1])
            alignment_length = int((line).split('\t')[options.ALIGNMENT_COLUMN - 1])

            subject_start = int((line).split('\t')[options.SUBJECT_START_COLUMN - 1])
            subject_end = int((line).split('\t')[options.SUBJECT_END_COLUMN - 1])

            do_reverse = subject_start > subject_end

            if query not in db_seq_length_dic: continue
            if subject not in db_seq_length_dic: continue

            query_tupla = (
                db_seq_length_dic[query],
                int(line.split('\t')[options.QUERY_START_COLUMN - 1]),
                int(line.split('\t')[options.QUERY_END_COLUMN - 1])
            )
            subject_tupla = (
                db_seq_length_dic[subject],
                int(line.split('\t')[options.SUBJECT_START_COLUMN - 1]),
                int(line.split('\t')[options.SUBJECT_END_COLUMN - 1])
            )

            length_list_to_check = [query_tupla, subject_tupla]
            length_list_to_check = sorted(length_list_to_check, key=lambda length: length[0], reverse=False)

            # identity:
            if options.IDENTITY_VALUE > 0:
                if identity < float(options.IDENTITY_VALUE) * 100.0:
                    continue

            if query == subject:
                continue

            pair_counter += 1

            str_pair = query + "\t" + subject
            str_pair_reverse = subject + "\t" + query

            alignment_pair = AlignmentPairObject(str_pair, alignment_length, length_list_to_check[0][0], do_reverse)

            if str_pair_reverse not in pairs_dic:
                if str_pair not in pairs_dic:
                    pairs_dic[str_pair] = alignment_pair
                else:
                    ap = pairs_dic[str_pair]
                    ap.align_length += alignment_pair.align_length
                    pairs_dic[str_pair] = ap

            if linie_counter % 10000 == 0:
                sys.stdout.write("\nlines: " + str(linie_counter))
        except:
            print "Unexpected error:", os.sys.exc_info()[0]
            raise
    hl_in.close()

    sys.stdout.write("\nunique pairs_dic: " + str(pair_counter))

    # Pairs filtering   ##########################################

    pairs_black_list = []
    for k, v in pairs_dic.iteritems():
        align_frac_test = float(v.align_length) / float(v.shorter_length)

        if options.ALIGN2SHORTEST_VALUE > align_frac_test:
            pairs_black_list.append(k)

    for p in pairs_black_list:
        del pairs_dic[p]

    # Relations tree   ##########################################

    for s_p, ap in pairs_dic.iteritems():
        q = s_p.split('\t')[0]
        s = s_p.split('\t')[1]

        set_relations_dic(ap, q, s)
        set_relations_dic(ap, s, q)

    # Clustering   ##########################################

    cluster_numerator = 0
    contigs_clusters = {}  # number - list fo contigs

    pairs_counter = 0
    no_pairs = len(pairs_dic)
    for s, n in pairs_dic.iteritems():

        pairs_counter += 1
        if pairs_counter % 1000 == 0:
            sys.stdout.write("\n" + str(pairs_counter) + " of " + str(no_pairs))

        s_arrary = s.split("\t")

        s1 = s_arrary[0]
        s2 = s_arrary[1]

        if len(contigs_clusters) == 0:
            cluster_numerator += 1
            contigs_clusters[cluster_numerator] = {s1: 0, s2: 0}
            continue

        create_new_cluster_for_contigs = True
        add_contigs_to_existing_cluster = False
        contigs_cluster_number_in_dict = -1

        # check if contig is already in some cluster
        # if not, create new one
        for n, c_list in contigs_clusters.iteritems():

            if s1 not in c_list and s2 not in c_list:
                continue

            if s1 in c_list:
                add_contigs_to_existing_cluster = True
                contigs_cluster_number_in_dict = n
                create_new_cluster_for_contigs = False
                break
            if s2 in c_list:
                add_contigs_to_existing_cluster = True
                contigs_cluster_number_in_dict = n
                create_new_cluster_for_contigs = False
                break

        if add_contigs_to_existing_cluster:
            cluster = contigs_clusters[contigs_cluster_number_in_dict]

            if s1 not in cluster:
                cluster[s1] = 0

            if s2 not in cluster:
                cluster[s2] = 0

        if create_new_cluster_for_contigs:
            cluster_numerator += 1
            contigs_clusters[cluster_numerator] = {s1: 0, s2: 0}

    sys.stdout.write("\nNo.clusters: " + str(len(contigs_clusters)))

    # Merging clustering   ##########################################

    contigs_dics = get_contigs_list_as_dic(contigs_clusters)
    m_contigs_lookup = contigs_lookup(contigs_dics, contigs_clusters)

    merge_counter = 0
    sys.stdout.write("\ncontigs before merging: " + str(len(contigs_clusters)))

    while len(m_contigs_lookup) > 0:
        merge_counter += 1
        sys.stdout.write("\nmerging: " + str(merge_counter) + ", clusters to merge: " + str(len(m_contigs_lookup)))

        contigs_clusters = merge_contigs(m_contigs_lookup, contigs_clusters)
        m_contigs_lookup = contigs_lookup(contigs_dics, contigs_clusters)

    sys.stdout.write("\ncontigs after merging: " + str(len(contigs_clusters)))

    # Logging     ##########################################

    log_info_hlr.write("\nLogowanie zawartosci pairs_dic:")
    for k, v in pairs_dic.iteritems():
        log_info_hlr.write("\n" + str(k) + "\t" + str(v.align_length) + "\t" + str(v.shorter_length) + "\t" + str(float(v.align_length) / float(v.shorter_length)))

    # Create multi-fasta files     ##########################################

    # multifasta for each cluster:
    sys.stdout.write("\nGenerating cluster seq in fasta...")

    seqName_seq_dic = {}
    for seq in SeqIO.parse(open(options.QUERY_FASTA, 'r'), "fasta"):
        seqName_seq_dic[seq.name] = seq

    # generating cluster files::
    contigs_used_in_all_clusters = []
    for k, v in contigs_clusters.iteritems():

        sequences_list = []
        for c in v:
            sequences_list.append(seqName_seq_dic[c])
            contigs_used_in_all_clusters.append(c)

        # sorting
        sequences_list = sort_seq_in_list(sequences_list)

        # ordered
        ordered_seq_list = get_true_stranded_seq_list(sequences_list)

        cluster_file_name = "cluster_" + str(k) + ".fasta"
        file_hlr_cluster_seqs = open(options.OUTPUT_FOLDER + os.sep + cluster_file_name, "w")
        SeqIO.write(ordered_seq_list, file_hlr_cluster_seqs, "fasta")
        file_hlr_cluster_seqs.close()

        # logging
        log_info_hlr.write("\nCluster [name & count]: " + cluster_file_name + "\t" + str(len(ordered_seq_list)))
        for s in ordered_seq_list:
            log_info_hlr.write("\n" + s.name)

    # generating fasta file for contigs with hit only to itself (non cluster)
    sequences_list = []
    file_hlr_all_non_cluster_seqs = open(options.OUTPUT_FOLDER + os.sep + "all_non_cluster.fasta", "w")
    for name, seq in seqName_seq_dic.iteritems():
        if name not in contigs_used_in_all_clusters:
            sequences_list.append(seq)

    SeqIO.write(sequences_list, file_hlr_all_non_cluster_seqs, "fasta")
    file_hlr_all_non_cluster_seqs.close()

    # closing program    ##########################################
    t_end = time.time()
    sys.stdout.write("\nTime [s]: " + str(t_end - t_st))
    log_info_hlr.close()


if __name__ == "__main__":
    main(sys.argv[1:])

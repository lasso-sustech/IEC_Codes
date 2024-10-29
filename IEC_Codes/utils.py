# This file contains the utilities used for this project.
import copy
import numpy as np
from tqdm import tqdm
import pickle
import time
import Levenshtein as Lvst

map_list = {'A': 0, 'T': 1, 'G': 2, 'C': 3}


# Map numbers to bases
def Num2Base(Num):
    Base = []
    for n in Num:
        Base.append(Num2Base_single(n))
    return Base


# Map the number to the corresponding base
def Num2Base_single(n):
    if n == 0:
        return 'A'
    elif n == 1:
        return 'T'
    elif n == 2:
        return 'G'
    elif n == 3:
        return 'C'


# Divide DNA sequences into two categories according to whether the length is correct or not
def seqLenClassify(DNAseqfile, rightlenfile, wronglenfile, l_seq, deltal):
    f_DNA = open(DNAseqfile, 'r')
    f_rightlen = open(rightlenfile, 'w')
    f_wronglen = open(wronglenfile, 'w')
    seq_num = -1

    rseqind_list = []
    wseqind_list = []

    seq_list = f_DNA.readlines()

    for i in tqdm(seq_list):
        seq_num += 1
        seq = seq_list[seq_num]
        seq_ind = seq[0: seq.find(' ')]
        seq_str = seq[(seq.find(' ') + 1): (seq.find('\t'))]
        if len(seq_str) == l_seq:

            f_rightlen.write(seq)
            rseqind_list.append(int(seq_ind))

        elif abs(len(seq_str) - l_seq) <= deltal:

            f_wronglen.write(seq)
            wseqind_list.append(int(seq_ind))

    f_DNA.close()
    f_rightlen.close()
    f_wronglen.close()
    return seq_list, rseqind_list, wseqind_list


# Sequence clustering based on Levenshtein distance
def LvstClustering(seq_list, rseqind_list, wseqind_list, l_int, maxdis):
    clusterind = []

    time.sleep(0.1)
    for i in tqdm(rseqind_list):
        seq = seq_list[i - 1]
        seq_ind = int(seq[0: seq.find(' ')])
        seq_str = seq[seq.find(' ') + 1: seq.find('\t')]
        seq_head = seq_str[:l_int]
        seq_tail = seq_str[-l_int:]
        flag = 0
        for j in range(len(clusterind)):
            seq_compare_ind = clusterind[j][0]
            seq_compare = seq_list[seq_compare_ind - 1]
            seq_compare_str = seq_compare[seq_compare.find(' ') + 1: seq_compare.find('\t')]
            seq_compare_head = seq_compare_str[:l_int]

            Lvst_dis_head = Lvst.distance(seq_compare_head, seq_head)
            if Lvst_dis_head <= maxdis:
                seq_compare_tail = seq_compare_str[-l_int:]
                Lvst_dis_tail = Lvst.distance(seq_compare_tail, seq_tail)
                if Lvst_dis_tail <= maxdis:
                    flag = 1
                    clusterind[j].append(seq_ind)
                    break
        if not flag:
            clusterind.append([seq_ind])
    for i in tqdm(wseqind_list):
        seq = seq_list[i - 1]
        seq_ind = int(seq[0: seq.find(' ')])
        seq_str = seq[seq.find(' ') + 1: seq.find('\t')]
        seq_head = seq_str[:l_int]
        seq_tail = seq_str[-l_int:]
        for j in range(len(clusterind)):
            seq_compare_ind = clusterind[j][0]
            seq_compare = seq_list[seq_compare_ind - 1]
            seq_compare_str = seq_compare[seq_compare.find(' ') + 1: seq_compare.find('\t')]
            seq_compare_head = seq_compare_str[:l_int]

            Lvst_dis_head = Lvst.distance(seq_compare_head, seq_head)
            if Lvst_dis_head <= maxdis:
                seq_compare_tail = seq_compare_str[-l_int:]
                Lvst_dis_tail = Lvst.distance(seq_compare_tail, seq_tail)
                if Lvst_dis_tail <= maxdis:
                    clusterind[j].append(seq_ind)
                    break
    return clusterind


# Obtain the common subsequences between two sequences
def getCommonSubseq(seq1, seq2, corr_len):
    lseq1 = len(seq1)
    lseq2 = len(seq2)
    record = [[0 for i in range(lseq2 + 1)] for j in range(lseq1 + 1)]
    sub_seq_len = 0
    ind = 0

    for i in range(lseq1):
        for j in range(lseq2):
            if seq1[i] == seq2[j]:
                record[i + 1][j + 1] = record[i][j] + 1
                if record[i + 1][j + 1] > sub_seq_len:
                    sub_seq_len = record[i + 1][j + 1]
                    ind = i + 1
                    if sub_seq_len >= corr_len:
                        return seq1[ind - sub_seq_len:ind]
    return seq1[ind - sub_seq_len:ind]


# Calculate the hamming distance between two sequences
def hammingDist(seq1, seq2):
    return sum(seq11 != seq22 for seq11, seq22 in zip(seq1, seq2))


# Do indel error detection (a.k.a sequence alignment) for sequences with wrong lengths
def clusterWseqIndelDetection(seq_list, cluster_ind, l_seq, shift, window_size, common_len, common_threshold):
    indelseq_list = copy.copy(seq_list)
    rseq_num_indel_total = 0
    wseq_num_indel_total = 0
    for i in tqdm(range(len(cluster_ind))):
        cluster_temp = cluster_ind[i]
        rseq_num = 0
        for j in range(len(cluster_temp)):
            seq_ind = cluster_temp[j]
            seq = seq_list[seq_ind - 1]
            seq = seq[(seq.find(' ') + 1):seq.find('\t')]
            if len(seq) != l_seq:
                break
            rseq_num += 1
        seq_repre = majorityVotingForIndelDetection(seq_list, cluster_temp[:rseq_num], l_seq)
        for j in range(rseq_num):
            seq_ind = cluster_temp[j]
            seq = seq_list[seq_ind - 1]
            index = seq[0:(seq.find(' ') + 1)]
            duplication = seq[seq.find('\t'):]
            seq = seq[(seq.find(' ') + 1):seq.find('\t')]
            if hammingDist(seq, seq_repre) >= 5:
                rseq_num_indel_total += 1
                seq_noindel = slicingHammingDistBasedIndelDetection(seq, seq_repre, l_seq, shift, window_size, common_len, common_threshold)
                if hammingDist(seq_noindel, seq_repre) < 5:
                    rseq_num_indel_total -= 1
                seq = index + seq_noindel + duplication
                indelseq_list[seq_ind - 1] = seq
        for j in range(len(cluster_temp) - 1, -1, -1):

            seq_ind = cluster_temp[j]
            seq = seq_list[seq_ind - 1]
            index = seq[0:(seq.find(' ') + 1)]
            duplication = seq[seq.find('\t'):]
            seq = seq[(seq.find(' ') + 1):seq.find('\t')]
            if len(seq) == l_seq:
                break
            wseq_num_indel_total += 1
            seq_noindel = slicingHammingDistBasedIndelDetection(seq, seq_repre, l_seq, shift,window_size, common_len, common_threshold)
            if hammingDist(seq_noindel, seq_repre) < 5:
                wseq_num_indel_total -= 1
            seq = index + seq_noindel + duplication
            indelseq_list[seq_ind - 1] = seq
    return indelseq_list


# Hamming distance calculation based on slicing window for indel error detection
def slicingHammingDistBasedIndelDetection(seq_indel, seq_repre, l_seq, shift, window_size, common_len, common_threshold):
    while common_len >= common_threshold:
        res = getCommonSubseq(seq_indel, seq_repre, common_len)
        a_ind = seq_indel.find(res)
        b_ind = seq_repre.find(res)
        if seq_indel.count(res) == 1 and seq_repre.count(res) == 1 and abs(a_ind - b_ind) <= shift:
            break
        else:
            common_len -= 1
    if len(res) < common_threshold or common_len >= common_threshold * 3:
        return 'Z' * l_seq
    a = list(seq_indel[a_ind:] + seq_indel[:a_ind])
    b = list(seq_repre[b_ind:] + seq_repre[:b_ind])
    k = common_len - 1
    while k < min(len(a), l_seq):
        if a[k] == b[k]:
            k += 1
            continue
        hammingDist_list_repre = []
        hammingDist_list_indel = []
        for d in range(shift + 1):
            start_ind = k
            end_ind = min(window_size + k, len(a), l_seq)
            if end_ind + d <= l_seq:
                seq_indel_int = str(a[start_ind: end_ind])
                seq_repre_int = str(b[(start_ind + d): (end_ind + d)])
                hammingDist_list_repre.append(hammingDist(seq_indel_int, seq_repre_int))
            if end_ind + d <= len(a):
                seq_indel_int = str(a[(start_ind + d): (end_ind + d)])
                seq_repre_int = str(b[start_ind: end_ind])
                hammingDist_list_indel.append(hammingDist(seq_indel_int, seq_repre_int))
        if hammingDist_list_repre.index(min(hammingDist_list_repre)) == hammingDist_list_indel.index(
                min(hammingDist_list_indel)) == 0:
            k += 1
            continue
        if min(hammingDist_list_repre) < min(hammingDist_list_indel):
            a.insert(k, 'Z')
        elif min(hammingDist_list_repre) > min(hammingDist_list_indel):
            del a[k]
        else:
            k += 1
    if len(a) < l_seq:
        for k in range(l_seq - len(a)):
            a.append('Z')
    if len(a) > l_seq:
        for k in range(len(a) - l_seq):
            del a[-1]
    if len(a) != l_seq or a.count('Z') >= 5:
        a = ['Z'] * l_seq
    if b_ind != 0:
        seq_noindel = a[l_seq - b_ind:] + a[:l_seq - b_ind]
    else:
        seq_noindel = a
    seq_noindel = ''.join(seq_noindel)
    return seq_noindel


# Calculate the cumulative score in modified majority voting
def votingScoreCalculation(seq_str1, seq_str2):
    score = 0
    for i in range(len(seq_str1)):
        if seq_str1[i] == seq_str2[i]:
            score += 1
    return score


# Do modified majority voting to output the representative sequences in each cluster for indel error detection
def majorityVotingForIndelDetection(seq_list, rseq_inds, l_seq):
    voting_counter = [[0] * 4 for q in range(l_seq)]
    Z_ind = []
    Z_base = []

    for j in range(len(rseq_inds)):
        seq_ind = rseq_inds[j]
        seq = seq_list[seq_ind - 1]
        seq_str = seq[(seq.find(' ') + 1): seq.find('\t')]
        repeat_num = int(seq[(seq.find('\t') + 1): (seq.find('\n'))])
        for k in range(l_seq):
            if seq_str[k] == 'A':
                voting_counter[k][0] += repeat_num
            elif seq_str[k] == 'T':
                voting_counter[k][1] += repeat_num
            elif seq_str[k] == 'G':
                voting_counter[k][2] += repeat_num
            elif seq_str[k] == 'C':
                voting_counter[k][3] += repeat_num
    voting_max_result = np.argmax(voting_counter, 1)
    voting_str_list = Num2Base(voting_max_result)

    for k in range(l_seq):
        temp = sorted(voting_counter[k])
        if temp[-1] == temp[-2] == temp[-3] == temp[-4]:
            voting_str_list[k] = 'Z'
            Z_ind.append(k)
            Z_base.append(Num2Base(list(np.argsort(voting_counter[k])[-4:])))
        elif temp[-1] == temp[-2] == temp[-3]:
            voting_str_list[k] = 'Z'
            Z_ind.append(k)
            Z_base.append(Num2Base(list(np.argsort(voting_counter[k])[-3:])))
        elif temp[-1] == temp[-2]:
            voting_str_list[k] = 'Z'
            Z_ind.append(k)
            Z_base.append(Num2Base(list(np.argsort(voting_counter[k])[-2:])))
    voting_str_list_temp = copy.deepcopy(voting_str_list)
    for z in range(len(Z_ind)):
        z_ind = Z_ind[z]
        score = [0] * 4
        for j in range(len(rseq_inds)):
            seq_ind = rseq_inds[j]
            seq = seq_list[seq_ind - 1]
            seq_str = seq[(seq.find(' ') + 1): seq.find('\t')]
            repeat_num = int(seq[(seq.find('\t') + 1): (seq.find('\n'))])

            if seq_str[z_ind] in Z_base[z]:
                base_ind = map_list[seq_str[z_ind]]
                score[base_ind] += repeat_num * votingScoreCalculation(seq_str, ''.join(voting_str_list_temp))
        voting_max_result[z_ind] = np.argmax(score)
    voting_str = ''.join(Num2Base(voting_max_result))
    return voting_str


# Do modified majority voting to correct substitution errors and output the results.
def MajorityVoting(seq_list, clusterind, l_seq, votingfile):
    voting_list = []
    cluster_num = [0] * len(clusterind)
    for i in tqdm(range(len(clusterind))):
        voting_counter = [[0] * 4 for q in range(l_seq)]  # A, T, G, C
        Z_ind = []
        Z_base = []

        for j in range(len(clusterind[i])):
            seq_ind = clusterind[i][j]
            seq = seq_list[seq_ind - 1]
            seq_str = seq[(seq.find(' ') + 1): seq.find('\t')]
            repeat_num = int(seq[(seq.find('\t') + 1): (seq.find('\n'))])
            for k in range(l_seq):
                # try:
                if seq_str[k] == 'A':
                    voting_counter[k][0] += repeat_num
                elif seq_str[k] == 'T':
                    voting_counter[k][1] += repeat_num
                elif seq_str[k] == 'G':
                    voting_counter[k][2] += repeat_num
                elif seq_str[k] == 'C':
                    voting_counter[k][3] += repeat_num
            cluster_num[i] += repeat_num
        voting_max_result = np.argmax(voting_counter, 1)
        voting_str_list = Num2Base(voting_max_result)

        for k in range(l_seq):
            temp = sorted(voting_counter[k])
            if temp[-1] == temp[-2] == temp[-3] == temp[-4]:
                voting_str_list[k] = 'Z'
                Z_ind.append(k)
                Z_base.append(Num2Base(list(np.argsort(voting_counter[k])[-4:])))
            elif temp[-1] == temp[-2] == temp[-3]:
                voting_str_list[k] = 'Z'
                Z_ind.append(k)
                Z_base.append(Num2Base(list(np.argsort(voting_counter[k])[-3:])))
            elif temp[-1] == temp[-2]:
                voting_str_list[k] = 'Z'
                Z_ind.append(k)
                Z_base.append(Num2Base(list(np.argsort(voting_counter[k])[-2:])))
        voting_str_list_temp = copy.deepcopy(voting_str_list)
        for z in range(len(Z_ind)):
            z_ind = Z_ind[z]
            score = [0] * 4
            for j in range(len(clusterind[i])):
                seq_ind = clusterind[i][j]
                seq = seq_list[seq_ind - 1]
                seq_str = seq[(seq.find(' ') + 1): seq.find('\t')]
                repeat_num = int(seq[(seq.find('\t') + 1): (seq.find('\n'))])

                if seq_str[z_ind] in Z_base[z]:
                    base_ind = map_list[seq_str[z_ind]]
                    score[base_ind] += repeat_num * votingScoreCalculation(seq_str, ''.join(voting_str_list_temp))
            voting_max_result[z_ind] = np.argmax(score)
        voting_str = ''.join(Num2Base(voting_max_result))
        voting_list.append(voting_str)
    zip_a_b = zip(cluster_num, voting_list)
    sorted_zip = sorted(zip_a_b, key=lambda x: x[0], reverse=True)
    cluster_num_ordered, voting_list_ordered = zip(*sorted_zip)
    f = open(votingfile, 'w')
    for i in range(len(clusterind)):
        seq = voting_list_ordered[i] + '\n'
        f.write(seq)
    f.close()
    return voting_list_ordered


# Save the cluster file.
def clusterSave(seq_list, clusterind, clusterfile):
    f_r = open(clusterfile, 'w')
    for i in range(len(clusterind)):
        rcluster_index = 'cluster' + str(i + 1) + '\n'
        f_r.write(rcluster_index)
        for j in range(len(clusterind[i])):
            seq = seq_list[clusterind[i][j] - 1]
            f_r.write(seq)
    f_r.close()
    return clusterind


# Compare the output file and the file containing original encoding sequences to evaluate the performance (accuracy,
# coverage, redundancy). RS indicates that the sequences is considered right if the errors are within the error
# correction capability of the RS.

def accuracyCoverageRedundancyMeasurement(votingfile, oriseqfile, RS=0):
    f_ori = open(oriseqfile, 'r')
    ori_seq_list = f_ori.readlines()
    f_ori.close()
    f_v = open(votingfile, 'r')
    vseq_list = f_v.readlines()
    f_v.close()
    l = 160
    acc_list = [0] * len(ori_seq_list)
    for i in tqdm(range(len(vseq_list))):
        seq = vseq_list[i].strip('\n')
        if len(seq) == l:
            for j in range(len(ori_seq_list)):
                seq_ori = ori_seq_list[j].strip('\n')
                byte_ind = set()
                for k in range(len(seq)):
                    if seq[k] != seq_ori[k]:
                        byte_ind.add(k // 4)
                    if len(byte_ind) > RS // 2:
                        break

                if len(byte_ind) <= RS // 2:
                    acc_list[j] += 1
                    break
        else:
            continue
    accuracy = sum(acc_list) / len(vseq_list)
    coverage = np.count_nonzero(acc_list) / len(ori_seq_list)
    redundancy = max(0.0, len(vseq_list) / len(ori_seq_list) - 1)
    print('Accuracy = ' + str(accuracy))
    print('Coverage = ' + str(coverage))
    print('Redundancy = ' + str(redundancy))
    acc_np = np.array(acc_list)
    acc_index = np.where(acc_np == 0)[0]
    acc_ind = [m + 1 for m in acc_index]
    print('These sequences are not covered.' + str(acc_ind))
    return accuracy, coverage, redundancy


# The following functions are designed to save and load objects, which can be called on demand.
def objectSave(obj, filename):
    f = open(filename, 'wb')
    pickle.dump(obj, f)
    f.close()


def objectLoad(filename):
    f = open(filename, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj

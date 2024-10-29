import time
from math import floor
import numpy as np
from tqdm import tqdm


def set_rank(a_dict):
    a_sort_list = sorted(a_dict.items(), key=lambda x: x[1], reverse=True)
    a_sort_dict = dict(a_sort_list)
    # for n, s in a_sort_list:
    #     a_sort_dict[n] = s
    return a_sort_dict

def accuracyCoverageRedundancyMeasurement(votingfile, oriseqfile,RS):
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
                        byte_ind.add(k//4)
                    if len(byte_ind) > RS//2:
                        break

                if len(byte_ind) <= RS//2:
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


# a_set = [0.1, 0.5, 1, 2, 3, 4]
# a_set = [0.1, 0.2, 0.5, 1, 1.5]
# a_set = [0.1, 0.3, 0.5, 0.7, 1, 1.5]
# for cc in a_set:
# file = 'R2_q30.txt'
# file_repeat = 'DNA_mix_order_RS6.txt'
#     cc=5
#     file = 'DNA_mix'+str(cc)+'.txt'
#     file_repeat = 'DNA_mix'+str(cc)+'_order.txt'
# f = open(file, 'r')
# f2 = open(file_repeat, 'w')
# seq_list = f.readlines()
# l = 160
# seq_revised = dict()
# for i in tqdm(range(0, len(seq_list))):
#     seq = seq_list[i].strip('\n')
#     if seq_revised.get(seq) is not None:
#         seq_revised[seq] += 1
#     else:
#         seq_revised[seq] = 1
#
# seq_nums = 0
# seq_revised = set_rank(seq_revised)
# for seq, repeat in seq_revised.items():
#     if len(seq) == l:
#         f2.write(seq + '\n')
#
# accuracyCoverageRedundancyMeasurement(file_repeat,'medenc6.txt',6)
#
# f.close()
# f2.close()

file_repeat = 'P1.txt'
file_order = 'P1_order.txt'
f_r = open(file_repeat,'r')
f_o = open(file_order, 'w')
seq_list = f_r.readlines()
for seq in seq_list:
    seq_str = seq[:seq.find('\t')]
    if len(seq_str) == 160:
        f_o.write(seq_str + '\n')
f_r.close()
f_o.close()
accuracyCoverageRedundancyMeasurement(file_order,'rs2.txt',0)


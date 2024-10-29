import random
import scipy
from tqdm import tqdm
from DNAHandler import DNAHandler


def sub_file_append(seq_seq_list, seq_errfile, sub_prob):
    d = DNAHandler()
    seq_err_list = d.add_rand_sub_new(seq_seq_list, sub_prob, True)
    f_err = open(seq_errfile, 'w')
    for seq_err in seq_err_list:
        f_err.write(seq_err + '\n')


def ins_file_append(seq_seq_list, seq_errfile, ins_prob):
    d = DNAHandler()
    seq_err_list = d.add_rand_ins_new(seq_seq_list, ins_prob, True)
    f_err = open(seq_errfile, 'w')
    for seq_err in seq_err_list:
        f_err.write(seq_err + '\n')


def del_file_append(seq_seq_list, seq_errfile, del_prob):
    d = DNAHandler()
    seq_err_list = d.add_rand_del_new(seq_seq_list, del_prob, True)
    f_err = open(seq_errfile, 'w')
    for seq_err in seq_err_list:
        f_err.write(seq_err + '\n')


def mix_file_append(seq_seq_list, seq_errfile, sub_prob, ins_prob, del_prob):
    d = DNAHandler()
    l = len(seq_seq_list)

    seq_seq_list1 = seq_seq_list[:l//2]
    seq_err_list1 = d.add_rand_ins_new(seq_seq_list1, ins_prob * 2, True)
    seq_seq_list2 = seq_seq_list[l//2:]
    seq_err_list2 = d.add_rand_del_new(seq_seq_list2, del_prob * 2, True)
    seq_err_list = seq_err_list1 + seq_err_list2
    seq_err_list = d.add_rand_sub_new(seq_err_list, sub_prob, True)
    random.shuffle(seq_err_list)
    # seq_err_list = d.add_rand_sub_new(seq_seq_list, sub_prob, True)
    # seq_err_list = d.add_rand_indel_new(seq_err_list, ins_prob, del_prob, True)
    f_err = open(seq_errfile, 'w')
    for seq_err in seq_err_list:
        f_err.write(seq_err + '\n')


if __name__ == '__main__':
    seq_file = 'medenc6.txt'
    seq_sub_file = 'DNA_sub5.txt'
    seq_ins_file = 'DNA_ins5.txt'
    seq_del_file = 'DNA_del.txt'

    seq_seq_file = 'DNA_seq_RS6.txt'

    f_seq = open(seq_seq_file, 'w')
    l = 160
    f = open(seq_file, 'r')
    seq_list = f.readlines()
    seq_seq_list = []

    poisson_mean = 1000
    for seq in tqdm(seq_list):
        seq_num = scipy.stats.poisson.rvs(poisson_mean)
        seq = seq.strip('\n')
        for i in range(seq_num):
            seq_seq_list.append(seq)
    random.shuffle(seq_seq_list)
    for seq in seq_seq_list:
        f_seq.write(seq + '\n')

    # seq_seq_list = f_seq.readlines()
    # for i in range(len(seq_seq_list)):
    #     seq_seq_list[i] = seq_seq_list[i].strip('\n')

    # exit()
    # sub_prob = .05
    # ins_prob = .0
    # del_prob = .0
    # a_set = [0.1, 0.5, 1, 2, 3, 4]
    # a_set = [0.1, 0.3, 0.5, 0.7, 1, 1.5]
    # for a in a_set:
        # a = 5
    # seq_mix_file = 'DNA_mix' + str(a) + '.txt'
    # sub_prob = .02
    # ins_prob = .005
    # del_prob = .01 * a
    seq_mix_file = 'DNA_mix_RS6.txt'
    sub_prob = .05
    ins_prob = .005
    del_prob = .005
    # sub_file_append(seq_seq_list, seq_sub_file, sub_prob)
    # ins_file_append(seq_seq_list, seq_ins_file, ins_prob)
    # del_file_append(seq_seq_list, seq_del_file, del_prob)
    mix_file_append(seq_seq_list, seq_mix_file, sub_prob, ins_prob, del_prob)




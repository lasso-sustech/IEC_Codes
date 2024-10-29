import time
from math import floor

from tqdm import tqdm


def set_rank(a_dict):
    a_sort_list = sorted(a_dict.items(), key=lambda x: x[1], reverse=True)
    a_sort_dict = dict(a_sort_list)
    # for n, s in a_sort_list:
    #     a_sort_dict[n] = s
    return a_sort_dict

# cc = 5
# a_set = [0.1, 0.5, 1, 2, 3, 4]
# a_set = [0.1, 0.2, 0.5, 1, 1.5]
a_set = [0.1, 0.3, 0.5, 0.7, 1, 1.5]

# for cc in a_set:
#     file = 'DNA_mix'+str(cc)+'.txt'
#     file_repeat = 'DNA_mix'+str(cc)+'_repeat.txt'
file = 'q30.txt'
file_repeat = 'q30_repeat.txt'

f = open(file, 'r')
f2 = open(file_repeat, 'w')
seq_list = f.readlines()

seq_revised = dict()
for i in tqdm(range(0, len(seq_list))):
    seq = seq_list[i].strip('\n')
    if seq_revised.get(seq) is not None:
        seq_revised[seq] += 1
    else:
        seq_revised[seq] = 1

seq_nums = 0
seq_revised = set_rank(seq_revised)
for seq, repeat in seq_revised.items():
    seq_str = seq + '\t' + str(repeat) + '\n'
    seq_nums = seq_nums + repeat
    f2.write(seq_str)

f.close()
f2.close()



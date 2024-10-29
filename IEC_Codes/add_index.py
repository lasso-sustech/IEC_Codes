from tqdm import tqdm

# f = open('DNA_repeat.txt','r')
# f_ind = open('DNA_index.txt','w')
# a_set = [0.1, 0.5, 1, 2, 3, 4]
# a_set = [0.1, 0.2, 0.5, 1, 1.5]
a_set = [0.1, 0.3, 0.5, 0.7, 1, 1.5]
# for cc in a_set:
# cc=5
seq_num = 0
    # f = open('DNA_mix'+str(cc)+'_repeat.txt','r')
    # f_ind = open('DNA_mix'+str(cc)+'_index.txt','w')
f = open('P1.txt','r')
f_ind = open('P1_index.txt','w')
seq = f.readlines()
index = range(1, len(seq) + 1)
index_seq_list = []
for i in tqdm(index):
    index_seq = str(i) + ' ' + seq[i - 1]
    index_seq_list.append(index_seq)
f_ind.writelines(index_seq_list)
f.close()
f_ind.close()
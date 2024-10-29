from utils import *
import argparse

parser = argparse.ArgumentParser(description='Test for argparse', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--length', '-l', help='length of the encoding sequences', required=True, type=int)
parser.add_argument('--deltal', '-d', help='maximum length difference from the correct length', default=1)
parser.add_argument('--alpha', '-a', help='value of alpha', default=0.3, type=float)
parser.add_argument('--lint', '-li', help='value of interval length', default=20, type=int)

parser.add_argument('--shift', '-s', help='maximum shift of the window', default=1, type=int)
parser.add_argument('--commonlen', '-cl', help='length of the common sequence', default=5, type=int)
parser.add_argument('--commonth', '-ct', help='threshold of the length of the common sequence', default=3, type=int)
parser.add_argument('--windowsize', '-ws', help='window size', default=5, type=int)


parser.add_argument('--input', '-f', help='name of the input file', required=True)
parser.add_argument('--rlfile', '-r', help='name of the file containing sequences with the right length',
                    default='seq_rightlen.txt', type=str)
parser.add_argument('--wlfile', '-w', help='sequences with the wrong length', default='seq_wronglen.txt', type=str)
parser.add_argument('--cluster', '-c', help='clusters', default='cluster.txt', type=str)
parser.add_argument('--idcluster', '-idc', help='clusters after indel error detection and sequence alignment',
                    default='idcluster.txt', type=str)
parser.add_argument('--original', '-e', help='original encoding sequences', default=None, type=str)
parser.add_argument('--output', '-o', help='name of the output file', default='output.txt', type=str)

args = parser.parse_args()

if __name__ == '__main__':
    DNAseqfile = args.input
    rightlenfile = args.rlfile
    wronglenfile = args.wlfile
    clusterfile = args.cluster
    idclusterfile = args.idcluster
    outputngfile = args.output
    oriseqfile = args.original

    l_seq = args.length
    deltal = args.deltal

    print('IEC begins!')
    print('Preprocessing: Input file loading and Length classification------------------')

    seq_list, rseqind_list, wseqind_list = seqLenClassify(DNAseqfile, rightlenfile, wronglenfile, l_seq,
                                                          deltal)

    print('The first phase: Clustering------------------')

    alpha = args.alpha
    l_int = args.lint
    maxdis = alpha * l_int


    clusterind = LvstClustering(seq_list, rseqind_list, wseqind_list, l_int, maxdis)
    clusterSave(seq_list, clusterind, clusterfile)
    print('The second phase: Indel error detection------------------')

    shift = args.shift
    window_size = args.windowsize
    common_len = args.commonlen
    common_threshold = args.commonth
    indelseq_list = clusterWseqIndelDetection(seq_list, clusterind, l_seq, shift, window_size, common_len, common_threshold)
    clusterSave(indelseq_list, clusterind, idclusterfile)
    print('The third phase: Majority voting------------------')
    voting_list = MajorityVoting(indelseq_list, clusterind, l_seq, outputngfile)
    if oriseqfile:
        print('Evaluating the performance------------------')
        accuracy, coverage, redundancy = accuracyCoverageRedundancyMeasurement(outputngfile, oriseqfile, 0)


    print('IEC finishes! Please refer to the output file.')

import os

import numpy as np
import pandas as pd
from multiprocessing import Pool
from generate_pssm import *
# from Bio import SeqIO

aa_alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
root_path = ''

class Seq():
    def __init__(self, id, seq):
        self.name = id[1:]
        self.seq = seq
        self.id =id[1:]

def parser_fasta(fasta):
    ids, seqs = load_seqs(fasta)
    res = []
    for i in range(len(ids)):
        res.append(Seq(ids[i], seqs[i]))
    return res

def read_fasta(fn):
    """
    read fasta file
    :param fn: file name
    :return: AA sequences, AA names
    """
    seq_names = []
    seqs = []
    with open(fn, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0] == '>':
                # 去掉'|'后的标签，因为可能是程序自己加的
                if '|' in line:
                    line = line.split('|')[0]
                seq_names.append(line)

            else:
                seqs.append(line)
    return seqs, seq_names

def gen_features(fn, pt, feature_name, pssm_file = None):
    """

    :param fn: file name
    :param pt: peptide type
    :param feature_name: target feature name
    :param pssm_file: pssm file path, some features need generate pssm matrix first
    :return: generate feature files under user folder
    """
    print('--------------generating feature', feature_name,'--------------')
    if feature_name == 'AAC':
        return AAC(fn, pt)
    # elif feature_name == 'BIT20':
    #     return BIT20NT4(fn, pt)
    elif feature_name == 'Kmer':
        return Kmer(fn, pt)
    elif feature_name == 'Top_n_gram':
        return Top_n_gram(fn, pt)
    elif feature_name == 'DP':
        return DP(fn, pt)
    elif feature_name == 'DT':
        return DT(fn, pt)
    elif feature_name == 'DR':
        return DR(fn, pt)
    elif feature_name == 'PseAAC':
        return PseAAC(fn, pt)
    else:
        feature = Feature()
        feature.get_dbt_feature(fn , feature_name, pt, pssm_file)


def AAC(fn, pt):
    """
    AAC: Amino acid component
    :param fn: user folder name
    :param pt: peptide name
    :return: list of AAC
    """
    seqs, _ = read_fasta(fn + '/test.fasta')
    res = []
    aa_dict = {}
    for i in range(len(aa_alphabet)):
        aa_dict[aa_alphabet[i]] = 0
    for seq in seqs:
        tmp = aa_dict.copy()
        n = len(seq)
        for a in seq:
            tmp[a] += 1
        for a in aa_alphabet:
            tmp[a] /= n
        res.append(list(tmp.values()))

    df = pd.DataFrame(res)
    df.to_csv(fn + '/' + pt + '_AAC.txt', header=None, index=None, sep=',')

    return res


def Top_n_gram(fn, pt):
    """
    creating Top-n-gram matirx under user dir
    """
    cmd = 'python ' + root_path + 'tool/BioSeq-Analysis2/profile.py ' + \
          fn + '/test.fasta -method Top-n-gram -f csv -out ' + \
          fn + '/' + pt + '_Top-n-gram.txt -n 2'
    os.system(cmd)


def DP(fn, pt):
    """
    creating DP matirx under user dir
    """
    cmd = 'python ' + root_path + 'tool/BioSeq-Analysis2/nac.py ' + \
        fn + '/test.fasta Protein DP -f csv -out ' + \
        fn + '/' + pt + '_DP.txt'
    os.system(cmd)

def DT(fn, pt):
    """
    creating DT matirx under user dir
    """
    cmd = 'python ' + root_path + 'tool/BioSeq-Analysis2/profile.py ' + \
          fn + '/test.fasta -method DT -f csv -out ' + \
          fn + '/' + pt + '_DT.txt -max_dis 3'
    os.system(cmd)

def DR(fn, pt):
    """
    creating DR matirx under user dir
    """
    cmd = 'python ' + root_path + 'tool/BioSeq-Analysis2/nac.py ' + \
        fn + '/test.fasta Protein DR -f csv -out ' + \
        fn + '/' + pt + '_DR.txt'
    os.system(cmd)

def Kmer(fn, pt):
    """
    creating Kmer matirx under user dir
    """
    cmd = 'python ' + root_path + 'tool/BioSeq-Analysis2/nac.py ' + \
        fn + '/test.fasta Protein Kmer -f csv -out ' + \
        fn + '/' + pt + '_Kmer.txt -k 2'
    os.system(cmd)

def PseAAC(fn, pt):
    """
    creating PseAAC matirx under user dir
    """
    cmd = 'python ' + root_path + 'tool/BioSeq-Analysis2/pse.py ' + fn + '/test.fasta Protein -method PC-PseAAC-General ' + \
        '-f csv -out ' + fn + '/' + pt + '_PseAAC.txt -lamada 2 -w 0.3 -labels -1'

    os.system(cmd)

def gen_peps(fn, pt):
    """
    creating Bit20NTCT, CTD, GAP5, Ngram(N=1)1, OVNT(5)84
    """
    print('--------------generating feature', 'BIT CTD GAP5 Ngram OVNT(5)84', '--------------')
    cmd = 'python ' + root_path + 'gen_pepred_features.py -src ' + fn + ' -pt ' + pt
    os.system(cmd)


class Feature():
    """
    Generate features including PPCT, PSSM_DT and PSFM_DBT and save them under user idr
    """
    co = 0
    def __init__(self):
        co = 0
        self.seq_len = 700
        self.coding = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L',
        'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

        self.blosum62 = {}
        # TODO
        # blosum_reader = open('./tool/psiblast/blosum62', 'r')
        blosum_reader = open(root_path + 'tool/psiblast/blosum62', 'r')
        count = 0
        for line in blosum_reader:
            count = count + 1
            if count <= 7:
                continue
            line = line.strip('\r').split()
            self.blosum62[line[0]] = [float(x) for x in line[1:21]]

    def get_protein_blosum(self, protein):
        seq_str = str(protein.seq).replace('X', '')
        seq_str = seq_str.replace('B', '')
        seq_str = seq_str.replace('Z', '')
        seq_str = seq_str.replace('U', '')
        seq_str = seq_str.replace('J', '')
        seq_str = seq_str.replace('O', '')
        protein_lst = []
        for aa in seq_str:
            aa = aa.upper()
            protein_lst.append(self.blosum62[aa])
        return np.array(protein_lst)

    def generate_protein_psfm(self, protein):
        seq_str = str(protein.seq).replace('X', '')
        seq_str = seq_str.replace('B', '')
        seq_str = seq_str.replace('Z', '')
        seq_str = seq_str.replace('U', '')
        seq_str = seq_str.replace('J', '')
        seq_str = seq_str.replace('O', '')
        mat=np.zeros([len(seq_str),20])
        i=0
        for aa in seq_str:
            aa = aa.upper()
            indx = self.coding.index(aa)
            mat[i][indx]=1.0
            i += 1
        return mat

    def one_hot(self, file_path):
        # seq_record = list(SeqIO.parse(file_path, 'fasta'))
        seq_record = parser_fasta(file_path)
        N=len(seq_record)
        mats = np.asarray(
            np.zeros([N, self.seq_len, 20]))
        i=0
        for prot in seq_record:
            if len(prot.seq)<self.seq_len:
                seqLen=len(prot.seq)
            else:
                seqLen=self.seq_len
            for j in range(seqLen):
                e=list(prot.seq)[j]
                indx=self.coding.index(e)
                mats[i][j][indx]=1.0
            i+=1
        return mats

    def read_pssm(self, pssm_file):
        with open(pssm_file, 'r') as f:
            lines = f.readlines()
            lines = lines[3:-6]
            pro_seq=[]
            mat = []
            for line in lines:
                tmp = line.strip('\n').split()
                if len(tmp)==0:
                    break
                pro_seq.append(tmp[1])
                tmp = tmp[2:22]
                mat.append(tmp)
            mat = np.array(mat)
            mat = mat.astype(float)
            return pro_seq, mat

    def read_psfm(self, pssm_file):
        with open(pssm_file, 'r') as f:
            lines = f.readlines()
            lines = lines[3:-6]
            pro_seq=[]
            mat = []
            for line in lines:
                tmp = line.strip('\n').split()
                if len(tmp)==0:
                    break
                pro_seq.append(tmp[1])
                tmp = tmp[22:42]
                mat.append(tmp)
            mat = np.array(mat)
            mat = mat.astype(float)
            mat = np.divide(mat, 100)
            return pro_seq,mat

    def average(self, matrixSum, seqLen):
        matrix_array = np.array(matrixSum)
        matrix_array = np.divide(matrix_array, seqLen)
        matrix_array_shp = np.shape(matrix_array)
        matrix_average = [(np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1],)))]
        return matrix_average

    def sigmoid(self, x):
        s = 1 / (1 + np.exp(-x))
        return s

    def preHandleColumns(self, PSFM, PSSM, STEP, feature):
        PSSM=np.asarray(PSSM, float)
        PSFM=np.asarray(PSFM, float)
        mat = np.zeros((20,20),float)

        seq_cn = np.shape(PSSM)[0]
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    if feature == 'PSSM_DT' : mat[i][j] += (PSSM[k][i] * PSSM[k + STEP][j])
                    elif feature == 'PSFM_DBT': mat[i][j] += (PSFM[k][i] * PSFM[k + STEP][j])
                    elif feature == 'PPCT': mat[i][j] += (PSSM[k][i] * PSFM[k + STEP][j] + PSFM[k][i]* PSSM[k+STEP][j] +PSFM[k][i] * PSFM[k + STEP][j] +PSSM[k][i] * PSSM[k + STEP][j])

        return mat

    def psfm_dbt(self, PSFM, PSSM, end, feature):
        seq_cn = float(np.shape(PSSM)[0])
        vector = []

        for i in range(0, end + 1):
            matrix = self.preHandleColumns(PSFM, PSSM, i, feature)
            ksb_vector = self.average( matrix, float(seq_cn - i))
            vector += list(ksb_vector[0])

        return vector

    #PSSM-DT PSFM-DBT PPCT
    def get_dbt_feature(self, source_path, feature, dataset, pssm_dir=None):
        file_path = source_path + '/test.fasta'
        print('Reading:', file_path)
        save_path = source_path + '/' + dataset + '_' + feature + '.txt'
        if feature == 'PSSM_DT': end = 5
        elif feature == 'PSFM_DBT': end = 4
        elif feature == 'PPCT': end = 4
        print('save_path:', save_path)

        read_count = 0
        generate_count = 0
        if os.path.isfile(save_path):
            mats = np.loadtxt(save_path, dtype = np.float)
            return mats
        else:
            # seq_record = list(SeqIO.parse(file_path, 'fasta'))
            seq_record = parser_fasta(file_path)
            mats = []
            n = len(seq_record)
            count = 0
            for prot in seq_record:
                n -= 1
                if '|' in prot.name:
                    pssm_path = pssm_dir + '/' + (prot.name).split('|')[0]+'_'+(prot.name).split('|')[1] + '.pssm'
                else:
                    pssm_path = pssm_dir + '/' + prot.name + '.pssm'
                if os.path.isfile(pssm_path):
                    read_count += 1
                    pro_seq, pssm_profile=self.read_pssm(pssm_path)
                    pro_seq, psfm_profile = self.read_psfm(pssm_path)
                    print('read pssm '+ str(prot.name)  + ': ' +str(n))
                else:
                    generate_count += 1
                    print('generate pssm '+prot.name  + ': ' +str(n))
                    pssm_profile = self.get_protein_blosum(prot)
                    psfm_profile = self.generate_protein_psfm(prot)
                vec = self.psfm_dbt(psfm_profile, pssm_profile,end,feature)
                mats.append(vec)

            print('mats', np.array(mats).shape)
            np.savetxt(save_path, mats, fmt ='%.12f')

            print('read',read_count,'generate',generate_count)

if __name__ == '__main__':
    fn = 'test/'
    pt = 'AAP'
    gen_peps(fn,  pt)

    #
    # #test Kmer
    # gen_features(fn, pt, 'Kmer')
    #
    # #test DR
    # gen_features(fn, pt, 'DR')
    #
    # test PseAAC
    # gen_features(fn, pt, 'PseAAC')
    #
    # # test Top-n-gram
    # gen_features(fn, pt, 'Top_n_gram')
    #
    # # test DT
    # gen_features(fn, pt, 'DT')
    #
    # # test DP
    # gen_features(fn, pt, 'DP')

    # gen_pssm(fn, fn + '/profile')
    #
    # pssm_path = '../tmp/test/profile/pssm'
    # # test PPCT
    # gen_features(fn, pt, 'PPCT', pssm_path)
    # # test PSSM_DT
    # gen_features(fn, pt, 'PSSM_DT', pssm_path)
    # # test PSFM_DBT
    # gen_features(fn, pt, 'PSFM_DBT', pssm_path)


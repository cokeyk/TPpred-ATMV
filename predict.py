import shutil

import numpy as np
import zipfile
import argparse
import os
from multiprocessing import Pool
from generate_features import gen_features, read_fasta, gen_peps
from generate_pssm import gen_pssm
import matlab
import matlab.engine
import warnings
from scipy.io import loadmat
from utils import check_rewrite_data

warnings.filterwarnings("ignore")

model_path = 'matlab_models'

def get_features(source_folder, pt):
    """
    Get 10 features
    :param args:
    :return:
    """
    # pool = Pool(3)
    fea1 = ['DP','Top_n_gram','Kmer','DT','DR','PseAAC']
    for fea in fea1:
        # pool.apply_async(gen_features, (source_folder, pt, fea,))
        gen_features(source_folder, pt, fea)

    gen_pssm(source_folder, source_folder + '/profile')
    fea2 = ['PPCT','PSFM_DBT','PSSM_DT']
    pssm_file = source_folder + '/profile/pssm'
    for fea in fea2:
        # pool.apply_async(gen_features, (source_folder, pt, fea,pssm_file))
        gen_features(source_folder, pt, fea,pssm_file)
    # pool.close()
    # pool.join()

    gen_peps(source_folder ,pt)

def predict_py(ttdata, W, alpha):
    """
    :param ttdata: (n_view, n, m)
    :param W:   (n_view, out, m)
    :param alpha: n_view
    :return:
    """
    n_samples = ttdata[0].shape[0]

    pre_label = np.zeros((2, n_samples))
    for i in range(len(ttdata)):

        pre_label = pre_label + alpha[i] * np.matmul(W[i], ttdata[i].transpose())

    pos_score = []
    # softmax
    for i in range(n_samples):
        tmp = np.max(pre_label[:, i])
        pos_score.append(np.exp(pre_label[0, i] - tmp) / (np.exp(pre_label[0, i] - tmp) + np.exp(pre_label[1, i] - tmp)))
    return np.array(pos_score)

def predict(args):

    source_folder = args.src  # current user folder
    fasta_file = args.fasta   # fasta name
    pt = args.type  # peptide type
    th = args.th # threshold
    result_folder = args.o

    if source_folder[-1] != '/': source_folder += '/'
    if result_folder[-1] != '/': result_folder += '/'


    with open(source_folder + fasta_file, 'r') as f:
        file_seqs = f.readlines()
    check_rewrite_data(file_seqs, source_folder +'test.fasta' )

    # if not os.path.exists(source_folder + 'test.fasta'):
    #     shutil.copyfile(source_folder + fasta_file, source_folder + 'test.fasta')
    # get sequence names

    seqs, seq_names = read_fasta(source_folder + fasta_file)

    # Generate features and save them in the path source_folder
    get_features(source_folder, pt)

    # load data by matlab
    eng = matlab.engine.start_matlab()
    X = eng.read_data(source_folder + '/', pt)
    X = [np.array(x) for x in X]

    # get matlab model parameters
    Mat = loadmat(model_path + '/' + pt + '/' + pt + '.mat')
    W = Mat['W'].tolist()[0]
    alpha = np.squeeze(Mat['alpha'], axis=1).tolist()

    # predict score
    pos_score = predict_py(X, W, alpha)

    # predict class
    pred_clss = np.zeros(len(pos_score))
    pred_clss[pos_score >= th] = 1

    print('Writing results!')


    # save result files under result folder

    if not os.path.exists(result_folder):
        os.makedirs(result_folder)

    with open(result_folder + 'name.txt', 'w') as f:
        for i in seq_names:
            f.write(str(i))
            f.write('\n')

    with open(result_folder + 'prob.txt', 'w') as f:
        for i in pos_score:
            f.write(str(float('%.3f' % i)))
            f.write('\n')
    with open(result_folder + 'class.txt', 'w') as f:
        for i in pred_clss:
            f.write(str(int(i)))
            f.write('\n')

    # compress result files

    print('Writing over ! Results can be found in folder', result_folder)


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-src', type=str, help='The input folder which contains fasta format target file.'
                       , default= 'test')
    parse.add_argument('-fasta', type=str, help='Target file name in fasta format.'
                       , default='AAPT.txt')

    parse.add_argument('-type', type=str, help='The peptide type to be predicted.', default='AAP')
    parse.add_argument('-th', type=float, help='The threshold for predicting.', default=0.5)
    parse.add_argument('-o', type=str, help='Output file', default='result')
    args = parse.parse_args()
    print('args:', args)
    predict(args)

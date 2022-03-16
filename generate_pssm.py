# -*- coding: utf-8 -*-
import sys, os, subprocess
import argparse
from multiprocessing import Pool

# TODO
BLAST = os.getcwd() + '/tool/psiblast/psiblast'
BLAST_DB = os.getcwd() +'/tool/psiblast/nrdb90/nrdb90'

complet_n=0

def load_seqs(fn):
    # 加载序列数据
    ids = []
    seqs = []
    id = 0

    with open(fn, 'r') as f:
        lines = f.readlines()
        seq_tmp = ""

        for i, line in enumerate(lines):
            line = line.strip()
            if line[0] == '>':
                id = line.replace('|','_')
                id = id.split(' ')[0]
            elif i < len(lines) -1 and lines[i+1][0] != '>':
                seq_tmp += line.strip()
            else:
                seq_tmp += line.strip()
                seqs.append(seq_tmp)
                ids.append(id)
                id = 0
                seq_tmp = ""

    return ids, seqs

def generateFasta(input_file, profile_home):
    print('Generating fasta files...')
    fasta_dir = profile_home+'/fasta'
    if not os.path.isdir(fasta_dir):
        os.makedirs(fasta_dir)
    ids, seqs = load_seqs(input_file)
    with open(profile_home + '/id_list.txt', 'w') as f:
        for id in ids:
            f.write(id[1:] + '\n')

    for i, id in enumerate(ids):
        name = id[1:]
        fasta_file = fasta_dir + '/' + name + '.fasta'

        with open(fasta_file, 'w') as wf:
            wf.write('>' + name + '\n')
            wf.write(seqs[i] + '\n')


def run_search(fdz):
    fd = fdz[0]
    profile_home = fdz[1]
    protein_name = fd.split('.')[0]
    global complet_n
    complet_n += 1
    print('Processing:%s---%d' % (protein_name, complet_n*6))
    outfmt_type = 5
    num_iter = 3
    # evalue_threshold = 0.05
    evalue_threshold = 0.001

    fasta_file = profile_home + '/fasta/' + protein_name + '.fasta'
    xml_file = profile_home + '/xml/' + protein_name + '.xml'
    pssm_file = profile_home + '/pssm/' + protein_name + '.pssm'
    if os.path.isfile(pssm_file):
        pass
    else:
        cmd = ' '.join([BLAST,
                        '-query ' + fasta_file,
                        '-db ' + BLAST_DB,
                        '-out ' + xml_file,
                        '-evalue ' + str(evalue_threshold),
                        '-num_iterations ' + str(num_iter),
                        '-outfmt ' + str(outfmt_type),
                        '-out_ascii_pssm ' + pssm_file,  # Write the pssm file
                        '-num_threads ' + '2']
                       )
        return_code = subprocess.call(cmd, shell=True)

def run_blast(profile_home):
    print('Generating PSSM:')
    fasta_dir = profile_home + '/fasta'
    seq_DIR = os.listdir(fasta_dir)
    pssm_dir = profile_home + '/pssm'
    if not os.path.isdir(pssm_dir):
        os.makedirs(pssm_dir)
    xml_dir = profile_home + '/xml'
    if not os.path.isdir(xml_dir):
        os.makedirs(xml_dir)

    profile_homes = [profile_home for i in seq_DIR]
    pool = Pool(4)
    results = pool.map(run_search, zip(seq_DIR, profile_homes))
    pool.close()
    pool.join()

def gen_pssm(input_path, save_path):
    profile_home = save_path
    fasta_file = input_path + '/test.fasta'
    if not os.path.isdir(profile_home):
        os.makedirs(profile_home)
    generateFasta(fasta_file, profile_home)
    run_blast(profile_home)

def main(args):
    file_path = args.fasta_file
    profile_home = os.path.split(file_path)[0] + '/' + os.path.split(file_path)[1].split('.')[0]
    print(profile_home)
    if not os.path.isdir(profile_home):
        os.makedirs(profile_home)
    generateFasta(file_path, profile_home)
    run_blast(profile_home)


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--fasta_file', type=str, default='../tmp/test/test.fasta')
    # args = parser.parse_args()
    #
    # print('Start!', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    # main(args)
    # print('End!', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    gen_pssm('../tmp/test','../tmp/test/profile')

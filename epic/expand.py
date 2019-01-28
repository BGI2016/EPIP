import pandas as pd
import os
import subprocess
import pkg_resources
import shutil
from sklearn.externals import joblib
import pickle


package = "epic"
tmpdir = "{}/{}".format(os.getcwd(), "tmp")

def formatPSSM(df, allele):
    df['HUMAN1'] = 'HUMAN'
    df['HUMAN2'] = 'HUMAN'
    df['class'] = 'class-1'
    df['HLA'] = allele
    df['YES'] = 'YES'
    df['length'] = df.Peptide.str.len()
    return df[['HUMAN2','HUMAN2','class','HLA','YES','Peptide','length']]

def train_PSSM(para):
    global tmpdir
    os.makedirs(tmpdir, exist_ok=True)

    allele = para.allele
    peptide = pd.read_csv(para.input, header=None, names=['Peptide'])
    for l in peptide.Peptide.str.len().drop_duplicates():
        peptide_ls = peptide[peptide.Peptide.str.len() == l]
        if peptide_ls.shape[0] <= 100:
            print("The data size of training {} {}mer PSSM is less than 100, which might lead to less "
                  "accurate predictions.".format(allele, l))

    pssm = formatPSSM(peptide, allele)
    pssm_train = "{}/pssm_train.txt".format(tmpdir)
    pssm.to_csv(pssm_train, sep="\t", header=None, index=False)

    pssm_dir = pkg_resources.resource_filename(package, 'model/PSSM/PSSMx/')
    train_script = pkg_resources.resource_filename(package, 'MHCItrain.pl')
    train_command = "perl %s %s HUMAN NA %s" % (train_script, pssm_train, pssm_dir)
    subprocess.check_output(train_command, shell=True)

    shutil.rmtree(tmpdir)

def learn_length_dist(para):
    length_dist_pkl = 'model/allele_length_distribution.pkl'
    length_dist_file = pkg_resources.resource_filename(package, length_dist_pkl)
    length_dist = joblib.load(length_dist_file)

    peptide = pd.read_csv(para.input, header=None, names=['Peptide'])
    peptide = peptide.drop_duplicates('Peptide')

    allele = para.allele.split(r"-")[1]
    length_dist[allele] = {l:5 for l in range(8,16)}

    for l in range(8,16):
        length_dist[allele][l] = length_dist[allele].get(l) + \
                                          peptide[peptide.Peptide.str.len() == l].shape[0]

        length_dist[allele] = {l: ((length_dist[allele][l]) /
                                            sum(length_dist[allele].values())) for l in range(8, 16)}

    pickle.dump(length_dist, open(length_dist_file, "wb"))

def retrain(para):
    train_PSSM(para)
    learn_length_dist(para)

#Copyright (c) 2018, BGI-Shenzhen

import os
import numpy as np
import pandas as pd
from sklearn.externals import joblib
import pkg_resources
from collections import defaultdict
import uuid
import shutil


package = "epip"
tmpdir = "{}/{}".format(os.getcwd(), "tmp")


def valid_peptide(df):
    if df.isnull().values.any():
        raise Exception("Error: peptide file contains NaN values")

    if df.shape[1] > 1:
        raise Exception("Error: more than one column in peptide file.")

    invalid_peptide_index = df.iloc[:, 0].str.contains(r'[a-z]|U|B|J|O|Z|X')
    if sum(invalid_peptide_index) != 0:
        invalid_peptides = df.iloc[:, 0][invalid_peptide_index].tolist()
        raise Exception("Error: invalid peptides found. %s" % (','.join(invalid_peptides)))

    return True


def valid_exp(df):
    if df.isnull().values.any():
        raise Exception("Error: expression file contains NaN values")

    if df.shape[1] > 1:
        raise Exception("Error: more than one column in expression file.")

    return True


def collect_unaccepted_peptides(input, allele, accept_lens):
    filter_peptides = input[~input.Peptide.str.len().isin(accept_lens)]
    filter_peptides['Allele'] = allele
    filter_peptides['Score'] = "NA"
    filter_peptides['Present'] = "NA"
    return filter_peptides


def run_PSSM(para, allele):
    global package
    global tmpdir

    os.makedirs(tmpdir, exist_ok=True)

    input = para.input
    accept_lens = para.allele_length.get(allele.split("-")[1])

    input = pd.read_table(input, header = None, names = ['Peptide'])
    fail_collection = collect_unaccepted_peptides(input, allele, accept_lens)


    def use_PSSMx(hla, l):
        excluded_hla = ["HLA-" + hla for hla in ["A1101", "A0207", "A0201", "A2402", "A0203", "A0101", "A0301",
                                                 "A3101", "B3501", "B4001", "B4403", "B5101", "B5401", "B0702",
                                                 "A3201", "B1501", "A6802", "A0204", "A2902", "B4402", "B5701"]]

        excluded_hla_lr = {hla: range(9, 12) for hla in excluded_hla}
        excluded_hla_lr['B5101'] = range(8, 11)

        if (hla in excluded_hla):
            if l in excluded_hla_lr.get(hla):
                return False
            else:
                return True
        else:
            return True

    if valid_peptide(input):
        predictions = pd.DataFrame()
        for length in accept_lens:
            input2 = input[input.Peptide.str.len() == length]
            if input2.shape[0] == 0:
                continue

            if not use_PSSMx(allele, length):
                tmp_pep = "{}/{}mer_peptide_{}.txt".format(tmpdir, length, uuid.uuid4())
                input2.Peptide.to_csv(tmp_pep, index = False, header = None)

                peptide = tmp_pep
                prediction_score = np.array([0])

                for i in range(1, 6):
                    pssm_file_resource = 'model/PSSM/stack{}/pssm_file.lst'.format(i)
                    pssm_file = pkg_resources.resource_filename(package, pssm_file_resource)
                    pssm_prediction = PSSM(peptide, length, allele, pssm_file)
                    prediction_score = prediction_score + np.array(pssm_prediction.Score)

                prediction_merge = pd.DataFrame()
                prediction_merge['Peptide'] = pssm_prediction.Peptide.copy()
                prediction_merge['Allele'] = pssm_prediction.Allele.copy()
                prediction_merge['Score'] = list(prediction_score / 5)

                predictions = pd.concat([predictions, prediction_merge])
            else:
                tmp_pep = "{}/{}mer_peptide_{}.txt".format(tmpdir, length, uuid.uuid4())
                input2.Peptide.to_csv(tmp_pep, index=False, header=None)
                peptide = tmp_pep

                pssm_file_resource = 'model/PSSM/PSSMx/pssm_file.lst'
                pssm_file = pkg_resources.resource_filename(package, pssm_file_resource)
                pssm_prediction = PSSM(peptide, length, allele, pssm_file)
                predictions = pd.concat([predictions, pssm_prediction[['Peptide','Allele','Score']]])

        shutil.rmtree(tmpdir)
        return predictions, fail_collection

def PSSM(peptide, length, hla, pssm_file_list):
    global package

    pssm_dict = {}
    with open(pssm_file_list, 'r') as f:
        for line in f:
            line = line.strip()
            tmp = line.split("\t")
            pssm_file_resource = "model/{}".format(tmp[0])
            pssm_file = pkg_resources.resource_filename(package, pssm_file_resource)
            pssm_dict["%s:%s"%(tmp[1],tmp[2])] = pssm_file

    pssm_key = hla + ":" + str(length)
    if pssm_dict.get(pssm_key) != None:
        pssm_result = scorePeptides_PSSM(peptide, hla, length, pssm_dict.get(pssm_key))
    else:
        raise Exception("Error: No PSSM found for HLA %s with length %d" % (hla, length))

    return pssm_result


def scorePeptides_PSSM(peptide, allele, length, pssm):
    lib = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    pssm_dict = defaultdict(dict)

    global package
    global tmpdir

    with open(pssm, 'r') as f:
        for group in get_groups(f, ">"):
            group[0] = group[0].replace(r">","")
            if len(group) != 21:
                raise Exception("Error:Wrong PSSM, number of amino acid unequal to 20.")

            for i in range(1, len(group)):
                pssm_dict[group[0]][lib[i - 1]] = group[i].split("\t")

    pssm_output = "{}/pssm_score_{}.txt".format(tmpdir,uuid.uuid4())
    with open(pssm_output, "w") as fw:
        with open(peptide, 'r') as f:
            for line in f:
                line = line.strip()
                pssm_id = allele + " " + str(length)
                if pssm_dict.get(pssm_id) != None:
                    score = 0
                    for i in range(len(line)):
                        char = line[i]
                        score += np.float64(pssm_dict[pssm_id][char][i])

                    score = score / length
                    fw.write("%s,%s,%.15f\n" % (line, allele, score))
                else:
                    raise Exception("Error: No PSSM found for HLA %s with length %d" % (allele, length))


    pssm_result = pd.read_csv(pssm_output, header = None)
    pssm_result.columns = ['Peptide','Allele','Score']
    os.remove(pssm_output)
    return pssm_result


def get_groups(seq, group_by):
    data = []
    for line in seq:
        line = line.strip()
        if line.startswith(group_by):
            if data:
                yield data
                data = []
        data.append(line)

    if data:
        yield data


def run_mode1(para, allele):
    pssm_output, fail_collection = run_PSSM(para, allele)
    return pssm_output, fail_collection

def run_mode2(para, allele):
    global package
    global collections

    pssm, fail_collection = run_PSSM(para, allele)
    pssm = pssm.reset_index(drop=True)

    ###load length distribution
    length_dist_resource = 'model/allele_length_distribution.pkl'
    length_dist_file = pkg_resources.resource_filename(package, length_dist_resource)
    length_dist = joblib.load(length_dist_file)

    epic_input = build_epic_input(pssm, para, length_dist.get(allele.split(r"-")[1]))
    # allele = para.allele.split(r"-")[1]
    ###load scaler
    scaler_resource = 'model/scaler.pkl'
    scaler_file = pkg_resources.resource_filename(package, scaler_resource)
    scaler = joblib.load(scaler_file)
    ###load predictor
    predictor_resource = 'model/EPIC.pkl'
    predictor_file = pkg_resources.resource_filename(package, predictor_resource)
    predictor = joblib.load(predictor_file)
    scaled_df = scaler.transform(epic_input[['Score','TPM','length']])

    proba = pd.DataFrame(predictor.predict_proba(scaled_df)[:,1],index = epic_input.Peptide, columns=['Score'])
    proba['Allele'] = allele
    return proba, fail_collection


def build_epic_input(pssm, para, allele_length_dist = None):
    exp_file = para.exp
    exp = pd.read_table(exp_file, skip_blank_lines=False, header=None, names = ['TPM'])
    peptide_file = para.input
    peptide = pd.read_table(peptide_file, skip_blank_lines=False, header=None, names=['Peptide'])

    if valid_exp(exp) and valid_peptide(peptide):
        peptide_exp = pd.concat([peptide, exp], axis=1)
        peptide_exp.TPM = np.log2(peptide_exp.TPM + 1)
        epic_input = pd.merge(pssm, peptide_exp, on = 'Peptide')
        epic_input['length'] = epic_input.Peptide.str.len()
        epic_input.length = epic_input.length.replace(allele_length_dist)

        return epic_input

def run(para, allele, mode):
    result = None
    if mode == 1:
        result = run_mode1(para, allele)
    elif mode == 2:
        result = run_mode2(para, allele)

    return result





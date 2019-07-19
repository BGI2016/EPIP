#Copyright (c) 2018, BGI-Shenzhen

from epic.parameter import Parameter
from epic.predict import run
from epic.expand import retrain
import pandas as pd


def start():
    para = Parameter().parse()
    results = pd.DataFrame()
    fail_collections = pd.DataFrame()
    if para.mode == 1 or para.mode == 2:
        print("running...")
        for allele in para.allele.split(","):
            result, fail_collection = run(para, allele, para.mode)
            results = pd.concat([results, result])
            fail_collections = pd.concat([fail_collections, fail_collection])
        print("writing output...")
        write_output(results, fail_collections, para)
        print("writing output finished.")
    elif para.mode == 3:
        print("add %s to EPIC" % (para.allele))
        retrain(para)


def write_output(results, fail_collections, para):
    results.loc[results.Score > para.threshold, 'Present'] = 1
    results.Present = results.Present.fillna(0)
    if para.mode == 2:
        results['Peptide'] = results.index

    output = pd.concat([results[['Peptide', 'Allele', 'Score', 'Present']], fail_collections])

    if para.sort_output:
        output.sort_values('Score', ascending=False).to_csv(para.output, index = False)
    else:
        output.to_csv(para.output, index = False)
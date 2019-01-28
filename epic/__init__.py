#Copyright (c) 2018, BGI-Shenzhen

from epic.parameter import Parameter
from epic.predict import run
from epic.expand import retrain
import pandas as pd


def start():
    para = Parameter().parse()
    results = pd.DataFrame()
    if para.mode == 1 or para.mode == 2:
        print("running...")
        for allele in para.allele.split(","):
            result = run(para, allele, para.mode)
            results = pd.concat([results, result])
        print("writing output...")
        write_output(results, para)
        print("writing output finished.")
    elif para.mode == 3:
        print("add %s to EPIC" % (para.allele))
        retrain(para)


def write_output(df, para):
    df.loc[df.Score > para.threshold, 'Present'] = 1
    df.Present = df.Present.fillna(0)
    if para.mode == 2:
        df['Peptide'] = df.index

    output = df[['Peptide', 'Allele', 'Score', 'Present']]

    if para.sort_output:
        output.sort_values('Score', ascending=False).to_csv(para.output, index = False)
    else:
        output.to_csv(para.output, index = False)
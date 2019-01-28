#Copyright (c) 2018, BGI-Shenzhen

import argparse
from sklearn.externals import joblib
import pkg_resources
import pickle

class Parameter():
    def __init__(self):
        self._parser = argparse.ArgumentParser()
    
    def parse(self):
        self._parser.add_argument("-v", "--version", help = "EPIC version", action = "store_true")
        self._parser.add_argument("-a", "--allele", help = 'HLA alleles, separated by ,. Use "-p" to check out supported HLA alleles')
        self._parser.add_argument("-l", "--length", help = "lengths of peptides, separated by ','.")
        self._parser.add_argument("-f", "--file", help="input file that holds the peptide list")
        self._parser.add_argument("-e", "--exp", help="expression value corresponding to each peptide in the input file."
                                                      "Should be represented by TPM.")
        self._parser.add_argument("-m", "--mode", help="prediction mode. \nMode 1: using PSSM method\nMode 2: using EPIC to "
                                                       "predict based on PSSM score and expression value.\nMode 3:add additional "
                                                       "supported alleles to EPIC Default:1", type = int, choices = [1,2,3],
                                  default=1)
        self._parser.add_argument("-p", "--supportAlleles", help = "HLA alleles supported by EPIC", action = "store_true")
        self._parser.add_argument("-t", "--threshold", help = "threshold for predicting positive outcome.\nMode 1, the threshold ranges "
                                                              "from -0.8 to 0.8.\n Mode 2, the threshold ranges from 0 to 1. Default threshold"
                                                              "is 0 for mode 1 and 0.5 for mode 2.")
        self._parser.add_argument("-s", "--sort", help = "sort peptides according to the predicted score.",
                                   action = "store_true" )
        self._parser.add_argument("-o", "--output", help = "output csv file path")

        args = self._parser.parse_args()

        if args.supportAlleles:
            print("Alleles supported by EPICv1.0 so far:\n"
                  "HLA-A0101\nHLA-A0201\nHLA-A0203\nHLA-A0204\nHLA-A0207\nHLA-A0301\n"
                "HLA-A1101\nHLA-A2301\nHLA-A2402\nHLA-A2501\nHLA-A2601\nHLA-A2901\nHLA-A2902\nHLA-A3001\nHLA-A3004\nHLA-A3101\n"
                "HLA-A3201\nHLA-A6801\nHLA-A6802\nHLA-B0702\nHLA-B0801\nHLA-B1301\nHLA-B1402\nHLA-B1501\nHLA-B1511\nHLA-B15186\n"
                "HLA-B1801\nHLA-B1803\nHLA-B2701\nHLA-B2705\nHLA-B3501\nHLA-B3503\nHLA-B3508\nHLA-B3701\nHLA-B3801\nHLA-B3901\n"
                "HLA-B3906\nHLA-B3924\nHLA-B4001\nHLA-B4002\nHLA-B4101\nHLA-B4402\nHLA-B4403\nHLA-B4501\nHLA-B4901\nHLA-B5001\n"
                "HLA-B5101\nHLA-B5201\nHLA-B5401\nHLA-B5501\nHLA-B5601\nHLA-B5701\nHLA-B7301\nHLA-C0102\nHLA-C0202\nHLA-C0301\n"
                "HLA-C0303\nHLA-C0304\nHLA-C0401\nHLA-C0501\nHLA-C0602\nHLA-C0701\nHLA-C0702\nHLA-C0802\nHLA-C1203\nHLA-C1204\n"
                "HLA-C1402\nHLA-C1502\nHLA-C1505\nHLA-C1601\nHLA-C1701\nHLA-G0101\n")
            exit(1)

        if args.version:
            print("EPIC version 1.0")
            exit(1)

        self.mode = args.mode
        if self.mode == 1 or self.mode == 2:
            self.check_prediction_required_parameters(args)
        elif self.mode == 3:
            self.check_retrain_required_parameters(args)
            return self

        if self.mode == 1 and args.threshold is None:
            self.threshold = self.assign_default_threshold(self.allele, self.mode)
        elif self.mode == 1 and args.threshold is not None:
            threshold = float(args.threshold)
            self.assign_threshold(-0.8, 0.8, threshold)
        elif self.mode == 2 and args.threshold is None:
            self.threshold = self.assign_default_threshold(self.allele, self.mode)
        elif self.mode == 2 and args.threshold is not None:
            threshold = float(args.threshold)
            self.assign_threshold(0, 1, threshold)

        self.len = args.length.split(",")
        self.sort_output = args.sort

        return self

    def check_prediction_required_parameters(self, args):
        self.input = args.file
        self.output = args.output
        self.length = args.length
        self.allele = args.allele

        if self.input is None:
            raise Exception("Error: input file is required")
        elif self.output is None:
            raise Exception("Error: output path is required")
        elif self.allele is None:
            raise Exception("Error: HLA allele is required")
        elif self.length is None:
            raise Exception("Error: peptides length is required")

        if args.mode == 2 and args.exp is None:
            raise Exception("Error: expression file is required by mode 2")
        else:
            self.exp = args.exp

        supported_allels_pkl = pkg_resources.resource_filename("epic", "model/supported_alleles.pkl")
        self.supported_allele_list = joblib.load(supported_allels_pkl)

        for allele in self.allele.split(","):
            if allele not in self.supported_allele_list:
                raise Exception(
                    'Error: HLA allele %s is not supported by EPICv1.0. Use "-p" to check out supported alleles' % (
                        allele))

    def check_retrain_required_parameters(self, args):
        self.input = args.file
        self.allele = args.allele

        if self.input is None:
            raise Exception("Error: input file is required")
        elif self.allele is None:
            raise Exception("Error: HLA allele is required")

        supported_allels_pkl = pkg_resources.resource_filename("epic", "model/supported_alleles.pkl")
        self.supported_allele_list = joblib.load(supported_allels_pkl)
        self.supported_allele_list.append(self.allele)

        supported_allels_pkl = pkg_resources.resource_filename("epic", "model/supported_alleles.pkl")
        pickle.dump(self.supported_allele_list, open(supported_allels_pkl, "wb"))

    def assign_threshold(self, min, max, threshold):
        if threshold > max or threshold < min:
            raise Exception("Error: Invalid threshold for mode 1. Range:(%f, %f)" % (min, max))
        else:
            self.threshold = threshold

    def assign_default_threshold(self, allele, mode):
        mode1_default_threshold = 0

        mode2_default_threshold = 0.5

        if mode == 1:
            return mode1_default_threshold
        elif mode == 2:
            return mode2_default_threshold
        else:
            raise Exception("Unsupported mode.")



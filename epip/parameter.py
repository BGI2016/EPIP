#Copyright (c) 2018, BGI-Shenzhen

import argparse
from sklearn.externals import joblib
import pkg_resources
import pickle
import collections

class Parameter():
    def __init__(self):
        self._parser = argparse.ArgumentParser()


    def parse(self):
        self._parser.add_argument("-v", "--version", help = "EPIC version", action = "store_true")
        self._parser.add_argument("-a", "--allele", help = 'HLA alleles, separated by ,. Use "-p" to check out supported HLA alleles')
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

        self.allele_length = {"A0101":[9,10,11], "A0201":[9,10,11], "A0203":[9,10], "A0204":[9,10,11],
                         "A0207":[9,10,11], "A0301":[9,10,11], "A2402":[9,10,11], "A2902":[9,10],
                         "A3101":[9,10,11], "A6802":[9,10,11], "B3501":[9,10,11], "B4402":[9,10,11],
                         "B4403":[9,10,11], "B5101":[8,9], "B5401":[9,10,11], "B5701":[9,10,11],
                         "A1101":[9,10,11], "A3201":[9], "B0702":[9,10,11], "B1501":[9,10,11], "B4001":[9,10,11],
                         "C0102":[8,9,10], "C0202":[8,9,10], "C0303":[8,9,10], "C0304":[8,9,10], "C0401":[8,9,10],
                         "C0602":[8,9], "C0701":[9], "C0702":[8,9], "C1203":[8,9],
                         "C1402":[8,9,10], "C1502":[8,9], "C1601":[8,9,10], "C1701":[8,9], "A2301":[9,10], "A2501":[9,10],
                         "A2601":[9], "A2901":[9], "A6801":[9,10,11], "B0801":[8,9], "B1402":[9], "B1518":[8,9,10],
                         "B1801":[8,9,10], "B2701":[9,10], "B2705":[9,10,11], "B3503":[9,10,11], "B3508":[9,10],
                         "B3701":[9], "B3801":[9], "B3901":[8,9,10],"B3906":[8,9,10],"B3924":[8,9],"B4002":[9,10,11],
                         "B4101":[9,10], "B4501":[9,10,11], "B4901":[9], "B5001":[9], "B5501":[9], "B5601":[9],"B7301":[9],
                         "A3301":[8,9],"B1302":[9,10],"B1503":[9],"A3303":[9,10,11]}

        if args.supportAlleles:
            print("Alleles supported by EPICv1.0 so far:\n %s\n" % ('\n'.join(['HLA-' + a for a in self.allele_length.keys()])))
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
            self.threshold = self.assign_default_threshold(self.mode)
        elif self.mode == 1 and args.threshold is not None:
            threshold = float(args.threshold)
            self.assign_threshold(-0.8, 0.8, threshold)
        elif self.mode == 2 and args.threshold is None:
            self.threshold = self.assign_default_threshold(self.mode)
        elif self.mode == 2 and args.threshold is not None:
            threshold = float(args.threshold)
            self.assign_threshold(0, 1, threshold)

        self.sort_output = args.sort

        return self


    def check_prediction_required_parameters(self, args):
        self.input = args.file
        self.output = args.output
        self.allele = args.allele

        if self.input is None:
            raise Exception("Error: input file is required")
        elif self.output is None:
            raise Exception("Error: output path is required")
        elif self.allele is None:
            raise Exception("Error: HLA allele is required")


        if args.mode == 2 and args.exp is None:
            raise Exception("Error: expression file is required by mode 2")
        else:
            self.exp = args.exp

        supported_allels_pkl = pkg_resources.resource_filename("epip", "model/supported_alleles.pkl")
        self.supported_allele_list = joblib.load(supported_allels_pkl)

        for allele in self.allele.split(","):
            if allele not in self.supported_allele_list:
                print(
                    'Warning: HLA allele %s is not supported by EPICv1.0. Use "-p" to check out supported alleles' % (
                        allele))


    def check_retrain_required_parameters(self, args):
        self.input = args.file
        self.allele = args.allele

        if self.input is None:
            raise Exception("Error: input file is required")
        elif self.allele is None:
            raise Exception("Error: HLA allele is required")

        supported_allels_pkl = pkg_resources.resource_filename("epip", "model/supported_alleles.pkl")
        self.supported_allele_list = joblib.load(supported_allels_pkl)
        self.supported_allele_list.append(self.allele)

        supported_allels_pkl = pkg_resources.resource_filename("epip", "model/supported_alleles.pkl")
        pickle.dump(self.supported_allele_list, open(supported_allels_pkl, "wb"))


    def assign_threshold(self, min, max, threshold):
        if threshold > max or threshold < min:
            raise Exception("Error: Invalid threshold for mode 1. Range:(%f, %f)" % (min, max))
        else:
            self.threshold = threshold


    def assign_default_threshold(self, mode):
        mode1_default_threshold = 0

        mode2_default_threshold = 0.5

        if mode == 1:
            return mode1_default_threshold
        elif mode == 2:
            return mode2_default_threshold
        else:
            raise Exception("Unsupported mode.")


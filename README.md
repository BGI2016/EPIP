# EPIC
EPIC is standing for Epitope Presentation Integrated prediCtion, a tool aims to find out tumor specific neoantigens

#Copyright (c) 2019, BGI-Shenzhen

__EPIC version 1.0__

__Installation:__

1. The software runs in the python3.x environment. Download the source code and decompress it, go into the decompressed directory and run the command bwlow:

```pip install .```

2. Make sure that PERL is installed and added to the PATH environmental vairable.

__EPIC supports 3 modes by now:__  

1. mode1: predict epitope presentation based on PSSM only  
2. mode2: predict epitope presentation using the full EPIC model, which integrates PSSM, peptide expression and length  
3. mode3: add addition alleles to the EPIC supported allele repository. EPIC will build PSSM and obtain length distribution upon the provided peptide list and the corresponding allele.

__Before using the software, try running the example codes below. Normally, running the codes should output nothing:__  

```
epic-predict -m 1 -l 9,10,11 -a HLA-A1101 -f EPIC/test/test_peptide.txt -o EPIC/test/my_test_peptide_mode1_A1101_myout.txt
diff EPIC/test/my_test_peptide_mode1_A1101_myout.txt EPIC/test/test_peptide_mode1_A1101_output.txt
```
```
epic-predict -m 2 -l 9,10,11 -a HLA-A1101 -f EPIC/test/test_peptide.txt -e EPIC/test/test_peptide_exp.txt -o EPIC/test/my_test_peptide_mode2_A1101_myout.txt\n
diff EPIC/test/my_test_peptide_mode2_A1101_myout.txt EPIC/test/test_peptide_mode2_A1101_output.txt
```
```
epic-predict -m 2 -l 9,10,11 -a HLA-B0801 -f EPIC/test/test_peptide.txt -e EPIC/test/test_peptide_exp.txt -o EPIC/test/my_test_peptide_mode2_B0801_myout.txt\n
diff EPIC/test/my_test_peptide_mode2_B0801_myout.txt EPIC/test/test_peptide_mode2_B0801_output.txt
```
__Citation__  
If you find EPIC is useful to your research or work, please cite <article name>

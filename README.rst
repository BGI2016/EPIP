#Copyright (c) 2018, BGI-Shenzhen

EPIC version 2.1

Changes:

1. EPIC v2.1 can predict multiple lengths of peptides at the same time
2. EPIC v2.1 now supports prediction for up to 72 alleles, use "epic-predict -p" to see more details.

Installation:

The software runs in the python3.x environment. Download the source code and uncompress it, go into the 'epic' directory
 and run the command bwlow:

pip install .

Try running the example code below to make sure EPIC runs properly:

epic-predict -m 1 -l 9,10,11 -a HLA-A1101 -f EPIC/test/test_peptide.txt -o EPIC/test/my_test_peptide_mode1_output.txt
diff EPIC/test/my_test_peptide_mode2_output.txt EPIC/test/test_peptide_mode1_output.txt

epic-predict -m 2 -l 9,10,11 -a HLA-A1101 -f EPIC/test/test_peptide.txt -e EPIC/test/test_peptide_exp.txt -o EPIC/test/my_test_peptide_mode2_output.txt\n
diff EPIC/test/my_test_peptide_mode2_output.txt EPIC/test/test_peptide_mode2_output.txt

epic-predict -m 3 -l 9,10,11 -a HLA-B0801 -f EPIC/test/test_peptide.txt -e EPIC/test/test_peptide_exp.txt -o EPIC/test/my_test_peptide_mode3_B0801_output.txt\n
diff EPIC/test/my_test_peptide_mode3_B0801_output.txt EPIC/test/test_peptide_mode3_B0801_output.txt

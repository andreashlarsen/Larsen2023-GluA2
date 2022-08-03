# 
# MODELLER LOG
#
# use scripts from: 
# https://salilab.org/modeller/wiki/Missing%20residues
#
# IMPORT MODELLER
module load modeller/9v16/64
# STEP 1
#
# modify modeller_step1.py
mod9.16 modeller_step1.py
# this script generates *.seq
#
# now generate alignment file 
cp *.seq alignment.ali
# follow instructions on webpage (and use missing residue info in pdb file) to make alignment.ali
# https://salilab.org/modeller/wiki/Missing%20residues
#
# STEP 2
# modify modeller_step2.py
mod9.16 modeller_step2.py
#

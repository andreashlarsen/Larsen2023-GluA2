# extend by 10,000,000 ps = 10,000 ns = 10 us
gmx convert-tpr -s md.tpr -extend 10000000 -o md_extend.tpr -quiet

# overwrite old md.tpr
mv md_extend.tpr md.tpr

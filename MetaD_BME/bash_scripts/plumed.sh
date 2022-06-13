# for getting dist and weights, you need a mass file --mc NAME
# you get that by running a (short) sim with DUMPMASSCHARGE in the plumed file: 
# https://www.plumed.org/doc-v2.6/user-doc/html/_d_u_m_p_m_a_s_s_c_h_a_r_g_e.html
gmx mdrun -s md.tpr -nsteps 1 -plumed plumed_mc.dat -quiet

# get distance and weights
plumed driver --mf_xtc md_cent_dt4000_nojump.xtc --mc mcfile --plumed plumed_weight.dat

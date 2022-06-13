# flow file for running full-length glua2 sims
# in martini 3 beta
#
# 1) remove quisqalate and GSG1L from pdb 5vhz (5hvz_clean)
# 1) fix missing residues and missing atoms with modeller (5hvz_fill.pdb) - do not include whole N-term loop, only until res K823 
# 2) align to OPM 5vhz, so it is  placed right in bilayer (5vhz_fill_align_opm.pdb)
# 3) cg structure
# 4) place in box with bilayer

## ensure modules and other settings from batchrc file
source ~/.bashrc

# parallelisation with mdrun
pinoffset=0
noCPUs=4
noGPUs=1

# nb=auto to use gpu, nb=cpu to avoid use of gpu (if used by other jobs)
nb=auto

# define paths to scripts
py2=/sansom/s157/bioc1642/anaconda2/bin/python2.7
py3=/sansom/s157/bioc1642/anaconda3/bin/python3.7
martinize=/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/22/martini_3.0.b.3.2/martinize
insane=/sansom/s157/bioc1642/Desktop/Scripts/insane.py
cg2at=/sansom/s157/bioc1642/Desktop/Scripts/cg2at/cg2at-gmx2019-martini2.pl

# define paths to mdp files
min=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/minimization.mdp
eq=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/equilibration.mdp
prod=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/production.mdp

# define path to pdb file
protein=/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/Modeller_5vhz_short/5vhz_fill_align_opm.pdb
cp $protein . 

# define path to topology files
ffdir=/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/22/martini_3.0.b.3.2

# define path to plumed file
plumed_file=/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/22/metad.dat
cp $plumed_file .

## alter plumed file: smaller HEIGHT, longer PACE (more gentle metad)
sed -i -e 's/metad: METAD ARG=dist PACE=1000 HEIGHT=0.8 SIGMA=0.1 GRID_MIN=0 GRID_MAX=18 FILE=HILLS/metad: METAD ARG=dist PACE=5000 HEIGHT=0.4 SIGMA=0.1 GRID_MIN=0 GRID_MAX=18 FILE=HILLS/g' metad.dat

## CG protein  with martinize 3beta. -merge flag adjust elastic network
$py2 $martinize -f ${protein} -x CG.pdb -o topol.top -v -ff ${ffdir}/martini303v.partition -elastic -merge A,B -merge C,D -dssp dssp

## Build bilayer + protein with insane - no antifreeze water for martini3
box_xy=32
box_z=30
$py2 $insane -f CG.pdb -o Bilayer.gro -p topol.top -x $box_xy -y $box_xy -z $box_z -l POPC:100 -sol W -salt 0

## change typology file for system
sed -i -e 's/W  /WN /g' topol.top
sed -i -e 's/NA+ /TNA /g' topol.top
sed -i -e 's/#include "martini.itp"//g' topol.top
sed -i 's/Protein        1/Protein_A+Protein_B 1\nProtein_C+Protein_D 1/g' topol.top
cat << EOF > topol.add
#include "$ffdir/martini_v3.0.b.3.2.itp"
#include "$ffdir/martini_v3.0_ions.itp"
#include "$ffdir/martini_v3.0_phospholipids.itp"
#include "$ffdir/martini_v3.0_solvents.itp"
#include "Protein_A+Protein_B.itp"
#include "Protein_C+Protein_D.itp"

#ifdef POSRES
#include "posre.itp"
#endif
EOF

cat topol.add topol.top > tmp
rm topol.add
mv tmp topol.top

## change gro file
sed -i -e 's/W /WN/g' Bilayer.gro
sed -i -e 's/ W/WN/g' Bilayer.gro
sed -i -e 's/NA+ /NA  /g' Bilayer.gro # for some reason, it should be TNA i .top, but NA in .gro 
sed -i -e 's/ NA+/  NA/g' Bilayer.gro

## minimization
gmx grompp -f $min -c Bilayer.gro -p topol.top -o min.tpr -quiet
gmx mdrun -deffnm min -quiet -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb $nb -nsteps 1000

## check energies
echo 8 0 | gmx energy -f min.edr -o potential.xvg -quiet
#xmgrace potential.xvg

cat << EOF > index.input
14 | 15
name 20 SOL
13
name 21 LIP
q
EOF

## make index file
gmx make_ndx -f Bilayer.gro -quiet < index.input

## see overview of index file
# gmx make_ndx -f Bilayer.gro -n index.ndx -quiet

## generate position restraint for protein
echo 1 | gmx genrestr -f CG.pdb -fc 1000 1000 1000 -o posre.itp -quiet

## reduce time step in equilibration 
cp $eq .
sed -i -e 's/dt                       = 0.03/dt                       = 0.02/g' equilibration.mdp
eq=equilibration.mdp

## equilibrate, protein restrained, no metaD
gmx grompp -f $eq -c min.gro -r min.gro -p topol.top -o eq.tpr -quiet -n index.ndx
gmx mdrun -deffnm eq -v -quiet -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb $nb -nsteps 20000

## check temperature
echo 12 0 | gmx energy -f eq.edr -o temperature.xvg -quiet
#xmgrace temperature.xvg

## check pressure
echo 13 0 | gmx energy -f eq.edr -o pressure.xvg -quiet
#xmgrace pressure.xvg

## check density
echo 19 0 | gmx energy -f eq.edr -o density.xvg -quiet
#xmgrace density.xvg 

## reduce timestep from 35 fs to 20 fs (got too many LINCS warnings)
cp $prod prod.mdp
sed -i -e 's/dt                       = 0.035        ; 35 fs/dt                       = 0.020        ; 20 fs/g' prod.mdp

## production run (very short testrun), with metaD
gmx grompp -f prod.mdp -c eq.gro -r eq.gro -p topol.top -o md.tpr -quiet -n index.ndx
gmx mdrun -deffnm md -v -quiet -plumed metad.dat -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -nb $nb -nsteps 10000

## longer production run (performance 1 GPU, 4 CPU, pie1: ~600 ns/day, sams: ~400 ns/day)
continue=/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/31/continue.sh
cp $continue .
bash continue.sh

## calc free energy and evaluate convergence
#plumed sum_hills --hills HILLS
#plumed sum_hills --hills HILLS --stride 10000 --mintozero --min 0 --max 20 --bin 200

## center and fix pbc for visualization
#echo 1 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -pbc mol -ur compact -center -o md_cent.xtc -quiet

## center, fix pbc and reduce file size (original trajectory has 1 frame per 10 ps)
#echo 1 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -pbc mol -ur compact -center -dt 400 -o md_cent_dt400.xtc -quiet # 1 frame/400 ps
#echo 1 0 | gmx trjconv -f md_cent_dt400.xtc -s md.tpr -n index.ndx -pbc mol -ur compact -center -dt 1000 -o md_cent_dt1000.xtc -quiet # 1 frame/1000 ps
#echo 1 0 | gmx trjconv -f md_cent_dt1000.xtc -s md.tpr -n index.ndx -pbc mol -ur compact -center -dt 4000 -o md_cent_dt4000.xtc -quiet # 1 frame/4000 ps 
#echo 1 0 | gmx trjconv -f md_cent_dt4000.xtc -s md.tpr -n index.ndx -pbc mol -ur compact -center -dt 40000 -o md_cent_dt40000.xtc -quiet # 1 frame/40000 ps 
echo 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -pbc nojump -dt 1000 -o md_cent_dt1000_nojump.xtc -quiet # 1 frame/1000 ps, keep tetramer whole
echo 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -pbc nojump -dt 4000 -o md_cent_dt4000_nojump.xtc -quiet # 1 frame/4000 ps, keep tetramer whole
echo 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -pbc nojump -dt 40000 -o md_cent_dt40000_nojump.xtc -quiet # 1 frame/40000 ps, keep tetramer whole

## prepare BME
# rm WEIGHT # else, plumed.sh will append to this file
#bash plumed.sh #generate file: WEIGHT, with distances and log(weights), check what *.xtc file is used - how many frames to use?
#python plot_weigths.py # plot distances and weights and generate files: weight and distance, with weights (not log) and distances respectively
#check for zeros and infinities in weight file, adjust plot_weights accordingly
#check number of frames with: gmx check -f md_cent_dt4000_nojump.xtc -quiet
#bash cg2at.sh # run several instances in parallel 
#python pepsi.py
# may be a good idea to use taskset "--cpu-list 2-8", however, pepsi.py was faster without (5.5 sec for 10 reps vs 6.1 sec)
#python rg.py # generate file: Rg_pepsi.dat, with Rg

## run BME
#mkdir -p BME
#cd BME
#python run_BME.py
#python plot_BME.py (plot with chi2r vs S)
#select theta_opt
#python repeat_BME.py (plot with Rg histogram and plot with fit before/after BME)

## clean up
rm \#*
rm bck.*
cd ..	
rm \#*
rm bck.*


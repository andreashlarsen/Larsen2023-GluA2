source ~/.bashrc

## generate cg2at folder
folder=cg2at
mkdir -p $folder


first=2001
#last=2500
#first=1
last=2001
for i in $(seq $first $last)
#for i in 1 # testing
do

echo -------------------------------------
echo  PREPARE FOR CG2AT, frame $i of $last 
echo -------------------------------------

## extract ith frame from trajectory
cat << EOF > frame$i.ndx
[ frames ]

$i

EOF
#echo 1 | gmx trjconv -f md_cent_dt40000_nojump.xtc -s md.tpr -fr frame$i.ndx -o ${folder}/frame$i.pdb -quiet
echo 1 | gmx trjconv -f md_cent_dt4000_nojump.xtc -s md.tpr -fr frame$i.ndx -o ${folder}/frame$i.pdb -quiet

## delete index file 
rm frame$i.ndx

## go to directory
cd cg2at

echo -------------------------------------
echo  CG2AT, frame $i of $last
echo -------------------------------------

## delete previous CG2AT folder (if it exists)
rm -r CG2AT$i

## CG2AT - use old version of cg2at, which works better with martine 3beta
#/sansom/s157/bioc1642/Desktop/Scripts/cg2at_owen/cg2at/cg2at -c frame$i.pdb -ff charmm36-mar2019-updated -fg martini_3-0_charmm36 -w tip3p -silent -ncpus 2 -loc CG2AT$i
/sansom/s157/bioc1642/Desktop/Scripts/cg2at_owen/old_versions/cg2at_v1/cg2at/cg2at -c frame$i.pdb -ff charmm36-mar2019-updated -fg martini_3-0_charmm36 -w tip3p -silent -ncpus 1 -loc CG2AT$i

## delete ith frame
rm frame$i.pdb

## go to parent directory
cd ..

done


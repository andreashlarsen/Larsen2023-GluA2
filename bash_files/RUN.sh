#### steps:

#####################
#### PREPARATION ####
#####################

## ensure modules and other settings from batchrc file
source ~/.bashrc

## choose folder (specific settings)
folder=35 # see file "folders"
rep=1 # repeat number (for umbrella sampling) 

## set performance parameters
NCPU=4
pinoffset=0

## Script options
DOWNLOAD=0
PREPARE=0
PULL=0
ANALYSIS=1
HBOND=0

AXON=1

## fix missing atoms with modeller (see modeller folder) and do mutations in pymol, 
## truncate to only include chain B and D. Output:
if [ $folder -eq 16 ] || [ $folder -eq 99 ]
then
    pdb=4u2p_NTD_BD_H205A_2.pdb
elif [ $folder -eq 17 ]
then
    pdb=4u2p_NTD_BD_H205A_1.pdb
elif [ $folder -eq 35 ] || [ $folder -eq 36 ]
then
    pdb=5vhz_NTD_chainBD.pdb
else
    pdb=4u2p_NTD_BD.pdb
fi

echo pdb = $pdb
echo rep = $rep

## scripts and programs
py3=/sansom/s157/bioc1642/anaconda3/bin/python3.7
cp /sansom/s157/bioc1642/Desktop/Scripts/umbrella/extract_frames.py .
sed -i -e "s/step_size = 5  /step_size = 10  /g" extract_frames.py # change step size
sed -i -e "s/max_dist = 7.0/max_dist = 10.0/g" extract_frames.py # change max dist

##################################################
#### download once for all protonation states ####
##################################################

if [ $DOWNLOAD -eq 1 ]
then

    ## copy charmm36m ff params to folder (from McKarrels website)

    ## copy ions.mdp from http://www.mdtutorials.com/gmx/complex/Files/ions.mdp to the folder
    wget http://www.mdtutorials.com/gmx/complex/Files/ions.mdp

    ## copy min.mdp from http://www.mdtutorials.com/gmx/complex/Files/em.mdp to the folder
    wget http://www.mdtutorials.com/gmx/complex/Files/em.mdp
    mv em.mdp min.mdp

    ## copy nvt.mdp from http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp to the folder and change one tc group:
    wget http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp
    sed -i -e "s/tc-grps                 = Protein_JZ4/tc-grps                 = Protein/g" nvt.mdp

    ## copy npt.mdp from http://www.mdtutorials.com/gmx/complex/Files/npt.mdp to the folder and change one tc group:
    wget http://www.mdtutorials.com/gmx/complex/Files/npt.mdp
    sed -i -e "s/tc-grps                 = Protein_JZ4/tc-grps                 = Protein/g" npt.mdp

    ## copy md.mdp to the folder and change one tc group:
    wget http://www.mdtutorials.com/gmx/complex/Files/md.mdp
    sed -i -e "s/tc-grps                 = Protein_JZ4/tc-grps                 = Protein/g" md.mdp

    ## run pdb2gmx with default params, and chech His protonation states
    echo 1 | gmx pdb2gmx -f $pdb -o 4u2p_NTD_BD.gro -water tip3p -quiet
    # 37HIS: HSE; 84HIS: HSE; 98HIS: HSE; 205HIS: HSE; 210HIS: HSE; 265HIS: HSE; same in chain 2
fi

###############################################
#### run for individual protonation states ####
###############################################

## settings for simulation
umb=../../umbrella.mdp
mkdir -p $folder
mkdir -p ${folder}/rep${rep}

if [ "$folder" -lt 3 ]
then
	HIS205=$folder
	boxsize=1
elif [ "$folder" -lt 6 ]
then
	HIS205=$((folder - 3))
	boxsize=4
elif [ "$folder" -lt 9 ]
then
	HIS205=$((folder - 6))
	## choose if it is umbrella run
        UMBRELLA=1
elif [ "$folder" -lt 13 ]
then
	HIS205=$((folder - 10))
	UMBRELLA=1
	cp umbrella.mdp umbrella_long.mdp
	sed -i -e 's/nsteps                  = 5000000/nsteps                  = 10000000/g' umbrella_long.mdp
	umb=../umbrella_long.mdp
elif [ "$folder" -lt 16 ]
then
	HIS205=$((folder - 13))
	UMBRELLA=1
	cp umbrella.mdp umbrella_long.mdp
	sed -i -e 's/nsteps                  = 5000000/nsteps                  = 20000000/g' umbrella_long.mdp
        umb=../umbrella_long.mdp
elif [ "$folder" -lt 18 ] || [ $folder -eq 99 ]
then
	HIS205=1
	UMBRELLA=1
elif [ $folder -eq 35 ]
then
	HIS205=1
	UMBRELLA=1
elif [ $folder -eq 36 ]
then
        HIS205=2
        UMBRELLA=1
fi

echo HIS205   = $HIS205
echo boxsize  = $boxsize
echo UMBRELLA = $UMBRELLA

## axon directory
axon_dir=axon_umbrella_prot${folder}_rep${rep}_glua2_ph

## make input for protonation
cat << EOF > prot.input
1
1
1
1
$HIS205
1
1
1
1
1
$HIS205
1
1
EOF

# setup, minimization, equilibration.
if [ $PREPARE -eq 1 ]
then
    if [ $folder -eq 16 ] || [ $folder -eq 17 ] || [ $folder -eq 99 ]
    then
        gmx pdb2gmx -f $pdb -o $folder/rep$rep/prot.gro -water tip3p -quiet < prot.input
    else
        gmx pdb2gmx -f $pdb -o $folder/rep$rep/prot.gro -water tip3p -his -quiet < prot.input
    fi

    ## move to folder
    mv *.itp ${folder}/rep${rep}
    mv *.top ${folder}/rep${rep}
    cd ${folder}/rep${rep}

    ## change topology so it finds ff params in parent directory (i.e. change "./..." to "../...")
    sed -i -e 's/#include ".\/charmm36/#include "..\/..\/charmm36/g' topol.top

    ## make box
    if [ $UMBRELLA -eq 1 ]
    then
        gmx editconf -f prot.gro -o prot_newbox.gro -c -box 25 15 15 -quiet
    else
	gmx editconf -f prot.gro -o prot_newbox.gro -c -d $boxsize -bt cubic -quiet
    fi

    ## solvate
    gmx solvate -cp prot_newbox.gro -cs spc216.gro -o prot_solv.gro -p topol.top -quiet

    ## add ions
    gmx grompp -f ../../ions.mdp -c prot_solv.gro -p topol.top -o ions.tpr -quiet
    echo 13 | gmx genion -s ions.tpr -o prot_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -quiet

    ## minimisation
    gmx grompp -f ../../min.mdp -c prot_solv_ions.gro -p topol.top -o em.tpr -quiet
    gmx mdrun -v -deffnm em -quiet -pin on -pinoffset $pinoffset -ntmpi 1 -ntomp $NCPU -nb auto

    ## check energies
    echo 11 0 | gmx energy -f em.edr -o potential.xvg -quiet
    #xmgrace potential.xvg

    ## nvt equilibration
    gmx grompp -f ../../nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -quiet
    gmx mdrun -deffnm nvt -v -quiet -pin on -pinoffset $pinoffset -ntmpi 1 -ntomp $NCPU -nb auto -nsteps 10000

    ## check temperature
    echo 16 0 | gmx energy -f nvt.edr -o temperature.xvg -quiet
    #xmgrace temperature.xvg

    ## npt equilibration
    gmx grompp -f ../../npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -quiet
    gmx mdrun -deffnm npt -v -quiet -pin on -pinoffset $pinoffset -ntmpi 1 -ntomp $NCPU -nb auto -nsteps 10000

    ## check temperature
    echo 17 0 | gmx energy -f nvt.edr -o pressure.xvg -quiet
    #xmgrace pressure.xvg

    ## check density
    echo 23 0 | gmx energy -f npt.edr -o density.xvg -quiet
    #xmgrace density.xvg 

    cd ../..
fi

cd ${folder}/rep${rep}

## make input for index file (chain A and B)
if [ $folder -eq 35 ] || [ $folder -eq 36 ] 
then
cat << EOF > index.input
ri 1-376
name 18 ChainA
ri 377-752
name 19 ChainB
q
EOF
else
cat << EOF > index.input
ri 1-374
name 18 ChainA
ri 375-748
name 19 ChainB
q
EOF
fi

if [ $UMBRELLA -eq 1 ]
then

    ## pull simulation
    if [ $PULL -eq 1 ]
    then

        # generate index file
        gmx make_ndx -f npt.gro -o index.ndx -quiet < index.input

        #check index file
        #gmx make_ndx -f npt.gro -n index.ndx -o index.ndx -quiet

        ## pull run
        gmx grompp -f ../../pull.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o pull.tpr -quiet
        gmx mdrun -deffnm pull -v -cpi pull.cpt -px pull_pullx.xvg -pf pull_pullf.xvg -quiet -pin on -pinoffset $pinoffset -ntmpi 1 -ntomp $NCPU -nb auto -nsteps 5000000

        # extract frames
        $py3 ../../extract_frames.py # generates extract_frames.ndx
        echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -n -fr extract_frames.ndx -o extract_frames.xtc -quiet

        ## get number of frames (3 lines in file are not frames numbers)
        Nframes=-3; while read -r LINE; do (( Nframes++ )); done < extract_frames.ndx
        echo Nframes  = $Nframes

        # prepare for umbrella sampling, one frame at a time
        for frame in $(seq 1 $Nframes)
        do
            # make index file for one frame
	    cat << EOF > frames_step.ndx
 [ frames ]

$frame

EOF
            # extract one frame
            echo 0 | gmx trjconv -f extract_frames.xtc -s pull.tpr -fr frames_step.ndx -o step_${frame}.gro -quiet

            ## grompp umbrella production runs
            gmx grompp -f $umb -c step_${frame}.gro -r step_${frame}.gro -n index.ndx -p topol.top -o umbrella_step${frame}.tpr -maxwarn 1 -quiet

            ## add lines to summary files for WHAM
            echo "umbrella_step${frame}.tpr" >> tpr_files_umbrella.dat
            echo "umbrella_step${frame}_pullf.xvg" >> pullf_files_umbrella.dat 
        done

        ## run umbrella sampling on axon 
        if [ $AXON -eq 1 ]
        then
            ## make axon directory
            mkdir -p $axon_dir
	
	    ## alter axon flow file and move to axon dir
	    flow=Flow_UMBRELLA.sh
	    cp ../../$flow .
            sed -i -e "s/job-name=GLUA2_FOLDER/job-name=GLU_${folder}_${rep}/g" $flow
            Nframes=-3; while read -r LINE; do (( Nframes++ )); done < extract_frames.ndx
            sed -i -e "s/array=1-NFRAMES/array=1-$Nframes/g" $flow 
            mv $flow $axon_dir

	    ## move umbrella tpr files to axon dir
            mv umbrella_step*.tpr $axon_dir
        
	    ## copy to axon
            scp -r $axon_dir axon:/home/bioc1642/ 
        fi    
     
    fi

    ## do analysis, WHAM
    if [ $ANALYSIS -eq 1 ]
    then
	if [ $AXON -eq 1 ]
	then
            ## copy from axon
	    cp *_files_umbrella.dat $axon_dir
	    cd $axon_dir
            scp axon:/home/bioc1642/$axon_dir/umbrella*pullf.xvg .
        fi
        ## run WHAM (weighted histogram analysis method)
        #skip=$(( time*1000/5*0 )) # skip first fifth of frames (equilibration) for WHAM
        taskset --cpu-list 4 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -o -bsres -temp 323 -quiet #-b $skip #-nBootstrap 200 # outcomment "-nBootstrap 200" to get uncertainties  

        # check PMF curve
        #xmgrace profile.xvg
        #xmgrace bsResults.xvg
        #xmgrace bsProfs.xvg
        #xmgrace -nxy histo.xvg

	## go back to parent directory
	cd ..
    fi
    
    ## calculate H-bonds
    if [ $HBOND -eq 1 ]
    then
        ## get number of frames (3 lines in file are not frames numbers)
        Nframes=-3; while read -r LINE; do (( Nframes++ )); done < extract_frames.ndx
        echo Nframes  = $Nframes
        for step in $(seq 1 $Nframes)
        do
            cd $axon_dir
            if [ $AXON -eq 1 ]
            then
                scp axon:/home/bioc1642/$axon_dir/umbrella_step$step.xtc .
            fi
            echo 18 19 | taskset --cpu-list 2-6 gmx hbond -f umbrella_step$step.xtc -s umbrella_step$step.tpr -n ../index.ndx -quiet -num hbnum$step.xvg
            if [ $AXON -eq 1 ]
            then
                rm umbrella_step$step.xtc
            fi
            cd ..
        done
    fi
else
    ## test production run
    gmx grompp -f ../../md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -quiet    
    gmx mdrun -deffnm md -cpi md.cpt -v -quiet -pin on -pinoffset $pinoffset -ntmpi 1 -ntomp $NCPU -nb auto -nsteps 5000

    if [ $AXON -eq 1 ]
    then
        # generate axon folder
        axon_folder=glua2_${folder}_axon
	mkdir -p $axon_folder
	cp md.* $axon_folder

	# upload
	scp -r $axon_folder axon:/home/bioc1642

	# download (when finished)	
	scp -r axon:/home/bioc1642/$axon_folder .
    else
	## production, continuation (500000 = 1 ns, 50000000 = 100 ns, 100000000 = 200 ns)
	gmx mdrun -deffnm md -cpi md.cpt -v -quiet -pin on -pinoffset $pinoffset -ntmpi 1 -ntomp $NCPU -nb auto -nsteps 50000000
    fi 
    
    ## reduce size of mdp file and center
    echo 1 0 | gmx trjconv -s md.tpr -f md.xtc -o md_cent.xtc -pbc mol -center -ur compact -quiet -dt 100
fi

## clean up
rm \#*
rm step*
cd ..
rm \#*
cd ..
rm \#*

# extra:
#
# Hbonds
#
#for step in {28..38}
#do
#scp axon:/home/bioc1642/axon_umbrella_prot6_glua2_ph/umbrella_step$step.xtc .
#echo 18 19 | taskset --cpu-list 3-6 gmx hbond -f umbrella_step$step.xtc -s umbrella_step$step.tpr -n ../index.ndx -quiet -num hbnum$step.xvg
#rm umbrella_step$step.xtc
#done

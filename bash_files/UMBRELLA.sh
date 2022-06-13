#!/bin/sh

## ensure modules are loades etc from batch file
source ~/.bashrc

## set parameters

# protein folder
folder=6

# activate/deactivate modules of script
PULL=1
PUSH=1
PUSH_EXTRA=1
COPYPULL=0
GROMPP=1
MDRUN=0
AXON=1
ANALYSIS=0

# remove constraint on PIP2 in production runs
NO_CONSTRAINT=1

# define folder to copy preparation files from (if COPYPULL)
pull_folder=/sansom/s157/bioc1642/Desktop/prj_C2/25/PC_15/umbrella_t50_l5_k2000_n0

# define paths to mdp files
umb=/sansom/s157/bioc1642/Desktop/Scripts/umbrella/umbrella.mdp
pull=/sansom/s157/bioc1642/Desktop/Scripts/umbrella/pull.mdp

# timestep in umbrella sampling
dt=1000

# tolerence in WHAM (default is 1e-6)
tol=1e-6

# define paths to scripts 
py3=/sansom/s157/bioc1642/anaconda3/bin/python3.7
extract_frames=/sansom/s157/bioc1642/Desktop/Scripts/umbrella/extract_frames.py

# parallelisation with mdrun
pinoffset=0
noCPUs=4
noGPUs=1

# cpu for cpu only, auto for using gpu
nb=auto

# parameters for push and pull and umbrella 
step_size=5
# 50 = 0.5 nm
# 20 = 0.2 nm
# 10 = 0.1 nm
#  5 = 0.05 nm
#  2 = 0.02 nm
#  1 = 0.01 nm
umb_rep=1
spring_const=2000
time=1000
nsteps=$((time*1000*1000/dt))

# folder name
if [ $lip -eq 1 ]
then
  Folderprefix=PCPSPIP2PIP3
elif [ $lip -eq 2 ]
then
  Folderprefix=PCPSPIP2
elif [ $lip -eq 3 ]
then
  Folderprefix=PCPS
elif [ $lip -eq 4 ]
then
  Folderprefix=PC
fi
    
############ GENERATE FOLDERS ETC #######################################################
    
## define folder name
folder=${folder_prot}/${Folderprefix}_$rep
cd $folder

## go to folder (directory)
dir=umbrella_t${time}_l${step_size}_k${spring_const}_n${umb_rep}
mkdir -p $dir
cd $dir
    
############ PREPARE TO RUN SIMULATION ###################################################
      
## extract key frame from trajectory
cat << EOF > frame_index.ndx
[ frames ]

$key_frame

EOF
echo 0 | gmx trjconv -f ../md.xtc -s ../md.tpr  -o key_frame.gro -fr frame_index.ndx -quiet
       
## make index file
if [ $lip -eq 2 ]
then
  cat << EOF > index.input
r W WF ION
name 19 SOL
r POPC POPS POP2
name 20 LIP
q
EOF
elif [ $lip -eq 3 ]
then
  cat << EOF > index.input
r W WF ION
name 18 SOL
r POPC POPS
name 19 LIP
q
EOF
elif [ $lip -eq 4 ]
then
  cat << EOF > index.input
r W WF ION
name 17 SOL
r POPC
name 18 LIP
q
EOF
fi
gmx make_ndx -f ../md.gro -quiet < index.input
rm index.input

if [ $PULL -eq 1 ]
then
  ## make pull simulation
  cp $pull pull_slow.mdp
  sed -i -e "s/pull_coord1_rate        = 0.0001/pull_coord1_rate        = 0.00002/g" pull_slow.mdp
  gmx grompp -f pull_slow.mdp -c key_frame.gro -r key_frame.gro -p ../topol.top -n index.ndx -o pull.tpr -quiet -maxwarn 1
  gmx mdrun -deffnm pull -v -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet # it will run until it crashes (after roughly 20 min)
    
  ## check position vs time and force vs time
  #xmgrace pull_pullx.xvg
  #xmgrace pull_pullf.xvg  
   
  ## extract frames from steered MD trajectory 
  cp $extract_frames .
  sed -i -e "s/step_size = 5  /step_size = ${step_size}  /g" extract_frames.py # change step size
  $py3 extract_frames.py # generates extract_frames.ndx
  echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -n  -fr extract_frames.ndx -o extract_frames.xtc -quiet
	
  ## check number of extracted frames
  #gmx check -f extract_frames.xtc -quiet
fi

if [ $PUSH -eq 1 ]
then
  ## make push simulation
  cp $pull push_slow.mdp
  sed -i -e "s/pull_coord1_rate        = 0.0001/pull_coord1_rate        = -0.00002/g" push_slow.mdp
  gmx grompp -f push_slow.mdp -c key_frame.gro -r key_frame.gro -p ../topol.top -n index.ndx -o push.tpr -quiet -maxwarn 1
  
  if [ $PUSH_EXTRA -eq 1 ]
  then
    gmx mdrun -deffnm push -v -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet -nsteps 2500000
    # nsteps of 5000000 corresponds to 150 ns, and 3.0 nm with rate of 0.00002 nm/ps and dt of 30 fs.
    # nsteps of 3333333 corresponds to 100 ns, and 2.0 nm with rate of 0.00002 nm/ps and dt of 30 fs.
    # nsteps of 2500000 corresponds to 75  ns, and 1.5 nm with rate of 0.00002 nm/ps and dt of 30 fs.
    # nsteps of 1428572 corresponds to 50  ns, and 1.0 nm with rate of 0.00002 nm/ps and dt of 35 fs.
    # nsteps of 1428572 corresponds to 42.9  ns, and 0.86 nm with rate of 0.00002 nm/ps and dt of 30 fs.
  else
    gmx mdrun -deffnm push -v -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet -nsteps 714286
    # nsteps of 714286 corresponds to 25 ns, and 0.5 nm with rate of 0.00002 nm/ps and dt of 35 fs.
    # nsteps of 714286 corresponds to 21.4 ns, and 0.43 nm with rate of 0.00002 nm/ps and dt of 30 fs.
  fi

  ## check position vs time and force vs time
  #xmgrace pull_pullx.xvg
  #xmgrace pull_pullf.xvg

  ## extract frames from steered MD trajectory
  cp $extract_frames extract_frames_push.py
  sed -i -e "s/step_size = 5  /step_size = ${step_size}  /g" extract_frames_push.py # change step size
  sed -i -e "s/file_in = 'pull_pullx.xvg'/file_in = 'push_pullx.xvg'/g" extract_frames_push.py
  sed -i -e "s/file_out = 'extract_frames.ndx'/file_out = 'extract_frames_push.ndx'/g" extract_frames_push.py  
  sed -i -e "s/min_dist = dist\[0\]/min_dist = dist[-1]/g" extract_frames_push.py
  sed -i -e "s/max_dist = 7.0/max_dist = dist[0]/g" extract_frames_push.py
  sed -i -e "s/np.where(dist>d)/np.where(dist<d)/g" extract_frames_push.py

  $py3 extract_frames_push.py # generates extract_frames.ndx
  echo 0 | gmx trjconv -f push.xtc -s push.tpr -n  -fr extract_frames_push.ndx -o extract_frames_push.xtc -quiet

  ## check number of extracted frames
  #gmx check -f extract_frames.xtc -quiet
fi


# copy pull files from other directory
if [ $COPYPULL -eq 1 ]
then
  cp $pull_folder/pull.* .
  cp $pull_folder/pull_*.xvg .
      
  ## extract frames from steered MD trajectory 
  cp $extract_frames .
  sed -i -e "s/step_size = 5  /step_size = ${step_size}  /g" extract_frames.py # change step size
  $py3 extract_frames.py # generates extract_frames.ndx
  echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -fr extract_frames.ndx -o extract_frames.xtc -quiet
fi   

############ RUN SIMULATION ##############################################################
       
## get number of frames (3 lines in file are not frames numbers)
Nframes_pull=-3; while read -r LINE; do (( Nframes_pull++ )); done < extract_frames.ndx
echo $Nframes_pull

Nframes_push=-3; while read -r LINE; do (( Nframes_push++ )); done < extract_frames_push.ndx
echo $Nframes_push    

Nframes=$((Nframes_push+Nframes_pull))
echo $Nframes

if [ $GROMPP -eq 1 ]
then
  ## change production time and spring constant of umbrella run
  cp $umb umbrella.mdp
  sed -i -e "s/pull_coord1_k           = 1000/pull_coord1_k           = ${spring_const}/g" umbrella.mdp
      
  if [ $NO_CONSTRAINT -eq 1 ]
  then
    sed -i -e "s/define                   = -DPOSRES_POP2/;define                   = -DPOSRES_POP2/g" umbrella.mdp
  fi

  for frame in $(seq 1 $Nframes)
  do
    
    if [ $frame -le $Nframes_push ]
    then 
	cat << EOF > frames_step.ndx
 [ frames ]

$frame

EOF
	echo 0 | gmx trjconv -f extract_frames_push.xtc -s push.tpr -fr frames_step.ndx -o step_${frame}.gro -quiet
    else
	cat << EOF > frames_step.ndx
 [ frames ]

$((frame-Nframes_push))

EOF
	echo 0 | gmx trjconv -f extract_frames.xtc -s pull.tpr -fr frames_step.ndx -o step_${frame}.gro -quiet
    fi

    ## grompp umbrella production runs
    gmx grompp -f umbrella.mdp -c step_${frame}.gro -r step_${frame}.gro -n index.ndx -p ../topol.top -o umbrella_step${frame}.tpr -maxwarn 1 -quiet
    	
    ## add lines to summary files for WHAM
    echo "umbrella_step${frame}.tpr" >> tpr_files_umbrella.dat
    echo "umbrella_step${frame}_pullf.xvg" >> pullf_files_umbrella.dat
    
  done
fi	
    
if [ $MDRUN -eq 1 ]
then
  ## mdrun
  for frame in $(seq 1 $Nframes)
  do
    gmx mdrun -deffnm umbrella_step${frame} -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet -v -nb $nb -nsteps $nsteps
  done
fi

# cp summary files - ready for truncation
cp tpr_files_umbrella.dat tpr_files_umbrella_trunc.dat
cp pullf_files_umbrella.dat pullf_files_umbrella_trunc.dat

if [ $AXON -eq 1 ]
then
  ## scp tpr files to Axon
  axon_dir=axon_umbrella_prot${folder_prot}_lip${Folderprefix}_rep${rep}_t${time}_l${step_size}_k${spring_const}_n${umb_rep}
  mkdir -p $axon_dir
  cp umbrella_step*.tpr $axon_dir
  scp -r $axon_dir axon:/home/bioc1642/
  # remove folder after scp
  rm -r $axon_dir
fi
    
############ ANALYSIS ###################################################################
if [ $ANALYSIS -eq 1 ]
  then

  ## run WHAM (weighted histogram analysis method)
  skip=$(( time*1000/5*0 )) # skip first fifth of frames (equilibration) for WHAM
  #export OMP_NUM_THREADS=$noCPUs
  gmx wham -tol $tol -it tpr_files_umbrella_trunc.dat -if pullf_files_umbrella_trunc.dat -hist -o -bsres -temp 323 -quiet -b $skip #-nBootstrap 200 # outcomment "-nBootstrap 200" to get uncertainties  
      
  # check PMF curve
  #xmgrace profile.xvg
  #xmgrace bsResults.xvg
  #xmgrace bsProfs.xvg
  #xmgrace -nxy histo.xvg
fi

############ FINISH #####################################################################

## clean up
rm \#*

## navigate back to parent directory
cd ../../..


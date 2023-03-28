echo '#!/bin/bash --login' > $1.sh
echo '#SBATCH -J '$1 >> $1.sh
out=$1
file=-%J.o
outfile=$out$file
echo '#SBATCH -o '$outfile >> $1.sh
echo '#SBATCH --ntasks=1' >> $1.sh
echo '#SBATCH --cpus-per-task=8' >> $1.sh
echo '#SBATCH -N 1' >> $1.sh
echo '#SBATCH --mem-per-cpu=5G' >> $1.sh
echo '#SBATCH -p mendoza_q' >> $1.sh
time=$2
wall=:00:00
timewall=$time$wall
echo '#SBATCH -t '$timewall >> $1.sh
echo 'export JOB='$1 >> $1.sh

echo 'export DIR=$SLURM_SUBMIT_DIR
export scratch=$SCRATCH/crys22-beta

ml -* CRYSTAL/22-beta

echo "submit directory: "
echo $SLURM_SUBMIT_DIR

rm -fr $scratch/$JOB
mkdir  -p $scratch/${JOB}

# the following line need be modified according to where your input is located
cp $DIR/INPUT/${JOB}.d12  $scratch/${JOB}/INPUT
cp $DIR/INPUT/${JOB}.sopseud $scratch/${JOB}/SOPSEUD
#cp $DIR/INPUT/${JOB}.f9 $scratch/${JOB}/fort.20
cp $DIR/INPUT/${JOB}.f9 $scratch/${JOB}/fort.9

export FILE = $DIR/INPUT/${JOB}.f20
if [ -f "$FILE"]]; then
   #Check if f20 exists
   cp $DIR/INPUT/${JOB}.f20 $scratch/$JOB/fort.20
else
  #USE f9 file as substitute
   cp $DIR/INPUT/${JOB}.f9 $scratch/$JOB/fort.20
fi

cd $scratch/$JOB

# in the following, -np parameters should be equal to those specified above.
srun Pcrystal 2>&1 >& $DIR/${JOB}.out
cp fort.20 ${DIR}/${JOB}.f20
cp fort.9 ${DIR}/${JOB}.f9

# uncomment the next 5 lines if you want to remove the scratch directory
#if [ $? -eq 0 ]
#then
#    cd ${DIR}
#    rm -rf $scratch/${JOB}
#fi' >> $1.sh

sbatch $1.sh

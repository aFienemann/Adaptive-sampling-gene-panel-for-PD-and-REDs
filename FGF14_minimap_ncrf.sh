#! /bin/bash
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=shortterm
#  Use one node:
#  Request 10 cores (hard constraint):
#SBATCH -c 30
#  Request 10GB of memory (hard constraint):
#SBATCH --mem=100GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=1-0:0:0
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=500G
#  Notify me at job start and end:
#SBATCH --mail-type=ALL
#  Send the notifications to:
#SBATCH --mail-user=
#  Find your job easier with a name:
#SBATCH --job-name=ncrf_map

#Initialize the module system:
source /etc/profile.d/modules.sh
# Allow aliases (required by some modules):
shopt -s expand_aliases
# Load your necessary modules 
module load singularity/v3.5.3 
module load python3/3.9.2
module load minimap2/v2.22
module load bcftools/1.9 
module load nanostat/1.5.0 
module load samtools/1.15 
module load ncrf/v1.01.02 

#make directories
mkdir $SCRATCH/inputdata/

#copy files into inputdata
cp -a /Path/barcode95.bam $SCRATCH/inputdata/
cp -a /Path/FGF14.fa $SCRATCH/inputdata/
ls $SCRATCH/inputdata/

#make bam into fastq
samtools bam2fq $SCRATCH/inputdata/barcode95.bam > $SCRATCH/inputdata/all.fastq
ls $SCRATCH/inputdata/

cd $SCRATCH/inputdata/
#alignment with minimap2
minimap2 -a -x map-ont $SCRATCH/inputdata/FGF14.fa all.fastq > all.sam #aligns sequences with reference.fa
	samtools view -S -b all.sam > all.bam
	samtools sort all.bam -o all.sort.bam
	samtools view -h -F 4 all.sort.bam | perl -lane '$l = 0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l > 300 or /^@/' | samtools view -bS - > filt_FGF14.sorted.bam
	
	samtools sort filt_FGF14.sorted.bam -o filt_FGF14.sorted.bam
	samtools index filt_FGF14.sorted.bam
	cp -a $SCRATCH/inputdata/filt_FGF14.sorted.bam /Path/
	cp -a $SCRATCH/inputdata/filt_FGF14.sorted.bam.bai /Path/
	samtools fastq /Path/filt_FGF14.sorted.bam > lenfilt_FGF14.fastq
	
echo "finished with minimap"

	export d=barcode95
	echo $d
	#select a motif
	MOTIF1=GAA
	#output folder
	export JOB_OUTPUT_DIR="/Path"
	CONFIG=/data/neuro_mtnanopore/NCRF_AF/$d
	IDENTIFIER=$d
	mkdir $JOB_OUTPUT_DIR/$d/
	mkdir $SCRATCH/inputdata/$d/
	ls $SCRATCH/inputdata/

cat lenfilt_FGF14.fastq | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $SCRATCH/fastq/$d.fasta
	cp -a $SCRATCH/fastq/$d.fasta $CONFIG
ls $SCRATCH/inputdata/
	# Run for Motif 1
	SAMPLE=`basename $d`
	echo "ls $d"
	ls $d
	cd $d
	
	cat $SCRATCH/fastq/$d.fasta \
	| singularity exec $NCRF_CONTAINER NCRF $MOTIF1 --stats=events --positionalevents --maxnoise=80%  --minlength=10 \
	| singularity exec $NCRF_CONTAINER ncrf_sort.py --sortby=mratio  \
	| tee $CONFIG/"${MOTIF1}_raw_${IDENTIFIER}.summary" \
	| singularity exec $NCRF_CONTAINER ncrf_consensus_filter.py \
	| singularity exec $NCRF_CONTAINER ncrf_sort.py --sortby=mratio  \
	| tee $CONFIG/"${MOTIF1}_refined_${IDENTIFIER}.summary" \
	| singularity exec $NCRF_CONTAINER ncrf_summary.py \
	> $CONFIG/"${MOTIF1}_summary_${IDENTIFIER}.summary"

echo "--->NCRF-MOTIF1=${MOTIF1}-Run concluded"
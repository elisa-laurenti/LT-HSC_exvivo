#!/bin/bash
# set the number of nodes
#SBATCH -p partition
#SBATCH -A ACCOUNT-SL2-CPU
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=user.name@cam.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --time 72:00:00
#SBATCH --job-name star
#SBATCH --output star_%J.log

#########################################################
fastqdir=fastq_dir
outdir=output_dir

module load star
STAR --version

module load samtools
samtools --version

mkdir -p ${outdir}
cd ${outdir}
echo pwd

mkdir STARout
mkdir fcounts

###########################################################
### STAR ######################
genomeDir=STAR_genome_index/star_ref_GRCh38_and_ERCC92_overhang49

filetype=.fq.gz
if [[ ${filetype} =~ '.fq.gz' ]]
then
        readFilesCommand=zcat

elif [[ ${filetype} =~ '.fastq' ]]
then
        readFilesCommand=-
elif [[ ${filetype} =~ '.bz2' ]]
then
        readFilesCommand=bzcat
fi

runThreadN=`nproc`
genomeLoad=LoadAndKeep
outSAMtype=SAM
outFilterMultimapNmax=1
outSAMunmapped=Within

### featureCounts ###############
a=references/gtf/Homo_sapiens.GRCh38.81_and_ERCC92.gtf

T=24
t=exon
g=gene_id

##########################################################
for f in ${fastqdir}/*${filetype}
do
        fn=$(basename $f)

        if [[ ${fn} =~ 'lostreads' ]]
        then
                continue
        fi

        if [[ ${fn} =~ 'r_2.fq.gz' ]]
        then
                continue
        fi


        cell_id=${fn%.r_*.fq.gz}

        # echo "$f"
        echo "$fn"
        echo $cell_id

        read_1=${fastqdir}/${cell_id}.r_1.fq.gz
        read_2=${fastqdir}/${cell_id}.r_2.fq.gz

        if [[ -f $read_1 && -f $read_2 ]]; then
                echo 'both exist'
        else
                echo "read missing == "$read_1
                echo "read missing == "$read_2

                continue
        fi

        echo "Running STAR for ${cell_id} ..."
        STAR \
                --runThreadN ${runThreadN} \
                --genomeLoad ${genomeLoad} \
                --genomeDir ${genomeDir} \
                --outSAMtype ${outSAMtype} \
                --outFilterMultimapNmax ${outFilterMultimapNmax} \
                --outFileNamePrefix ${outdir}/STARout/${cell_id}_ \
                --outSAMunmapped ${outSAMunmapped} \
                --readFilesCommand ${readFilesCommand} \
                --readFilesIn $read_1 $read_2

        echo "Converting SAM files for ${cell_id} ..."
        samtools view -S -b -@24 ${outdir}/STARout/${cell_id}_Aligned.out.sam > ${outdir}/STARout/${cell_id}_Aligned.out.bam
        samtools sort -@ 24 ${outdir}/STARout/${cell_id}_Aligned.out.bam > ${outdir}/STARout/${cell_id}_Aligned.out.sort.bam
        samtools index ${outdir}/STARout/${cell_id}_Aligned.out.sort.bam

        echo "Cleaning up ..."
        rm ${outdir}/STARout/${cell_id}_Aligned.out.bam
        rm ${outdir}/STARout/${cell_id}_Aligned.out.sam
done

echo "Cleaning up STAR left overs ..."
STAR \
        --genomeLoad Remove \
        --genomeDir ${genomeDir}

rm Log.out
rm Log.progree.out
rm Aligned.out.sam
rm -r _STARtmp

##########################################################
echo "Running featureCounts for ${cell_id} ..."

bams=${outdir}/STARout/*_Aligned.out.sort.bam
path/to/featureCounts \
       -T ${T} \
       -t ${t} \
       -g ${g} \
       -a ${a} \
       -o ${outdir}/fcounts/fcounts.txt \
       ${bams}



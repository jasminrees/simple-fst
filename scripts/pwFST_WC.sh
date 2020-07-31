#$ -S /bin/bash #defines the used shell
#$ -cwd
#$ -q all.q
# load personal profile
#$ -V
#$ -m e
#$ -M joshua_schmidt@eva.mpg.de
#$ -l h_vmem=3G,virtual_free=3G
#$ -e /mnt/scratch/josh/sge_josh_logs/1000Genomes/fst
#$ -o /mnt/scratch/josh/sge_josh_logs/1000Genomes/fst


# Job Array with jobs corresponding to each chromosome 

#ARGV
pop1=${1}
pop2=${2}
chr=${SGE_TASK_ID}

#PATHs
scripts="/mnt/scratch/josh/chimpanzee_pop/scripts"
popIndvs="/mnt/scratch/josh/FST_test"
vcfs="/mnt/scratch/josh/1000Genomes/vcf_allpop1586"
outpath="/mnt/scratch/josh/1000Genomes/chimp_sample_size_comparison/chr_${chr}"

mkdir -p ${outpath}
log=${outpath}/${pop1}_${pop2}_fst_log.txt
echo $pop1 > ${log}
echo $pop2 >> ${log}
echo "perl ${scripts}/pairwise_FST_from_VCF_min.pl ${vcfs}/allpop1586_chr${chr}.vcf.gz ${popIndvs}/${pop1}.txt ${popIndvs}/${pop2}.txt ${chr} ${outpath}" >>${log}
perl ${scripts}/pairwise_FST_from_VCF_min.pl ${vcfs}/allpop1586_chr${chr}.vcf.gz ${popIndvs}/${pop1}.txt ${popIndvs}/${pop2}.txt ${chr} ${outpath} 2>> ${log}
exit 0

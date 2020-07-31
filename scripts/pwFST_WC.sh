#$ -S /bin/bash #defines the used shell
#$ -cwd
#$ -q all.q
# load personal profile
#$ -V
#$ -m e
#$ -l h_vmem=3G,virtual_free=3G
#$ -e /dev/null
#$ -o /dev/null


# Job Array with jobs corresponding to each chromosome 

#ARGV
pop1=${1}
pop2=${2}
chr=${SGE_TASK_ID}

#PATHs
scripts="/path/to/scripts"
popIDs="/path/to/ids"
vcfs="path/to/vcfs"
outpath="/path/to/output/chr_${chr}"

mkdir -p ${outpath}
log=${outpath}/${pop1}_${pop2}_fst_log.txt
echo $pop1 > ${log}
echo $pop2 >> ${log}
echo "perl ${scripts}/pairwise_FST_WC.pl ${vcfs}/allpop1586_chr${chr}.vcf.gz ${popIDs}/${pop1}.txt ${popIDs}/${pop2}.txt ${chr} ${outpath}" >>${log}
perl ${scripts}/pairwise_FST_WC.pl ${vcfs}/allpop1586_chr${chr}.vcf.gz ${popIDs}/${pop1}.txt ${popIDs}/${pop2}.txt ${chr} ${outpath} 2>> ${log}
exit 0

#PBS -l walltime=2036:00:00
#PBS -l nodes=1:ppn=10
#PBS -j oe
#PBS -o ${out}.results.txt
#PBS -V 
cd $PBS_O_WORKDIR

pts='-p 8 --bowtie2 --bowtie2-sensitivity-level very_sensitive --phred33-quals'
ropts='--no-bam-output --estimate-rspd'
rsemref='/public/software/genomics/genome/mm10/mm10v79_rsem_index/mm10_enst_v79'

rsem-calculate-expression $pts $ropts --paired-end ${p1} ${p2} $rsemref ${out}


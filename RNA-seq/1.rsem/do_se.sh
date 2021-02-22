#PBS -l walltime=336:00:00
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -o ${out}.results.txt
#PBS -V 
cd $PBS_O_WORKDIR

pts='-p 6 --bowtie2 --bowtie2-sensitivity-level very_sensitive --phred33-quals'
ropts='--no-bam-output --estimate-rspd'
rsemref='/public/software/genomics/genome/mm10/mm10v79_rsem_index/mm10_enst_v79'

rsem-calculate-expression $pts $ropts ${inp} $rsemref ${out}

 



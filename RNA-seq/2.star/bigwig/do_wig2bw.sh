#PBS -l walltime=1336:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -q default
#PBS -o ${out}.txt
#PBS -V 
cd $PBS_O_WORKDIR

wigToBigWig $inp mm10.chrom.sizes ${out}.bw


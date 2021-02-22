#PBS -l walltime=1336:00:00
#PBS -l nodes=1:ppn=20,mem=60gb
#PBS -j oe
#PBS -q cv2
#PBS -o s97_20191223-5-E13.0.log
#PBS -V
cd $PBS_O_WORKDIR
ulimit -n 2000

whitelist="--soloCBwhitelist /public/home/tools/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt --outSAMattributes NH HI AS nM CR CY UR UY"
mods="--soloType Droplet --soloFeatures Gene SJ --soloBarcodeReadLength 0 "
teopts=" --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --twopassMode Basic"
opts="--soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --runRNGseed 42 --runThreadN 20 --readFilesCommand zcat"
ref="--genomeDir /public/home/omicsdata/genome/mm10/star_index/star2.7.0.f/mm10"
path="../fqs"

/public/software/genomics/unstable/STAR-2.7/bin/Linux_x86_64_static/STAR $opts $teopts $mods $whitelist --outFileNamePrefix s97_20191223-5-E13.0 $ref --readFilesIn $path/20191223_5_1_BKDL192553839-1a-AK919_R2.fq.gz,$path/20191223_5_2_BKDL192553839-1a-AK920_R2.fq.gz,$path/20191223_5_3_BKDL192553839-1a-AK921_R2.fq.gz,$path/20191223_5_4_BKDL192553839-1a-AK922_R2.fq.gz $path/20191223_5_1_BKDL192553839-1a-AK919_R1.fq.gz,$path/20191223_5_2_BKDL192553839-1a-AK920_R1.fq.gz,$path/20191223_5_3_BKDL192553839-1a-AK921_R1.fq.gz,$path/20191223_5_4_BKDL192553839-1a-AK922_R1.fq.gz
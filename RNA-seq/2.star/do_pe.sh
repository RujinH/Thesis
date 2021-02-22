#PBS -l walltime=1336:00:00
#PBS -l nodes=1:ppn=15
#PBS -j oe
#PBS -q default
#PBS -o ${out}.results.txt
#PBS -V 
cd $PBS_O_WORKDIR

opt='--readFilesCommand zcat --outMultimapperOrder Random --runRNGseed 777 --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic --sjdbOverhang 100 --outWigType wiggle --outWigNorm RPM --outWigStrand Unstranded'
encode='--outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000'
ref='--genomeDir /public/home/genome/mm10/star_index/mm10'
splice=' --sjdbGTFfile /public/home/genome/mm10/star_index/mm10.ensemblv67.nopsuedo.gtf'

STAR $opt $encode $ref $splice --runThreadN 15 --readFilesIn ${p1} ${p2} --outFileNamePrefix ${out} --outTmpDir ${out}

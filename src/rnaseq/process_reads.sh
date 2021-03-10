# Abort script on any error
set -e
set -o pipefail

mkdir -p FASTQC
mkdir -p STAR
mkdir -p STAR/BAMs
mkdir -p STAR/BAMs/merged
mkdir -p STAR/BAMs/merged/sorted
mkdir -p STAR/Log.final.out
mkdir -p STAR/Logs

# collecting file names from unmapped bams

samples="$*"
echo $samples
if [ "$samples" == "" ]; then
	for f in $(cd unmapped_bams; ls *.bam); do
		samples="$samples ${f%%.*}"
	done
	echo "Will process:"
	echo $samples | xargs -n 1 echo -e "\t"
fi

# convert unmapped bams to fastqs and map reads (STAR version STAR 2.6.0a)

for sample in $samples; do
	echo unmapped_bams $sample

	if [ ! -e fastq/$sample.fastq.gz ]; then
		echo bam2fq unmapped_bams/$sample.bam
		bedtools bamtofastq -i unmapped_bams/$sample.bam -fq fastq/$sample.fastq
		echo gzip fastq/$sample.fastq
		gzip fastq/$sample.fastq
	else
		echo fastq/$sample.fastq.gz already exists
	fi



	if [ ! -e STAR/BAMs/$sample.Aligned.sortedByCoord.out.bam ]; then
		echo map fastq/$sample.fastq.gz

		STAR --runThreadN 6 --genomeDir STAR_index_mm10_Allelome --readFilesIn fastq/$sample.fastq.gz --outFilterType BySJout\
		--outFilterMultimapNmax 1 --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate\
		--outFileNamePrefix STAR/$sample. --readFilesCommand zcat
		rm -r STAR/$sample._STARtmp
		
		mv STAR/$sample.Log.final.out STAR/Log.final.out
		mv STAR/$sample.Aligned.sortedByCoord.out.bam STAR/BAMs
		mv STAR/$sample.Log.out STAR/Logs
		mv STAR/$sample.Log.progress.out STAR/Logs
		mv STAR/$sample.SJ.out.tab STAR/Logs

	else
		echo STAR/BAMs/$sample.Aligned.sortedByCoord.out.bam already exists
	fi
		
done

# new sample names to merge lanes of from sequencing

samples="$*"
if [ "$samples" == "" ]; then
	for f in $(cd STAR/BAMs; ls *.bam); do
		samples="$samples ${f%%_L00*}"
	done
	echo "Will process:"
	echo $samples | xargs -n 1 echo -e "\t"
fi

# merging bam files from different lanes (samtools version 1.9-92)

for sample in $samples; do
	echo STAR/BAMs $sample

	if [ ! -e STAR/BAMs/merged/$sample.bam ]; then
		echo Creating STAR/BAMs/merged/$sample.bam
		samtools merge STAR/BAMs/merged/$sample.bam STAR/BAMs/$sample*
	else
		echo STAR/BAMs/merged/$sample.bam already exists
	fi
done

# sorting and indexing merged bam files

samples="$*"
if [ "$samples" == "" ]; then
	for f in $(cd STAR/BAMs/merged; ls *.bam); do
		samples="$samples ${f%%.*}"
	done
	echo "Will process:"
	echo $samples | xargs -n 1 echo -e "\t"
fi

for sample in $samples; do
	if [ ! -e STAR/BAMs/merged/sorted/$sample.bam ]; then
		echo Creating STAR/BAMs/merged/sorted/$sample.bam
	        samtools sort STAR/BAMs/merged/$sample.bam -o STAR/BAMs/merged/sorted/$sample.bam
	else
		echo STAR/BAMs/merged/sorted/$sample.bam already exists
	fi


	if [ ! -e STAR/BAMs/merged/sorted/$sample.bam.bai ]; then
		echo Creating STAR/BAMs/merged/sorted/$sample.bam.bai
		samtools index STAR/BAMs/merged/sorted/$sample.bam
	else
		echo STAR/BAMs/merged/sorted/$sample.bam.bai already exists
	fi

	
done

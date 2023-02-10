#!/usr/bin/env nextflow

/*
*****************************
 * Pipeline - smallRNA-meth *
 *                          *
 * Thomas R Burkard         *
 * IMP/IMBA Bioinformatics  *
 ****************************
*/


log.info """\
         =============================
         Pipeline - ATAC-seq process
         =============================

         outdir: ${params.outdir}
         reads: ${params.reads}
         splitted_readPair: ${params.splitted_readPair}
         genome: ${params.genome}
         adapter: ${params.adapter}
         tss: ${params.tss}

         """
         .stripIndent()

/*
 * Input parameters validation
 */

if (params.genome)
{
  genome_file = file(params.genome)
  //if( !genome_file.exists() ) exit 1, "Genome file doesn't exist: ${genome_file}"

} else {
  exit 1, "Missing genome file"
}

if (params.tss)
{
  tss_file = file(params.tss)
  if( !tss_file.exists() ) exit 1, "tss file doesn't exist: ${tss_file}"
} else {
  tss_file = file("NA")
}

/*
 * Validate input files
 */
multiqc_config = file(params.multiqc_config)

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for read files
 */
if(!params.splitted_readPair)
{
    Channel
        .fromPath(params.reads)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .map { file -> tuple(file.baseName, file) }
        .set { read_files }

 } else {
    Channel
        //.from(params.reads)
        .fromFilePairs( params.reads )
        .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        //.view()
        .set { read_files}
        //.map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
        //.into { raw_reads_fastqc; raw_reads_trimgalore }

 }


/*
 * bam to fastq and pair splitting
 */
process splitfastq {
    tag "Channel: ${name}"

    publishDir "${params.outdir}/FASTQs", mode: 'copy', pattern: '*.fastq'

    input:
        set val(name), file(bam) from read_files

    output:
        set name, file("${name}_R*.fastq") into splited_fastq
        set name, file("${name}_R*.fastq") into fastq_4fastqc
        set name, file("cntTotal.txt") into cnt_total

    script:
    if( ! params.splitted_readPair)
      """
        ml load bedtools/2.25.0-foss-2018b
        ml load samtools/1.10-foss-2018b

        samtools view -c ${bam} > cntTotal.txt
        bamToFastq -i ${bam} -fq ${name}.fastq
        bash ${baseDir}/scripts/deinterleave_fastq.sh < ${name}.fastq ${name}_R1.fastq ${name}_R2.fastq

      """
    else
      """
        echo 'splitted already'
        mv *R1*.fastq.gz ${name}_R1.fastq.gz
        mv *R2*.fastq.gz ${name}_R2.fastq.gz
        gunzip -f ${name}_R1.fastq.gz
        gunzip -f ${name}_R2.fastq.gz
        cat ${name}_R1.fastq ${name}_R2.fastq > merged.fastq
        cat merged.fastq | paste - - - - | wc -l > cntTotal.txt
        rm merged.fastq

      """
}


/*
 *  FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from fastq_4fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    ml load fastqc/0.11.8-java-1.8
    fastqc -q $reads

    """
}

/*
 * cut adapter
 */
process  cutadapt {
    tag "Channel: ${name}"

    publishDir "${params.outdir}/cutadapt", mode: 'copy'

    input:
        set val(name), file(fastq) from splited_fastq

    output:
        set name, file("${name}_R*.trim.fastq") into fastq_trimmed
        set name, file("cutadapt.${name}.err") into stat_cutadapt
        set name, file("cnt_cutadapt_*") into cnt_cutadapt

    script:
    """
        module load python/2.7.15-gcccore-7.3.0-bare
        module load cutadapt/1.18-foss-2018b-python-2.7.15

        cutadapt --minimum-length ${params.minLength} --overlap ${params.overlapLength} -a ${params.adapter} -A ${params.adapter} \
        -o ${name}_R1.trim.fastq -p ${name}_R2.trim.fastq ${name}_R1.fastq ${name}_R2.fastq > cutadapt.${name}.err

        cat ${name}_R1.trim.fastq | paste - - - - | wc -l > cnt_cutadapt_R1.txt
        cat ${name}_R2.trim.fastq | paste - - - - | wc -l > cnt_cutadapt_R2.txt

    """
}

/*
 * Align with bowite2
 */
process align {

    tag "Channel: ${name}"
    publishDir "${params.outdir}/aligned_bams", mode: 'copy'

    input:
    set name, file(fastq) from fastq_trimmed

    output:
      set name, file("${name}.bam") into aligned_bam
      file "${name}.bam"

    script:
    """
    ml load bowtie2/2.3.5.1-foss-2018b
    ml load samtools/1.10-foss-2018b

    bowtie2 -q -p 32 --no-mixed -X 2000 --dovetail --no-discordant -x ${params.genome} -1 ${name}_R1.trim.fastq -2 ${name}_R2.trim.fastq | samtools view -bSu - > ${name}.unsorted.bam
    samtools sort -o ${name}.bam ${name}.unsorted.bam
    samtools index -c -m 14 ${name}.bam
    rm ${name}.unsorted.bam

    """
}


/*
 * filter algined bam, duplicated removal
 */
process filter_rmdup_bam {

    tag "Channel: ${name}"
    publishDir "${params.outdir}/filtered_bams", mode: 'copy'

    input:
    set name, file(bam) from aligned_bam

    output:
      set name, file("${name}_uniq.rmdup.bam") into filterd_bam, bams_peakcalling

      set name, file ("${name}_cnt_trimmed.txt") into trimmed_stat
      set name, file ("${name}_cnt_mapped.txt") into mapped_stat
      set name, file ("${name}_cnt_rmChrM.txt") into rmChrM_stat
      set name, file ("${name}_cnt_uniq.txt") into unique_stat
      set name, file ("${name}_cnt_uniq.rmdup.txt") into usable_stat

      file "${name}_uniq.rmdup.bam"
      file "${name}_uniq.rmdup.bam.csi"
      file '*.pdf'

    script:
    """
      ml load samtools/1.10-foss-2018b
      ml load picard/2.20.6--0-biocontainers

      samtools sort --threads 16 -o ${name}.total.bam ${bam}
      samtools index -c -m 14 -@ 16 ${name}.total.bam

      samtools view --threads 16 -h ${name}.total.bam | grep -v chrM | samtools sort -O bam -o ${name}.rmChrM.bam

      samtools view --threads 16 -h -q 30 ${name}.rmChrM.bam > ${name}.rmChrm.unsorted.bam

      samtools sort --threads 16 -o ${name}_uniq.bam ${name}.rmChrm.unsorted.bam
      samtools index -c -m 14 -@ 16 ${name}_uniq.bam
      rm ${name}.rmChrm.unsorted.bam

      picard MarkDuplicates INPUT=${name}_uniq.bam OUTPUT=${name}_uniq.rmdup.unsorted.bam METRICS_FILE=${name}_picard.rmDup.txt \
      ASSUME_SORTED=true REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.2

      samtools sort --threads 16 -o ${name}_uniq.rmdup.bam ${name}_uniq.rmdup.unsorted.bam
      samtools index -c -m 14 -@ 16 ${name}_uniq.rmdup.bam
      rm ${name}_uniq.rmdup.unsorted.bam

      picard CollectInsertSizeMetrics HISTOGRAM_FILE=${name}_fragsize.pdf OUTPUT=${name}_frag_size.txt \\
      METRIC_ACCUMULATION_LEVEL=ALL_READS INCLUDE_DUPLICATES=false INPUT=${name}_uniq.rmdup.bam

      samtools view --threads 16 -c ${name}.total.bam > ${name}_cnt_trimmed.txt
      samtools view --threads 16 -c -F 4 ${name}.total.bam > ${name}_cnt_mapped.txt
      samtools view --threads 16 -c ${name}.rmChrM.bam > ${name}_cnt_rmChrM.txt
      samtools view --threads 16 -c ${name}_uniq.bam > ${name}_cnt_uniq.txt
      samtools view --threads 16 -c ${name}_uniq.rmdup.bam > ${name}_cnt_uniq.rmdup.txt

    """
}


/*
 * make bigwig with deeptools
 */
process bigwig_deeptools {

    tag "Channel: ${name}"
    publishDir "${params.outdir}/bigwigs", mode: 'copy'

    input:
    set name, file(bam) from filterd_bam

    output:
      file "*.bw"

    shell:
    """
    ml load samtools/1.10-foss-2018b

    samtools sort --threads 16 -o !{name}.sorted.bam !{bam}
    mv !{name}.sorted.bam !{name}.bam
    samtools index -c -m 14  -@ 16 !{name}.bam

    singularity exec --no-home --home /tmp /groups/tanaka/People/current/jiwang/local/deeptools_master.sif bamCoverage \
    -b !{name}.bam -o !{name}.bw --outFileFormat=bigwig --normalizeUsing CPM -p 16 --binSize 20

    """
}


process macs2_peak {

    tag "Channel: ${name}"
    publishDir "${params.outdir}/peaks_macs2", mode: 'copy'

    input:
    set name, file(bam) from bams_peakcalling

    output:
      file "*.{bed,xls,narrowPeak}"

    shell:
    """
    ml load macs2/2.2.5-foss-2018b-python-3.6.6
    macs2 callpeak -t $bam -n ${name} -f BAMPE -g 30000000000 -p 0.001

    """
}


/*
 * make statTable
 */
 cnt_total.concat(trimmed_stat, mapped_stat, rmChrM_stat, unique_stat, usable_stat)
     .groupTuple()
     .map{ stat1, stat2 -> [stat1, stat2[0], stat2[1], stat2[2], stat2[3], stat2[4], stat2[5]] }
     .set{ cntStat_files }

process statTable {

    tag "Channel: ${name}"

    publishDir "${params.outdir}/result/countStat", mode: 'copy'
    input:
        set name, file(cnt_total), file(trimmed), file(mapped), file(rmChrM_stat), file(unique_stat), file(usable_stat) from cntStat_files

    output:
        file "${name}.countStat.txt" into cnt_stat

    script:
    """
    echo -e "Name\tTotal\trimmed\taligned\trmChrM\tunique\tunique.rmdup" > ${name}.countStat.txt
    TOTAL=`cat ${cnt_total}`
    TRIMMED=`cat ${name}_cnt_trimmed.txt`
    MAPPED=`cat ${name}_cnt_mapped.txt`
    rmChrM=`cat ${name}_cnt_rmChrM.txt`
    UNIQUE=`cat ${name}_cnt_uniq.txt`
    USABLE=`cat ${name}_cnt_uniq.rmdup.txt`

    echo -e "${name}\t\$TOTAL\t\$TRIMMED\t\$MAPPED\t\$rmChrM\t\$UNIQUE\t\$USABLE" >> ${name}.countStat.txt

    """
}

/*
 *  stat table of all samples
 */
process countTable {

    publishDir "${params.outdir}/result", mode: 'copy'

    input:
	      file "countStat/*" from cnt_stat.collect()

    output:
	      file "countStatTable.txt"

    script:
    """
    ml load r/3.6.2-foss-2018b
    cp $baseDir/scripts/countTable.R ./countTable.R
    Rscript countTable.R

    """
}

/*
 * MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/result/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    ml load multiqc/1.9-foss-2018b-python-3.6.6
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m fastqc

    """
}


workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Execution status      : ${ workflow.success ? 'succeeded' : 'failed' }"
    println "Duration : ${workflow.duration}"
}

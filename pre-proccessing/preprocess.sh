gatk=$(grep '^gatk:' $CONFIG | awk '{print $2}')
picard=$(grep '^picard:' $CONFIG | awk '{print $2}')
bwa=$(grep '^bwa:' $CONFIG | awk '{print $2}')
reference=$(grep '^reference:' $CONFIG | awk '{print $2}')
output_dir=$(grep '^output_dir:' $CONFIG | awk '{print $2}')
threads=$(grep '^threads:' $CONFIG | awk '{print $2}')
platform=$(grep '^platform:' $CONFIG | awk '{print $2}')
library=$(grep '^library:' $CONFIG | awk '{print $2}')

known1=$(grep -A2 '^known_sites:' $CONFIG | sed -n '2p' | awk '{print $2}')
known2=$(grep -A2 '^known_sites:' $CONFIG | sed -n '3p' | awk '{print $2}')

OUTDIR="${output_dir}/${SAMPLE_ID}"
mkdir -p "$OUTDIR"

echo "[INFO] Processing $SAMPLE_ID"
echo "[INFO] Output: $OUTDIR"

java -jar $picard FastqToSam \
    F1=$FASTQ_R1 F2=$FASTQ_R2 \
    O=$OUTDIR/${SAMPLE_ID}.unmapped.bam \
    SM=$SAMPLE_ID PL=$platform LB=$library

java -jar $picard ValidateSamFile \
    I=$OUTDIR/${SAMPLE_ID}.unmapped.bam \
    MODE=SUMMARY

java -jar $picard SamToFastq \
    I=$OUTDIR/${SAMPLE_ID}.unmapped.bam \
    FASTQ=$OUTDIR/${SAMPLE_ID}_interleaved.fq

$bwa mem -t $threads $reference \
    $OUTDIR/${SAMPLE_ID}_interleaved.fq > $OUTDIR/${SAMPLE_ID}.aligned.sam

$gatk MergeBamAlignment \
    --UNMAPPED_BAM $OUTDIR/${SAMPLE_ID}.unmapped.bam \
    --ALIGNED_BAM $OUTDIR/${SAMPLE_ID}.aligned.sam \
    --REFERENCE_SEQUENCE $reference \
    --OUTPUT $OUTDIR/${SAMPLE_ID}.merged.bam

$gatk MarkDuplicates \
    -I $OUTDIR/${SAMPLE_ID}.merged.bam \
    -O $OUTDIR/${SAMPLE_ID}.dedup.bam \
    -M $OUTDIR/${SAMPLE_ID}_dup_metrics.txt

$gatk SortSam \
    -I $OUTDIR/${SAMPLE_ID}.dedup.bam \
    -O $OUTDIR/${SAMPLE_ID}.sorted.bam \
    --SORT_ORDER coordinate

$gatk SetNmMdAndUqTags \
    -I $OUTDIR/${SAMPLE_ID}.sorted.bam \
    -O $OUTDIR/${SAMPLE_ID}.NMUQ.bam \
    -R $reference

$gatk BaseRecalibrator \
    -I $OUTDIR/${SAMPLE_ID}.NMUQ.bam \
    -R $reference \
    --known-sites $known1 \
    --known-sites $known2 \
    -O $OUTDIR/${SAMPLE_ID}_recal.table

$gatk ApplyBQSR \
    -I $OUTDIR/${SAMPLE_ID}.NMUQ.bam \
    -R $reference \
    --bqsr-recal-file $OUTDIR/${SAMPLE_ID}_recal.table \
    -O $OUTDIR/${SAMPLE_ID}.final.bam

$gatk BuildBamIndex \
    -I $OUTDIR/${SAMPLE_ID}.final.bam

echo "[DONE] $SAMPLE_ID"

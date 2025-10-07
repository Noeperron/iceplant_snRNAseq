# ------------------------------
# Load modules
# ------------------------------
module load fastqc
module load multiqc/1.28

# ------------------------------
# Define paths and create output directories
# ------------------------------
FASTQ_DIR=/path/to/raw_fastq
OUTDIR=/path/to/output

QC_DIR=$OUTDIR/qc
TRIM_DIR=$OUTDIR/trimmed

mkdir -p "$QC_DIR" "$TRIM_DIR"
cd "$OUTDIR"

# ------------------------------
# Run FastQC on raw reads
# ------------------------------
for SAMPLE in AAA AMA MAA MMA NAA NMA
do
    echo "Running FastQC on raw reads for $SAMPLE..."
    fastqc -t $SLURM_CPUS_PER_TASK $FASTQ_DIR/${SAMPLE}.fastq.gz -o $QC_DIR
done

# Combine into one summary report
multiqc $QC_DIR -o $QC_DIR --filename initial_multiqc_report.html

# ------------------------------
# Load modules
# ------------------------------
module load fastqc
module load multiqc/1.28
module load trimmomatic/0.39

# ------------------------------
# Define paths and create output directories
# ------------------------------
FASTQ_DIR=/path/to/raw_fastq
OUTDIR=/path/to/output

QC_DIR=$OUTDIR/qc
TRIM_DIR=$OUTDIR/trimmed

mkdir -p "$QC_DIR" "$TRIM_DIR"
cd "$OUTDIR"

# -------------------------------
# Run Trimmomatic and FastQC to ensure adapters removed
# -------------------------------
for SAMPLE in AAA AMA MAA MMA NAA NMA
do
    echo "Trimming adapters and low-quality bases for $SAMPLE..."
    trimmomatic SE -threads $SLURM_CPUS_PER_TASK -phred33 \
        $FASTQ_DIR/${SAMPLE}.fastq.gz \
        $TRIM_DIR/${SAMPLE}_trimmed.fastq.gz \
        HEADCROP:12 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:30

    echo "Running FastQC on trimmed reads for $SAMPLE..."
    fastqc -t $SLURM_CPUS_PER_TASK $TRIM_DIR/${SAMPLE}_trimmed.fastq.gz -o $QC_DIR
done

# Combine into one summary report
multiqc $QC_DIR -o $QC_DIR --filename trimmed_multiqc_report.html

echo "02_trimming_and_qc.sh ran successfully!"

# ------------------------------
# Load modules
# ------------------------------
module purge
module load bowtie2/2.5.4
module load samtools/1.20

# ------------------------------
# Define paths and create output directories
# ------------------------------
AT_REF=/path/to/Arabidopsis_index/TAIR10_nuclear   # Bowtie2 index prefix
AT_GFF=/path/to/Arabidopsis_annotation/Athaliana_167_TAIR10.gene.gff3.gz

MC_REF=/path/to/IcePlant_Index/IcePlant            # Bowtie2 index prefix (Mesembryanthemum crystallinum)
MC_GFF=/path/to/IcePlantGenome/ICE_V1_polished_MAKER2_AllModels.gff

FASTQ_DIR=/path/to/raw_fastq
OUTDIR=/path/to/output

QC_DIR=$OUTDIR/qc
TRIM_DIR=$OUTDIR/trimmed
ALIGN_DIR=$OUTDIR/bam

mkdir -p "$ALIGN_DIR"

cd "$OUTDIR"

# ------------------------------
# Sample groups
# ------------------------------
MC_SAMPLES=("MMA" "AMA" "NMA")   # Ice Plant
AT_SAMPLES=("MAA" "AAA" "NAA")   # Arabidopsis


# ------------------------------
# Loop 1: Ice Plant samples
# ------------------------------
for SAMPLE in "${MC_SAMPLES[@]}"; do
    echo "$(date) - Aligning $SAMPLE to Ice Plant reference..."

    # align and sort
    bowtie2 -x "$MC_REF" \
        -U "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
        -p 24 2> "$QC_DIR/${SAMPLE}.bowtie2.log" | \
        samtools view -@ 24 -b -F 4 - | \
        samtools sort -@ 24 -o "$ALIGN_DIR/${SAMPLE}.sorted.bam"

    samtools index -@ 24 "$ALIGN_DIR/${SAMPLE}.sorted.bam"

    # filter to mapq>30 for peak calling
    samtools view -b -q 30 "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$ALIGN_DIR/${SAMPLE}.mapq30.bam"
    samtools index "$ALIGN_DIR/${SAMPLE}.mapq30.bam"

    # qc the alignment
    samtools flagstat "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.flagstat.txt"
    samtools stats "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.stats.txt"
    samtools idxstats "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.idxstats.txt"
    samtools coverage "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.coverage.txt"

    echo "$(date) - Finished $SAMPLE"
done

# ------------------------------
# Loop 2: Arabidopsis samples
# ------------------------------
for SAMPLE in "${AT_SAMPLES[@]}"; do
    echo "$(date) - Aligning $SAMPLE to Arabidopsis reference..."

    # align and sort
    bowtie2 -x "$AT_REF" \
        -U "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
        -p 24 2> "$QC_DIR/${SAMPLE}.bowtie2.log" | \
        samtools view -@ 24 -b -F 4 - | \
        samtools sort -@ 24 -o "$ALIGN_DIR/${SAMPLE}.sorted.bam"

    samtools index -@ 24 "$ALIGN_DIR/${SAMPLE}.sorted.bam"

    # filter to mapq value greater than or equal to 30
    samtools view -b -q 30 "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$ALIGN_DIR/${SAMPLE}.mapq30.bam"
    samtools index "$ALIGN_DIR/${SAMPLE}.mapq30.bam"

    # qc stats of the alignment
    samtools flagstat "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.flagstat.txt"
    samtools stats "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.stats.txt"
    samtools idxstats "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.idxstats.txt"
    samtools coverage "$ALIGN_DIR/${SAMPLE}.sorted.bam" > "$QC_DIR/${SAMPLE}.coverage.txt"

    echo "$(date) - Finished $SAMPLE"
done

echo "All alignments completed successfully!"

# ------------------------------
# Load modules
# ------------------------------
module purge
module load bedtools

# ------------------------------
# Define directories
# ------------------------------
BAM_DIR=/path/to/bam

# ------------------------------
# Convert only mapq30 BAM files to BED
# ------------------------------
for BAM in $BAM_DIR/*mapq30.bam; do
    BASENAME=$(basename $BAM .bam)
    BED_FILE=$BAM_DIR/${BASENAME}.bed
    echo "$(date) - Converting $BAM to $BED_FILE"
    bedtools bamtobed -i $BAM > $BED_FILE
done

echo "$(date) - All mapq30 BAM files converted to BED!"

# ------------------------------
# Load modules or environment
# ------------------------------
module purge
module load mamba
mamba activate GEM_env

# ------------------------------
# Define paths
# ------------------------------
ALIGN_DIR=/path/to/bam
OUTDIR=/path/to/peaks
mkdir -p $OUTDIR

# Genome and annotation files
GENOMES=(
    "/path/to/IcePlantAssembly_Polished.fasta"    # Ice Plant
    "/path/to/TAIR10_nuclear.fasta"              # Arabidopsis
)

GFFS=(
    "/path/to/ICE_V1_polished_MAKER2_AllModels.gff"  # Ice Plant
    "/path/to/TAIR10.gff3"                            # Arabidopsis
)

# Sample groups and negative controls
SAMPLES=(
    "MMA AMA"  # Ice Plant
    "MAA AAA"          # Arabidopsis
)

CONTROLS=(
    "NMA"  # Ice Plant negative control
    "NAA"  # Arabidopsis negative control
)

SPECIES=("IcePlant" "Arabidopsis")

# ------------------------------
# Loop over species
# ------------------------------
for i in ${!SPECIES[@]}; do
    SPEC=${SPECIES[$i]}
    GENOME=${GENOMES[$i]}
    GFF=${GFFS[$i]}
    CTRL=${CONTROLS[$i]}

    for TF in ${SAMPLES[$i]}; do
        echo "$(date) - Running GEM for $TF ($SPEC)..."
        mkdir -p $OUTDIR/${SPEC}/${TF}

        gem \
        --d /path/to/GEM/Read_Distribution_default.txt \
        --expt $ALIGN_DIR/${TF}.mapq30.bed \
        --ctrl $ALIGN_DIR/${CTRL}.mapq30.bed \
        --genome $GENOME \
        --gff $GFF \
        --k_min 6 \
        --k_max 20 \
        --k_seqs 600 \
        --outNP \
        --outMEME \
        --outJASPAR \
        --k_neg_dinu_shuffle \
        --t 24 \
        --out $OUTDIR/${SPEC}/${TF}
    done
done

echo "$(date) - GEM peak calling complete for all species!"

# ------------------------------
# Output folder
# ------------------------------
OUTDIR=/path/to/peaks/clean_sorted
mkdir -p "$OUTDIR"

# ------------------------------
# Files per species
# ------------------------------
FILES_AT=(
"/path/to/peaks/peaks_AT/AAA/AAA_outputs/AAA_1_GEM_events.narrowPeak"
"/path/to/peaks/peaks_AT/MAA/MAA_outputs/MAA_1_GEM_events.narrowPeak"
)

FILES_MC=(
"/path/to/peaks/peaks_MC/AMA/AMA_outputs/AMA_1_GEM_events.narrowPeak"
"/path/to/peaks/peaks_MC/MMA/MMA_outputs/MMA_1_GEM_events.narrowPeak"
)

# ------------------------------
# Process Arabidopsis files
# ------------------------------
for f in "${FILES_AT[@]}"; do
    base=$(basename "$f" .narrowPeak)
    out="$OUTDIR/${base}_clean_sorted.narrowPeak"

    echo "Processing AT file $f ..."

    awk 'BEGIN{OFS="\t"}{
        # Remove lowercase "chr" prefix
        sub(/^chr/,"",$1);

        # Extract numeric part of Chr1, Chr2, etc.
        match($1,/Chr([0-9]+)/,a);
        chrnum=a[1];

        # Print original line, prepending numeric value only for sort
        print chrnum,$0
    }' "$f" \
    | sort -k1,1n -k3,3n \
    | cut -f2- > "$out"
done

# ------------------------------
# Process Ice Plant files
# ------------------------------
for f in "${FILES_MC[@]}"; do
    base=$(basename "$f" .narrowPeak)
    out="$OUTDIR/${base}_clean_sorted.narrowPeak"

    echo "Processing MC file $f ..."

    awk 'BEGIN{OFS="\t"}{ $1=substr($1,4); print }' "$f" \
    | awk 'BEGIN{OFS="\t"}{
        split($1,a,"_");
        num=a[2];
        print num,$0
    }' \
    | sort -k1,1n -k3,3n \
    | cut -f2- > "$out"
done

echo "All TF samples processed. Cleaned files are in $OUTDIR"

# ------------------------------
# Load modules
# ------------------------------
module load bedtools/2.30.0

# ------------------------------
# Define paths
# ------------------------------
PEAK_DIR=/path/to/peaks/clean_sorted
OUT_DIR=/path/to/annot
mkdir -p "$OUT_DIR"

# Species-specific gene BED files
GENES_BED_AT=/path/to/Arabidopsis_genes.bed
GENES_BED_MC=/path/to/IcePlant_genes.bed

# ------------------------------
# Sample groups
# ------------------------------
SAMPLES_AT=("AAA" "MAA")
SAMPLES_MC=("MMA" "AMA")

# ------------------------------
# Annotate Arabidopsis peaks
# ------------------------------
for TF in "${SAMPLES_AT[@]}"; do
    PEAK_FILE=${PEAK_DIR}/${TF}_1_GEM_events_clean_sorted.narrowPeak
    OUT_FILE=${OUT_DIR}/${TF}_annotation.bed

    echo "Annotating $PEAK_FILE with nearest Arabidopsis gene..."
    bedtools closest -a "$PEAK_FILE" \
                     -b "$GENES_BED_AT" \
                     -d \
                     -t all \
                     > "$OUT_FILE"
done

# ------------------------------
# Annotate Ice Plant peaks
# ------------------------------
for TF in "${SAMPLES_MC[@]}"; do
    PEAK_FILE=${PEAK_DIR}/${TF}_1_GEM_events_clean_sorted_sorted.narrowPeak
    OUT_FILE=${OUT_DIR}/${TF}_annotation.bed

    echo "Annotating $PEAK_FILE with nearest Ice Plant gene..."
    bedtools closest -a "$PEAK_FILE" \
                     -b "$GENES_BED_MC" \
                     -d \
                     -t all \
                     > "$OUT_FILE"
done

echo "All TF narrowPeak files annotated with closest genes. Results in $OUT_DIR"

# ------------------------------
# Load modules
# ------------------------------
module load samtools
module load bedtools

# ------------------------------
# Directories
# ------------------------------
QC_DIR="/path/to/qc"
ALIGN_DIR="/path/to/bam"
GEM_DIR="/path/to/peaks/clean_sorted"

mkdir -p "$QC_DIR"

# ------------------------------
# Initialize QC summary files
# ------------------------------
echo -e "TF\tTotal_Reads\tReads_in_Peaks\tFRiP" > "$QC_DIR/frip_summary.txt"
echo -e "TF\tNum_Peaks" > "$QC_DIR/peak_counts.txt"
echo -e "TF\tMean_Peak_Length" > "$QC_DIR/peak_length_summary.txt"

# ------------------------------
# List of TFs to process
# ------------------------------
TF_LIST=("AAA" "AMA" "MAA" "MMA")

for TF in "${TF_LIST[@]}"; do
    echo "Running QC metrics for $TF..."

    BAM="$ALIGN_DIR/${TF}.mapq30.bam"
    PEAKS=$(ls "$GEM_DIR/${TF}"*.narrowPeak 2>/dev/null | head -n 1)

    if [[ ! -f "$BAM" ]]; then
        echo "Warning: BAM file $BAM not found, skipping $TF"
        continue
    fi

    if [[ ! -f "$PEAKS" ]]; then
        echo "Warning: Peak file $PEAKS not found, skipping $TF"
        continue
    fi

    # Total reads
    TOTAL_READS=$(samtools view -c "$BAM")

    # Reads in peaks
    READS_IN_PEAKS=$(bedtools intersect -abam "$BAM" -b "$PEAKS" -u | samtools view -c)

    # FRiP calculation
    FRIP=$(awk -v r="$READS_IN_PEAKS" -v t="$TOTAL_READS" 'BEGIN{if(t==0) print "0.0000"; else printf "%.4f", r/t}')

    echo -e "$TF\t$TOTAL_READS\t$READS_IN_PEAKS\t$FRIP" >> "$QC_DIR/frip_summary.txt"

    # Peak counts
    NUM_PEAKS=$(wc -l < "$PEAKS")
    echo -e "$TF\t$NUM_PEAKS" >> "$QC_DIR/peak_counts.txt"

    # Peak lengths
    awk '{print $3-$2}' "$PEAKS" > "$QC_DIR/${TF}_peak_lengths.txt"

    # Mean peak length
    MEAN_PEAK_LENGTH=$(awk '{sum+=$1} END {if(NR>0) print sum/NR; else print 0}' "$QC_DIR/${TF}_peak_lengths.txt")
    echo -e "$TF\t$MEAN_PEAK_LENGTH" >> "$QC_DIR/peak_length_summary.txt"
done

echo "QC metrics completed."

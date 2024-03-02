#!/bin/bash

#########################################################
# Script for standard GATK proccess for variant calling #
#########################################################

###################
# Loading Modules #
###################

module load bwa/0.7.12
module load gatk/4.1.9.0

###################################
# Preparation step  ~ These steps #
# needed to be done only once.    #
###################################

# Downloading Human Ref. Genome (hg38) files #

cd /set/path/to/Ref_Genome  # Path to all prep files 

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz


# Downloading known sites file for base quality score recalibration (bqsr) from GATK resource bundle #

wget -P ./ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ./ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


# Using gatk to create sequence dictionary file #

gatk CreateSequenceDictionary \
-R hg38.fa \
-O hg38_gatk.dict


##########################################################################
# Step 1 ~ Performing Alignment Using Burrows-Wheeler Aligner (BWA) tool #
##########################################################################

# Since we are dealing with multiple files, following loop was used to perform alignment # 

# Here the loop is used for Preterm alignment;replace P with T to perform Term data alignment #

# Define the path to the directory containing reference fasta file and the 'P*' directories #

Ref_Gen="/set/path/to/Ref_Genome/hg38.fa" # Path to hg38 ref genome fasta file

parent_directory="/path/to/WGS/RawData" # Parent RawData directory

known_site="/set/path/to/Ref_Genome/Homo_sapiens_assembly38.dbsnp138.vcf"


# Loop through directories starting with 'P' #

# Here the loop is used for Preterm alignment;replace P with T to perform Term data alignment #


for dir in "$parent_directory"/P*/; do 
    # Check if the directory contains paired-end FASTQ files
    if ls "$dir"/*_1.fq.gz >/dev/null 2>&1 && ls "$dir"/*_2.fq.gz >/dev/null 2>&1; then
        # Perform alignment using BWA
        cd "$dir" || exit 1
        aligned_folder="$dir"/aligned
        mkdir -p "$aligned_folder"

        # Extracting directory name and generating output file name
        dir_name=$(basename "$dir")
        output_sam="$aligned_folder"/"${dir_name}_paired.sam"

        # Running BWA alignment
        bwa mem -t 6 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" "$Ref_Gen" "$dir"/*_1.fq.gz "$dir"/*_2.fq.gz > "$output_sam"
    else
        echo "Paired-end FASTQ files not found in $dir"
    fi
done

########################################################
# Step 2 ~ Deduplication and Sorting .sam file to .bam #
########################################################


# Loop to check for .sam file and perform Duplicate Marking and sorting #

for dir in "$parent_directory"/P*/; do
    # Check if the 'aligned' directory exists inside the 'P' folder
    if [ -d "$dir/aligned" ]; then
        # Loop through SAM files in the 'aligned' directory
        for sam_file in "$dir"/aligned/*.sam; do
            # Check if SAM file exists
            if [ -f "$sam_file" ]; then
                # Get the directory name without extension
                dir_name=$(basename "$dir")
                output_bam="$dir/aligned/${dir_name}_sorted_dedup_reads.bam"
                # Mark duplicates and generate sorted, deduplicated BAM file
                gatk MarkDuplicatesSpark \
                -I "$sam_file" \
                -O "$output_bam"
            else
                echo "SAM file not found in $dir/aligned"
            fi
        done
    else
        echo "aligned directory not found in $dir"
    fi
done


##################################################
# Step 3 ~ Generating recal table and performing #
# bqsr on the deduplicated .bam files            #
##################################################


# Loop through WGS folders

for aligned_dir in "$parent_directory"/*/aligned/; do
    echo "Checking directory: $aligned_dir"

    # Check if the directory contains necessary files
    if [ -d "$aligned_dir" ]; then
        for bam_file in "$aligned_dir"*_sorted_dedup_reads.bam; do
            if [ -e "$bam_file" ]; then
                echo "BAM file found: $bam_file"
                output_recal="${bam_file%_sorted_dedup_reads.bam}_recal_data.table"

                # Building model using BaseRecalibrator
                gatk BaseRecalibrator \
                -I "$bam_file" \
                -R "$Ref_Gen" \
                --known-sites "$known_site" \
                -O "$output_recal"

                # Apply BQSR to adjust base quality score
                output_dedup="${bam_file%_sorted_dedup_reads.bam}_sorted_dedup_bqsr_reads.bam"
                gatk ApplyBQSR -I "$bam_file" \
                -R "$Ref_Gen" \
                --bqsr-recal-file "$output_recal" \
                -O "$output_dedup"
            else
                echo "Necessary files not found in $aligned_dir"
            fi
        done
    else
        echo "aligned directory not found in $aligned_dir"
    fi
done

############################
# Step 4 ~ Variant Calling #
############################

# Iterate through directories

for dir in /home/samaksh/SamakshServerMount/RawData/*/aligned; do
    # Extract P* from the directory path
    dirname=$(basename "$dir")
    value=$(basename $(dirname "$dir"))

    # Create results directory if it doesn't exist
    mkdir -p "$dir"/results

    # Run HaplotypeCaller
    for file in "$dir"/*_sorted_dedup_bqsr_reads.bam; do
        filename=$(basename "$file")
        base=${filename%.*}
        gatk HaplotypeCaller -I "$file" -O "$dir"/results/"$value"_raw_variants.vcf -R "$Ref_Gen"
    done
done

#######################################
# Step 5 ~ Separating SNPs and INDELs #
#######################################


# Loop through directories matching the pattern */aligned/results

for dir in "$parent_directory"/*/aligned/results; do

    # Loop through each file matching *_raw_variants.vcf inside the current directory

    for file in "$dir"/P*_raw_variants.vcf; do

        # Extract the sample name from the file path

        sample=$(basename "$(dirname "$file")")

        # Extract the base name of the file (without extension)

        base_name=$(basename "${file}" _raw_variants.vcf)

        # Debug information (optional)

        echo "Processing file: ${file}"
        echo "Sample: ${sample}"
        echo "Base name: ${base_name}"

        # Select SNPs

        gatk SelectVariants \
            -R "$Ref_Gen" \
            -V "$file" \
            --select-type SNP \
            -O "$dir"/"$base_name"_raw_snps.vcf

        # Select INDELs

        gatk SelectVariants \
            -R "$Ref_Gen" \
            -V "$file" \
            --select-type INDEL \
            -O "$dir"/"$base_name"_raw_indels.vcf

        # Check if output files are generated

        ls "$dir"
    done
done


###########################
# Step 6 ~ SNPs Filtering #
###########################

# Loop through directories matching the pattern */aligned/results #

for dir in ${parent_directory}/*/aligned/results; do

    # Loop through each file matching P*_raw_snps.vcf inside the current directory #

    for file in "${dir}"/P*_raw_snps.vcf; do

        # Run VariantFiltration

        gatk VariantFiltration \
            -R "${Ref_Gen}" \
            -V "${file}" \
            -O "${dir}/$(basename ${file%_raw_snps.vcf})_filtered_snps.vcf" \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 60.0" \
            -filter-name "MQ_filter" -filter "MQ < 40.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"

        # Select Variants that PASS filters

        gatk SelectVariants \
            --exclude-filtered \
            -V "${dir}/$(basename ${file%_raw_snps.vcf})_filtered_snps.vcf" \
            -O "${dir}/$(basename ${file%_raw_snps.vcf})_analysis_snps.vcf"

        # Exclude variants that failed genotype filters

        cat "${dir}/$(basename ${file%_raw_snps.vcf})_analysis_snps.vcf" | grep -v -E "DP_filter|GQ_filter" > "${dir}/$(basename ${file%_raw_snps.vcf})_final_snps.vcf"
    done
done

# Here *_final_snps.vcf will be used for further annotation and analysis #


#########################################
# Step 7 ~ Annotating the SNPs VCF file #
#########################################

ref="/set/path/to/Ref_Genome/hg38.fa"
main="/path/to/WGS/RawData"

for dir in "${main}"/*/aligned; do

    # Extract the sample name from the directory path

    sample_name=$(basename "$(dirname "${dir}")")

    # Loop through each VCF file inside the results directory

    for snp in "${dir}"/results/*_final_snps.vcf; do

        # Construct the correct BAM file path

        bam_file="${dir}/${sample_name}_sorted_dedup_bqsr_reads.bam"

        # Run VariantAnnotator

        gatk VariantAnnotator \
            -R "${ref}" \
            -I "${bam_file}" \
            -V "${snp}" \
            -O "${dir}/${sample_name}_anno.vcf" \
            --dbsnp /home/samaksh/Work/Ref_Genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz 
    done
done
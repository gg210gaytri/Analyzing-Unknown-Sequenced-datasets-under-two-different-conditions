# Import necessary libraries
import gzip
import os
import csv
import subprocess
from Bio import SeqIO


# List of gzipped fastq file paths
fastq_file_paths = [
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode01/condition1_sample0_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode02/condition1_sample1_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode03/condition1_sample2_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode04/condition1_sample3_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode05/condition1_sample4_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode01/condition2_sample0_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode02/condition2_sample1_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode03/condition2_sample2_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode04/condition2_sample3_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode05/condition2_sample4_merged.fastq.gz'
]

# Loop through each fastq file path
for fastq_file_path in fastq_file_paths:
    # Extract the file name
    fastq_file_name = os.path.basename(fastq_file_path)

    # Extract the contents of the gzipped fastq file
    with gzip.open(fastq_file_path, 'rt') as file:
        fastq_file_content = file.read()

# Read contents of true_positives.vcf
true_positives_vcf_path = '/home/gaytri/analyst_test_data/true_positives.vcf'
with open(true_positives_vcf_path, 'r') as file:
    true_positives_content = file.read()

# Function to filter true positives for a given condition
def filter_true_positives_by_condition(vcf_content, condition_name):
    lines = vcf_content.split('\n')
    header = [line for line in lines if line.startswith("#")]
    data_lines = [line for line in lines if not line.startswith("#")]

    # Filter true positives for the given condition
    condition_true_positives = [line for line in data_lines if condition_name in line]

    # Construct VCF content for the filtered true positives
    filtered_vcf_content = "\n".join(header + condition_true_positives)

    return filtered_vcf_content

# Filter true positives for condition1
condition1_true_positives = filter_true_positives_by_condition(true_positives_content, "condition1")

# Filter true positives for condition2
condition2_true_positives = filter_true_positives_by_condition(true_positives_content, "condition2")

# Function to calculate precision, recall, and F1 score
def calculate_precision_recall_f1_for_condition(vcf_content):
    lines = vcf_content.split('\n')
    total_records = len(lines) - 1 

    # Count true positives (PASS) and false positives (FAIL)
    tp = sum(1 for line in lines if "PASS" in line)
    fp = sum(1 for line in lines if "FAIL" in line)

    # Calculate precision, recall, and F1 score
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / total_records
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return precision, recall, f1_score

# Function to calculate true positives (PASS), false positives (FAIL), and false negatives (MISS) for a given condition
def calculate_true_false_positives_negatives_for_condition(vcf_content):
    lines = vcf_content.split('\n')
    total_records = len(lines) - 1  

    # Count true positives (PASS), false positives (FAIL), and false negatives (MISS)
    tp = sum(1 for line in lines if "PASS" in line)
    fp = sum(1 for line in lines if "FAIL" in line)
    fn = total_records - tp

    return tp, fp, fn
# Calculate precision, recall, and F1 score for condition1
precision_condition1, recall_condition1, f1_score_condition1 = calculate_precision_recall_f1_for_condition(condition1_true_positives)

# Calculate precision, recall, and F1 score for condition2
precision_condition2, recall_condition2, f1_score_condition2 = calculate_precision_recall_f1_for_condition(condition2_true_positives)


# Calculate true positives, false positives, and false negatives for condition1
tp_condition1, fp_condition1, fn_condition1 = calculate_true_false_positives_negatives_for_condition(condition1_true_positives)

# Calculate true positives, false positives, and false negatives for condition2
tp_condition2, fp_condition2, fn_condition2 = calculate_true_false_positives_negatives_for_condition(condition2_true_positives)

# Display the results for condition1
print("\nMetrics for Condition 1:")
print(f"True Positives: {tp_condition1}")
print(f"False Positives: {fp_condition1}")
print(f"False Negatives: {fn_condition1}")

# Display the results for condition2
print("\nMetrics for Condition 2:")
print(f"True Positives: {tp_condition2}")
print(f"False Positives: {fp_condition2}")
print(f"False Negatives: {fn_condition2}")

# Display the results for condition1
print("\nMetrics for Condition 1:")
print(f"Precision: {precision_condition1:.4f}")
print(f"Recall: {recall_condition1:.4f}")
print(f"F1 Score: {f1_score_condition1:.4f}")

# Display the results for condition2
print("\nMetrics for Condition 2:")
print(f"Precision: {precision_condition2:.4f}")
print(f"Recall: {recall_condition2:.4f}")
print(f"F1 Score: {f1_score_condition2:.4f}")

# Compare metrics between conditions
print("\nComparison between Conditions:")
print(f"Precision Difference: {precision_condition1 - precision_condition2:.4f}")
print(f"Recall Difference: {recall_condition1 - recall_condition2:.4f}")
print(f"F1 Score Difference: {f1_score_condition1 - f1_score_condition2:.4f}")

# Read contents of reference.fasta
with open('/home/gaytri/analyst_test_data/reference/reference.fasta', 'r') as file:
    reference_fasta_content = file.read()

# Read contents of reference.fasta.fai
with open('/home/gaytri/analyst_test_data/reference/reference.fasta.fai', 'r') as file:
    reference_fai_content = file.read()


# Function to run FastQC on a given fastq file
def run_fastqc(fastq_file_path, output_dir):
    command = f"fastqc -o {output_dir} {fastq_file_path}"
    subprocess.run(command, shell=True)

# Function to compare two FastQC reports and print differences
def compare_fastqc_reports(report1, report2):
    command = f"diff {report1}/summary.txt {report2}/summary.txt"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    print(f"\nDifferences in summary.txt between {report1} and {report2}:")
    print(result.stdout)

# Define output directories for FastQC results
condition1_fastqc_dir = '/home/gaytri/analyst_test_data/test_data/condition1/fastqc_results'
condition2_fastqc_dir = '/home/gaytri/analyst_test_data/test_data/condition2/fastqc_results'

# Create output directories if they don't exist
os.makedirs(condition1_fastqc_dir, exist_ok=True)
os.makedirs(condition2_fastqc_dir, exist_ok=True)

# Run FastQC on the first five files for condition1
for fastq_file_path in fastq_file_paths[:5]:
    run_fastqc(fastq_file_path, condition1_fastqc_dir)

# Run FastQC on the last five files for condition2
for fastq_file_path in fastq_file_paths[5:]:
    run_fastqc(fastq_file_path, condition2_fastqc_dir)

# Compare FastQC reports for condition1 and condition2
compare_fastqc_reports(condition1_fastqc_dir, condition2_fastqc_dir)

# Print a message indicating that FastQC has been run and compared
print("\nFastQC analysis and comparison completed.")


# Define paths to FASTQ files for each condition
condition1_fastq_file_paths = [
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode01/condition1_sample0_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode02/condition1_sample1_merged.fastq.gz', 
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode03/condition1_sample2_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode04/condition1_sample3_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition1/fastqs/barcode05/condition1_sample4_merged.fastq.gz',
]

condition2_fastq_file_paths = [
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode01/condition2_sample0_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode02/condition2_sample1_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode03/condition2_sample2_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode04/condition2_sample3_merged.fastq.gz',
    '/home/gaytri/analyst_test_data/test_data/condition2/fastqs/barcode05/condition2_sample4_merged.fastq.gz'
]

# Define the path to the reference FASTA file
reference_fasta_file_path = '/home/gaytri/analyst_test_data/reference/reference.fasta'

# Define output directories for read metrics, coverage, and sequencing quality
condition1_metrics_dir = '/home/gaytri/analyst_test_data/test_data/condition1/metrics'
condition2_metrics_dir = '/home/gaytri/analyst_test_data/test_data/condition2/metrics'

# Create output directories if they don't exist
os.makedirs(condition1_metrics_dir, exist_ok=True)
os.makedirs(condition2_metrics_dir, exist_ok=True)

# Function to run fastqc
def run_fastqc(fastq_file_path, output_dir):
    command = f"fastqc {fastq_file_path} -o {output_dir}"
    subprocess.run(command, shell=True)

# Function to calculate read metrics using samtools
def calculate_read_metrics(fastq_file_path, output_dir):
    # Create a temporary BAM file path
    temp_bam_path = os.path.join(output_dir, 'temp.bam')

    # Use bwa mem to align the reads and pipe the output to samtools view to create a BAM file
    subprocess.run(f"bwa mem {reference_fasta_file_path} {fastq_file_path} | samtools view -bS - > {temp_bam_path}", shell=True)

    # Calculate read metrics using samtools
    subprocess.run(f"samtools stats {temp_bam_path} > {output_dir}/read_metrics.txt", shell=True)

    # Remove temporary BAM file
    os.remove(temp_bam_path)

# Function to calculate coverage using bedtools
def calculate_coverage(bam_file_path, bed_file_path, output_dir):
    command = f"bedtools genomecov -ibam {bam_file_path} -g {bed_file_path} > {output_dir}/coverage.txt"
    subprocess.run(command, shell=True)

# Function to calculate coverage from FASTQ
def calculate_coverage_from_fastq(fastq_file_path, reference_fasta_file_path, output_dir):
    # Create a temporary BAM file path
    temp_bam_path = os.path.join(output_dir, 'temp.bam')

    # Use bwa mem to align the reads and pipe the output to samtools view to create a BAM file
    subprocess.run(f"bwa mem {reference_fasta_file_path} {fastq_file_path} | samtools view -bS - > {temp_bam_path}", shell=True)

    # Calculate coverage using bedtools
    calculate_coverage(temp_bam_path, reference_fasta_file_path, output_dir)

    # Remove temporary BAM file
    os.remove(temp_bam_path)

# Run read metrics and coverage calculations for condition1
for fastq_file_path in condition1_fastq_file_paths:
    run_fastqc(fastq_file_path, condition1_metrics_dir)
    calculate_read_metrics(fastq_file_path, condition1_metrics_dir)
    calculate_coverage_from_fastq(fastq_file_path, reference_fasta_file_path, condition1_metrics_dir)

# Run read metrics and coverage calculations for condition2
for fastq_file_path in condition2_fastq_file_paths:
    run_fastqc(fastq_file_path, condition2_metrics_dir)
    calculate_read_metrics(fastq_file_path, condition2_metrics_dir)
    calculate_coverage_from_fastq(fastq_file_path, reference_fasta_file_path, condition2_metrics_dir)

# Function to compare read metrics reports and print differences
def compare_read_metrics_reports(report1, report2):
    command = f"diff {report1}/read_metrics.txt {report2}/read_metrics.txt"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    print(f"\nDifferences in read_metrics.txt between {report1} and {report2}:")
    print(result.stdout)

# Compare read metrics reports for condition1 and condition2
compare_read_metrics_reports(condition1_metrics_dir, condition2_metrics_dir)

# Function to compare coverage reports and print differences
def compare_coverage_reports(report1, report2):
    command = f"diff {report1}/coverage.txt {report2}/coverage.txt"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    print(f"\nDifferences in coverage.txt between {report1} and {report2}:")
    print(result.stdout)

# Compare coverage reports for condition1 and condition2
compare_coverage_reports(condition1_metrics_dir, condition2_metrics_dir)

# Print a message indicating that read metrics and coverage analysis has been completed
print("\nRead metrics and coverage analysis completed.")




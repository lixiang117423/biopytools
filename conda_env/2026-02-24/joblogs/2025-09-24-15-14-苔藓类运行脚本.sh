#!/bin/bash

# Simplified WRKY Analysis Pipeline with gffread
# Author: Generated for WRKY transcription factor analysis
# Date: $(date)

# Configuration
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/98.pwrky/08.苔藓类"
FASTA_DIR="/share/org/YZWL/yzwl_lixg/project/98.pwrky/08.苔藓类/01.data/fasta"
GFF_DIR="/share/org/YZWL/yzwl_lixg/project/98.pwrky/08.苔藓类/01.data/gff/fixed"
HMM_DB="/share/org/YZWL/yzwl_lixg/project/98.pwrky/01.data/PF03106.hmm"
CDD_DB="/share/org/YZWL/yzwl_lixg/database/ncbicdd/db/ncbicdd"
WRKY_CDD_ID="460808"
EVALUE_THRESHOLD="1e-4"

# Create directories
mkdir -p "${WORK_DIR}"/{01.protein_seq,02.hmmsearch,03.ncbicdd,04.results}

# Log files
LOG_FILE="${WORK_DIR}/pipeline.log"
ERROR_SAMPLE_FILE="${WORK_DIR}/error_sample.txt"
SUMMARY_FILE="${WORK_DIR}/analysis_summary.txt"

# Logging functions
get_timestamp() {
    if command -v python3 >/dev/null 2>&1; then
        python3 -c "import datetime; print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S,%f')[:-3])"
    else
        date '+%Y-%m-%d %H:%M:%S'
    fi
}

log_info() {
    local message="$1"
    local timestamp=$(get_timestamp)
    echo "${timestamp} - INFO - ${message}" | tee -a "${LOG_FILE}"
}

log_error() {
    local message="$1"
    local timestamp=$(get_timestamp)
    echo "${timestamp} - ERROR - ${message}" | tee -a "${LOG_FILE}"
}

log_warn() {
    local message="$1"
    local timestamp=$(get_timestamp)
    echo "${timestamp} - WARN - ${message}" | tee -a "${LOG_FILE}"
}

# Initialize log
log_info "Simplified WRKY Analysis Pipeline with gffread Started"
log_info "Working directory: ${WORK_DIR}"
log_info "Script PID: $$"

# Clean error file
> "${ERROR_SAMPLE_FILE}"

# Function to extract sample names
get_samples() {
    find "${FASTA_DIR}" -name "*.fa" | sed 's/.*\///; s/\.fa$//' | sort | uniq
}

# Function to extract all protein sequences using gffread
extract_protein_sequences() {
    local sample="$1"
    local fa_file="${FASTA_DIR}/${sample}.fa"
    local gff_file="${GFF_DIR}/${sample}.primary.gff3"
    local output="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
    
    log_info "Extracting all protein sequences for ${sample} using gffread..."
    
    # Check if files exist
    if [[ ! -f "${fa_file}" ]]; then
        log_error "FASTA file not found: ${fa_file}"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
    
    if [[ ! -f "${gff_file}" ]]; then
        log_error "GFF file not found: ${gff_file}"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
    
    # Extract all protein sequences using gffread
    # -y: extract protein sequences
    # -g: genome FASTA file
    # -o: output file
    if gffread -y "${output}" -g "${fa_file}" "${gff_file}" 2>&1; then
        if [[ -s "${output}" ]]; then
            local seq_count=$(grep -c "^>" "${output}" 2>/dev/null || echo 0)
            log_info "Successfully extracted all protein sequences (${seq_count} sequences)"
            return 0
        else
            log_error "Output file is empty"
            echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
            return 1
        fi
    else
        log_error "Failed to extract protein sequences using gffread"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
}

# Function to run HMM search
run_hmm_search() {
    local sample="$1"
    local input="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
    local output="${WORK_DIR}/02.hmmsearch/${sample}.hmm.txt"
    local clean_output="${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt"
    
    log_info "Running HMM search for ${sample}..."
    
    # Check if input file exists
    if [[ ! -f "${input}" ]]; then
        log_error "Input file not found: ${input}"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
    
    if [[ ! -f "${HMM_DB}" ]]; then
        log_error "HMM database not found: ${HMM_DB}"
        return 1
    fi
    
    # Run HMM search
    local error_log="${WORK_DIR}/02.hmmsearch/${sample}.hmm.err"
    if hmmsearch --domtblout "${output}" --cut_tc --cpu 80 "${HMM_DB}" "${input}" > /dev/null 2>"${error_log}"; then
        if [[ -s "${output}" ]]; then
            # Clean results: extract protein_id, ali_from, ali_to, evalue
            grep -v "^#" "${output}" | awk '{
                if (NF >= 22) {
                    protein_id = $1
                    evalue = $7
                    ali_from = $18
                    ali_to = $19
                    print protein_id "\t" ali_from "\t" ali_to "\t" evalue
                }
            }' > "${clean_output}"
            
            local count=$(wc -l < "${clean_output}")
            log_info "HMM search completed successfully - found ${count} domains"
            return 0
        else
            log_warn "No HMM results found (empty output)"
            return 1
        fi
    else
        local exit_code=$?
        log_error "HMM search failed (exit code: ${exit_code})"
        if [[ -s "${error_log}" ]]; then
            log_error "Error details:"
            head -5 "${error_log}" | while read line; do
                log_error "   ${line}"
            done
        fi
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
}

# Function to run CDD search
run_cdd_search() {
    local sample="$1"
    local input="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
    local output="${WORK_DIR}/03.ncbicdd/${sample}.cdd.txt"
    local clean_output="${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt"
    local temp_input="${WORK_DIR}/03.ncbicdd/${sample}.temp.fa"
    
    log_info "Running CDD search for ${sample}..."
    
    # Check if input file exists
    if [[ ! -f "${input}" ]]; then
        log_error "Input file not found: ${input}"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
    
    # Create cleaned temporary file (replace periods with X)
    log_info "Creating temporary file with periods replaced by X..."
    
    local clean_result=$(python3 << EOF 2>&1
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    cleaned_sequences = []
    period_count = 0
    total_sequences = 0
    
    # Read and clean sequences
    for record in SeqIO.parse("${input}", "fasta"):
        total_sequences += 1
        original_seq = str(record.seq)
        
        # Count periods in this sequence
        seq_periods = original_seq.count('.')
        if seq_periods > 0:
            period_count += seq_periods
            # Replace periods with X (unknown amino acid)
            cleaned_seq = original_seq.replace('.', 'X')
            cleaned_record = SeqRecord(Seq(cleaned_seq), id=record.id, description=record.description)
            cleaned_sequences.append(cleaned_record)
        else:
            # No periods, keep original
            cleaned_sequences.append(record)
    
    # Write cleaned sequences to temporary file
    SeqIO.write(cleaned_sequences, "${temp_input}", "fasta")
    print(f"SUCCESS: Cleaned {total_sequences} sequences, replaced {period_count} periods with X")

except Exception as e:
    print(f"ERROR: Failed to create temporary file: {e}")
    exit(1)
EOF
)
    
    if [[ $? -ne 0 ]]; then
        log_error "Failed to create cleaned temporary file"
        log_error "${clean_result}"
        return 1
    fi
    
    log_info "${clean_result}"
    
    # Run CDD search
    local error_log="${WORK_DIR}/03.ncbicdd/${sample}.cdd.err"
    if rpsblast -query "${temp_input}" -outfmt 6 -evalue "${EVALUE_THRESHOLD}" \
                -db "${CDD_DB}" -out "${output}" -num_threads 80 > /dev/null 2>"${error_log}"; then
        # Clean up temporary file
        rm -f "${temp_input}"
        
        if [[ -s "${output}" ]]; then
            # Clean results: filter for WRKY CDD ID and extract relevant columns
            awk -v cdd_id="CDD:${WRKY_CDD_ID}" '$2 == cdd_id {
                print $1 "\t" $7 "\t" $8 "\t" $11
            }' "${output}" > "${clean_output}"
            
            local count=$(wc -l < "${clean_output}")
            log_info "CDD search completed successfully - found ${count} WRKY domains"
            return 0
        else
            log_warn "No CDD results found (empty output)"
            return 1
        fi
    else
        local exit_code=$?
        # Clean up temporary file even on failure
        rm -f "${temp_input}"
        
        log_error "CDD search failed (exit code: ${exit_code})"
        if [[ -s "${error_log}" ]]; then
            log_error "Error details:"
            head -5 "${error_log}" | while read line; do
                log_error "   ${line}"
            done
        fi
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
}

# Function to process single sample
process_sample() {
    local sample="$1"
    log_info ""
    log_info "Processing sample: ${sample}"
    log_info "========================================"
    
    # Step 1: Extract all protein sequences using gffread
    if ! extract_protein_sequences "${sample}"; then
        log_error "Skipping ${sample} due to protein extraction error"
        return 1
    fi
    
    local hmm_success=false
    local cdd_success=false
    
    # Step 2: Run HMM search
    if run_hmm_search "${sample}"; then
        hmm_success=true
        log_info "HMM analysis completed for ${sample}"
    else
        log_warn "HMM analysis failed for ${sample}"
    fi
    
    # Step 3: Run CDD search
    if run_cdd_search "${sample}"; then
        cdd_success=true
        log_info "CDD analysis completed for ${sample}"
    else
        log_warn "CDD analysis failed for ${sample}"
    fi
    
    # Generate simple results summary
    generate_sample_results "${sample}" "${hmm_success}" "${cdd_success}"
    
    if [[ "$hmm_success" == true ]] || [[ "$cdd_success" == true ]]; then
        log_info "Successfully completed ${sample}"
        return 0
    else
        log_error "Both HMM and CDD searches failed for ${sample}"
        return 1
    fi
}

# Function to generate simple results for each sample
generate_sample_results() {
    local sample="$1"
    local hmm_success="$2"
    local cdd_success="$3"
    local results_file="${WORK_DIR}/04.results/${sample}_results.txt"
    
    log_info "Generating results summary for ${sample}..."
    
    {
        echo "# WRKY Analysis Results for ${sample}"
        echo "# Generated at: $(date)"
        echo ""
        
        # Show protein sequence info
        local pep_file="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
        if [[ -f "${pep_file}" ]]; then
            local total_proteins=$(grep -c "^>" "${pep_file}" 2>/dev/null || echo 0)
            echo "## Protein Sequences"
            echo "Total proteins extracted: ${total_proteins}"
            echo ""
        fi
        
        if [[ "$hmm_success" == true ]]; then
            local hmm_file="${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt"
            local hmm_count=$(wc -l < "${hmm_file}" 2>/dev/null || echo 0)
            echo "## HMM Search Results"
            echo "HMM domains found: ${hmm_count}"
            if [[ ${hmm_count} -gt 0 ]]; then
                echo "Protein_ID	Start	End	E-value"
                cat "${hmm_file}"
            fi
            echo ""
        else
            echo "## HMM Search Results"
            echo "HMM search failed or no results found"
            echo ""
        fi
        
        if [[ "$cdd_success" == true ]]; then
            local cdd_file="${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt"
            local cdd_count=$(wc -l < "${cdd_file}" 2>/dev/null || echo 0)
            echo "## CDD Search Results"
            echo "CDD WRKY domains found: ${cdd_count}"
            if [[ ${cdd_count} -gt 0 ]]; then
                echo "Protein_ID	Start	End	E-value"
                cat "${cdd_file}"
            fi
            echo ""
        else
            echo "## CDD Search Results"
            echo "CDD search failed or no results found"
            echo ""
        fi
        
    } > "${results_file}"
    
    log_info "Results summary saved to: ${results_file}"
}

# Function to generate overall summary
generate_summary() {
    log_info ""
    log_info "Generating analysis summary..."
    
    {
        echo "# WRKY Analysis Summary"
        echo "======================="
        echo "Analysis completed at: $(date)"
        echo ""
        
        echo "## Sample Processing Results:"
        echo "----------------------------"
        
        local total_samples=0
        local successful_samples=0
        local error_samples=0
        local total_proteins=0
        local total_hmm_domains=0
        local total_cdd_domains=0
        
        if [[ -f "${ERROR_SAMPLE_FILE}" ]]; then
            error_samples=$(wc -l < "${ERROR_SAMPLE_FILE}" 2>/dev/null || echo 0)
        fi
        
        for sample in $(get_samples); do
            total_samples=$((total_samples + 1))
            if grep -q "${sample}" "${ERROR_SAMPLE_FILE}" 2>/dev/null; then
                echo "FAILED: ${sample}"
            else
                successful_samples=$((successful_samples + 1))
                
                # Count results
                local sample_proteins=0
                local hmm_count=0
                local cdd_count=0
                
                [[ -f "${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa" ]] && \
                    sample_proteins=$(grep -c "^>" "${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa" 2>/dev/null || echo 0)
                
                [[ -f "${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt" ]] && \
                    hmm_count=$(wc -l < "${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt" 2>/dev/null || echo 0)
                
                [[ -f "${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt" ]] && \
                    cdd_count=$(wc -l < "${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt" 2>/dev/null || echo 0)
                
                total_proteins=$((total_proteins + sample_proteins))
                total_hmm_domains=$((total_hmm_domains + hmm_count))
                total_cdd_domains=$((total_cdd_domains + cdd_count))
                
                echo "SUCCESS: ${sample} - Proteins: ${sample_proteins}, HMM domains: ${hmm_count}, CDD domains: ${cdd_count}"
            fi
        done
        
        echo ""
        echo "## Overall Statistics:"
        echo "---------------------"
        echo "Total samples: ${total_samples}"
        echo "Successful samples: ${successful_samples}"
        echo "Failed samples: ${error_samples}"
        if [[ ${total_samples} -gt 0 ]]; then
            echo "Success rate: $(( (successful_samples * 100) / total_samples ))%"
        fi
        echo ""
        echo "Total proteins extracted: ${total_proteins}"
        echo "Total HMM domains found: ${total_hmm_domains}"
        echo "Total CDD WRKY domains found: ${total_cdd_domains}"
        
        if [[ ${error_samples} -gt 0 ]]; then
            echo ""
            echo "## Failed samples:"
            cat "${ERROR_SAMPLE_FILE}" 2>/dev/null
        fi
        
    } > "${SUMMARY_FILE}"
    
    log_info "Summary saved to: ${SUMMARY_FILE}"
}

# Main execution function
main() {
    log_info "Starting simplified WRKY analysis pipeline with gffread..."
    
    # Get all samples
    samples=($(get_samples))
    log_info "Found ${#samples[@]} samples to process"
    log_info "Samples: ${samples[*]}"
    
    # Process samples
    local processed=0
    local successful=0
    local failed=0
    
    for sample in "${samples[@]}"; do
        processed=$((processed + 1))
        log_info ""
        log_info "Processing sample ${processed}/${#samples[@]}: ${sample}"
        
        if process_sample "${sample}"; then
            successful=$((successful + 1))
            log_info "Successfully completed ${sample}"
        else
            failed=$((failed + 1))
            log_warn "Failed to process ${sample}"
        fi
    done
    
    # Generate summary
    generate_summary
    
    log_info ""
    log_info "WRKY analysis pipeline completed!"
    log_info "Processing summary: ${processed} total, ${successful} successful, ${failed} failed"
    log_info "Check results in: ${WORK_DIR}"
    log_info "Summary: ${SUMMARY_FILE}"
    log_info "Individual results: ${WORK_DIR}/04.results/"
    log_info "Analysis finished at: $(date)"
}

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."
    
    local missing_deps=()
    
    command -v gffread >/dev/null 2>&1 || missing_deps+=("gffread")
    command -v hmmsearch >/dev/null 2>&1 || missing_deps+=("hmmsearch")
    command -v rpsblast >/dev/null 2>&1 || missing_deps+=("rpsblast")
    command -v python3 >/dev/null 2>&1 || missing_deps+=("python3")
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        log_error "Missing dependencies: ${missing_deps[*]}"
        exit 1
    fi
    
    # Check Python packages
    if ! python3 -c "import Bio" 2>/dev/null; then
        log_error "Missing Python package: biopython"
        exit 1
    fi
    
    log_info "All dependencies found"
}

# Run the pipeline
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Ensure log directory exists and is writable
    mkdir -p "$(dirname "${LOG_FILE}")"
    
    # Test log file writing
    if ! echo "Log test" > "${LOG_FILE}" 2>/dev/null; then
        echo "Cannot write to log file: ${LOG_FILE}"
        exit 1
    fi
    
    log_info "Starting simplified WRKY pipeline with gffread..."
    
    check_dependencies
    main "$@"
fi
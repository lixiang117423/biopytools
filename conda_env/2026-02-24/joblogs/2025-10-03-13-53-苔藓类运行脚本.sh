# #!/bin/bash

# # Simplified WRKY Analysis Pipeline with gffread
# # Author: Generated for WRKY transcription factor analysis
# # Date: $(date)

# # Configuration
# WORK_DIR="/share/org/YZWL/yzwl_lixg/project/98.wrky/08.苔藓类"
# FASTA_DIR="/share/org/YZWL/yzwl_lixg/project/98.wrky/08.苔藓类/01.data/fasta"
# GFF_DIR="/share/org/YZWL/yzwl_lixg/project/98.wrky/08.苔藓类/01.data/gff/fixed"
# HMM_DB="/share/org/YZWL/yzwl_lixg/project/98.wrky/01.data/PF03106.hmm"
# CDD_DB="/share/org/YZWL/yzwl_lixg/database/ncbicdd/db/ncbicdd"
# WRKY_CDD_ID="460808"
# EVALUE_THRESHOLD="1e-4"

# # Create directories
# mkdir -p "${WORK_DIR}"/{01.protein_seq,02.hmmsearch,03.ncbicdd,04.results}

# # Log files
# LOG_FILE="${WORK_DIR}/pipeline.log"
# ERROR_SAMPLE_FILE="${WORK_DIR}/error_sample.txt"
# SUMMARY_FILE="${WORK_DIR}/analysis_summary.txt"

# # Logging functions
# get_timestamp() {
#     if command -v python3 >/dev/null 2>&1; then
#         python3 -c "import datetime; print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S,%f')[:-3])"
#     else
#         date '+%Y-%m-%d %H:%M:%S'
#     fi
# }

# log_info() {
#     local message="$1"
#     local timestamp=$(get_timestamp)
#     echo "${timestamp} - INFO - ${message}" | tee -a "${LOG_FILE}"
# }

# log_error() {
#     local message="$1"
#     local timestamp=$(get_timestamp)
#     echo "${timestamp} - ERROR - ${message}" | tee -a "${LOG_FILE}"
# }

# log_warn() {
#     local message="$1"
#     local timestamp=$(get_timestamp)
#     echo "${timestamp} - WARN - ${message}" | tee -a "${LOG_FILE}"
# }

# # Initialize log
# log_info "Simplified WRKY Analysis Pipeline with gffread Started"
# log_info "Working directory: ${WORK_DIR}"
# log_info "Script PID: $$"

# # Clean error file
# > "${ERROR_SAMPLE_FILE}"

# # Function to extract sample names
# get_samples() {
#     find "${FASTA_DIR}" -name "*.fa" | sed 's/.*\///; s/\.fa$//' | sort | uniq
# }

# # Function to extract all protein sequences using gffread
# extract_protein_sequences() {
#     local sample="$1"
#     local fa_file="${FASTA_DIR}/${sample}.fa"
#     local gff_file="${GFF_DIR}/${sample}.primary.gff3"
#     local output="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
    
#     log_info "Extracting all protein sequences for ${sample} using gffread..."
    
#     # Check if files exist
#     if [[ ! -f "${fa_file}" ]]; then
#         log_error "FASTA file not found: ${fa_file}"
#         echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#         return 1
#     fi
    
#     if [[ ! -f "${gff_file}" ]]; then
#         log_error "GFF file not found: ${gff_file}"
#         echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#         return 1
#     fi
    
#     # Extract all protein sequences using gffread
#     # -y: extract protein sequences
#     # -g: genome FASTA file
#     # -o: output file
#     if gffread -y "${output}" -g "${fa_file}" "${gff_file}" 2>&1; then
#         if [[ -s "${output}" ]]; then
#             local seq_count=$(grep -c "^>" "${output}" 2>/dev/null || echo 0)
#             log_info "Successfully extracted all protein sequences (${seq_count} sequences)"
#             return 0
#         else
#             log_error "Output file is empty"
#             echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#             return 1
#         fi
#     else
#         log_error "Failed to extract protein sequences using gffread"
#         echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#         return 1
#     fi
# }

# # Function to run HMM search
# run_hmm_search() {
#     local sample="$1"
#     local input="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
#     local output="${WORK_DIR}/02.hmmsearch/${sample}.hmm.txt"
#     local clean_output="${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt"
    
#     log_info "Running HMM search for ${sample}..."
    
#     # Check if input file exists
#     if [[ ! -f "${input}" ]]; then
#         log_error "Input file not found: ${input}"
#         echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#         return 1
#     fi
    
#     if [[ ! -f "${HMM_DB}" ]]; then
#         log_error "HMM database not found: ${HMM_DB}"
#         return 1
#     fi
    
#     # Run HMM search
#     local error_log="${WORK_DIR}/02.hmmsearch/${sample}.hmm.err"
#     if hmmsearch --domtblout "${output}" --cut_tc --cpu 80 "${HMM_DB}" "${input}" > /dev/null 2>"${error_log}"; then
#         if [[ -s "${output}" ]]; then
#             # Clean results: extract protein_id, ali_from, ali_to, evalue
#             grep -v "^#" "${output}" | awk '{
#                 if (NF >= 22) {
#                     protein_id = $1
#                     evalue = $7
#                     ali_from = $18
#                     ali_to = $19
#                     print protein_id "\t" ali_from "\t" ali_to "\t" evalue
#                 }
#             }' > "${clean_output}"
            
#             local count=$(wc -l < "${clean_output}")
#             log_info "HMM search completed successfully - found ${count} domains"
#             return 0
#         else
#             log_warn "No HMM results found (empty output)"
#             return 1
#         fi
#     else
#         local exit_code=$?
#         log_error "HMM search failed (exit code: ${exit_code})"
#         if [[ -s "${error_log}" ]]; then
#             log_error "Error details:"
#             head -5 "${error_log}" | while read line; do
#                 log_error "   ${line}"
#             done
#         fi
#         echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#         return 1
#     fi
# }

# # Function to run CDD search
# run_cdd_search() {
#     local sample="$1"
#     local input="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
#     local output="${WORK_DIR}/03.ncbicdd/${sample}.cdd.txt"
#     local clean_output="${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt"
#     local temp_input="${WORK_DIR}/03.ncbicdd/${sample}.temp.fa"
    
#     log_info "Running CDD search for ${sample}..."
    
#     # Check if input file exists
#     if [[ ! -f "${input}" ]]; then
#         log_error "Input file not found: ${input}"
#         echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#         return 1
#     fi
    
#     # Create cleaned temporary file (replace periods with X)
#     log_info "Creating temporary file with periods replaced by X..."
    
#     local clean_result=$(python3 << EOF 2>&1
# try:
#     from Bio import SeqIO
#     from Bio.Seq import Seq
#     from Bio.SeqRecord import SeqRecord
    
#     cleaned_sequences = []
#     period_count = 0
#     total_sequences = 0
    
#     # Read and clean sequences
#     for record in SeqIO.parse("${input}", "fasta"):
#         total_sequences += 1
#         original_seq = str(record.seq)
        
#         # Count periods in this sequence
#         seq_periods = original_seq.count('.')
#         if seq_periods > 0:
#             period_count += seq_periods
#             # Replace periods with X (unknown amino acid)
#             cleaned_seq = original_seq.replace('.', 'X')
#             cleaned_record = SeqRecord(Seq(cleaned_seq), id=record.id, description=record.description)
#             cleaned_sequences.append(cleaned_record)
#         else:
#             # No periods, keep original
#             cleaned_sequences.append(record)
    
#     # Write cleaned sequences to temporary file
#     SeqIO.write(cleaned_sequences, "${temp_input}", "fasta")
#     print(f"SUCCESS: Cleaned {total_sequences} sequences, replaced {period_count} periods with X")

# except Exception as e:
#     print(f"ERROR: Failed to create temporary file: {e}")
#     exit(1)
# EOF
# )
    
#     if [[ $? -ne 0 ]]; then
#         log_error "Failed to create cleaned temporary file"
#         log_error "${clean_result}"
#         return 1
#     fi
    
#     log_info "${clean_result}"
    
#     # Run CDD search
#     local error_log="${WORK_DIR}/03.ncbicdd/${sample}.cdd.err"
#     if rpsblast -query "${temp_input}" -outfmt 6 -evalue "${EVALUE_THRESHOLD}" \
#                 -db "${CDD_DB}" -out "${output}" -num_threads 80 > /dev/null 2>"${error_log}"; then
#         # Clean up temporary file
#         rm -f "${temp_input}"
        
#         if [[ -s "${output}" ]]; then
#             # Clean results: filter for WRKY CDD ID and extract relevant columns
#             awk -v cdd_id="CDD:${WRKY_CDD_ID}" '$2 == cdd_id {
#                 print $1 "\t" $7 "\t" $8 "\t" $11
#             }' "${output}" > "${clean_output}"
            
#             local count=$(wc -l < "${clean_output}")
#             log_info "CDD search completed successfully - found ${count} WRKY domains"
#             return 0
#         else
#             log_warn "No CDD results found (empty output)"
#             return 1
#         fi
#     else
#         local exit_code=$?
#         # Clean up temporary file even on failure
#         rm -f "${temp_input}"
        
#         log_error "CDD search failed (exit code: ${exit_code})"
#         if [[ -s "${error_log}" ]]; then
#             log_error "Error details:"
#             head -5 "${error_log}" | while read line; do
#                 log_error "   ${line}"
#             done
#         fi
#         echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
#         return 1
#     fi
# }

# # Function to process single sample
# process_sample() {
#     local sample="$1"
#     log_info ""
#     log_info "Processing sample: ${sample}"
#     log_info "========================================"
    
#     # Step 1: Extract all protein sequences using gffread
#     if ! extract_protein_sequences "${sample}"; then
#         log_error "Skipping ${sample} due to protein extraction error"
#         return 1
#     fi
    
#     local hmm_success=false
#     local cdd_success=false
    
#     # Step 2: Run HMM search
#     if run_hmm_search "${sample}"; then
#         hmm_success=true
#         log_info "HMM analysis completed for ${sample}"
#     else
#         log_warn "HMM analysis failed for ${sample}"
#     fi
    
#     # Step 3: Run CDD search
#     if run_cdd_search "${sample}"; then
#         cdd_success=true
#         log_info "CDD analysis completed for ${sample}"
#     else
#         log_warn "CDD analysis failed for ${sample}"
#     fi
    
#     # Generate simple results summary
#     generate_sample_results "${sample}" "${hmm_success}" "${cdd_success}"
    
#     if [[ "$hmm_success" == true ]] || [[ "$cdd_success" == true ]]; then
#         log_info "Successfully completed ${sample}"
#         return 0
#     else
#         log_error "Both HMM and CDD searches failed for ${sample}"
#         return 1
#     fi
# }

# # Function to generate simple results for each sample
# generate_sample_results() {
#     local sample="$1"
#     local hmm_success="$2"
#     local cdd_success="$3"
#     local results_file="${WORK_DIR}/04.results/${sample}_results.txt"
    
#     log_info "Generating results summary for ${sample}..."
    
#     {
#         echo "# WRKY Analysis Results for ${sample}"
#         echo "# Generated at: $(date)"
#         echo ""
        
#         # Show protein sequence info
#         local pep_file="${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa"
#         if [[ -f "${pep_file}" ]]; then
#             local total_proteins=$(grep -c "^>" "${pep_file}" 2>/dev/null || echo 0)
#             echo "## Protein Sequences"
#             echo "Total proteins extracted: ${total_proteins}"
#             echo ""
#         fi
        
#         if [[ "$hmm_success" == true ]]; then
#             local hmm_file="${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt"
#             local hmm_count=$(wc -l < "${hmm_file}" 2>/dev/null || echo 0)
#             echo "## HMM Search Results"
#             echo "HMM domains found: ${hmm_count}"
#             if [[ ${hmm_count} -gt 0 ]]; then
#                 echo "Protein_ID	Start	End	E-value"
#                 cat "${hmm_file}"
#             fi
#             echo ""
#         else
#             echo "## HMM Search Results"
#             echo "HMM search failed or no results found"
#             echo ""
#         fi
        
#         if [[ "$cdd_success" == true ]]; then
#             local cdd_file="${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt"
#             local cdd_count=$(wc -l < "${cdd_file}" 2>/dev/null || echo 0)
#             echo "## CDD Search Results"
#             echo "CDD WRKY domains found: ${cdd_count}"
#             if [[ ${cdd_count} -gt 0 ]]; then
#                 echo "Protein_ID	Start	End	E-value"
#                 cat "${cdd_file}"
#             fi
#             echo ""
#         else
#             echo "## CDD Search Results"
#             echo "CDD search failed or no results found"
#             echo ""
#         fi
        
#     } > "${results_file}"
    
#     log_info "Results summary saved to: ${results_file}"
# }

# # Function to generate overall summary
# generate_summary() {
#     log_info ""
#     log_info "Generating analysis summary..."
    
#     {
#         echo "# WRKY Analysis Summary"
#         echo "======================="
#         echo "Analysis completed at: $(date)"
#         echo ""
        
#         echo "## Sample Processing Results:"
#         echo "----------------------------"
        
#         local total_samples=0
#         local successful_samples=0
#         local error_samples=0
#         local total_proteins=0
#         local total_hmm_domains=0
#         local total_cdd_domains=0
        
#         if [[ -f "${ERROR_SAMPLE_FILE}" ]]; then
#             error_samples=$(wc -l < "${ERROR_SAMPLE_FILE}" 2>/dev/null || echo 0)
#         fi
        
#         for sample in $(get_samples); do
#             total_samples=$((total_samples + 1))
#             if grep -q "${sample}" "${ERROR_SAMPLE_FILE}" 2>/dev/null; then
#                 echo "FAILED: ${sample}"
#             else
#                 successful_samples=$((successful_samples + 1))
                
#                 # Count results
#                 local sample_proteins=0
#                 local hmm_count=0
#                 local cdd_count=0
                
#                 [[ -f "${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa" ]] && \
#                     sample_proteins=$(grep -c "^>" "${WORK_DIR}/01.protein_seq/${sample}.all.pep.fa" 2>/dev/null || echo 0)
                
#                 [[ -f "${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt" ]] && \
#                     hmm_count=$(wc -l < "${WORK_DIR}/02.hmmsearch/${sample}.hmm.clean.txt" 2>/dev/null || echo 0)
                
#                 [[ -f "${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt" ]] && \
#                     cdd_count=$(wc -l < "${WORK_DIR}/03.ncbicdd/${sample}.cdd.clean.txt" 2>/dev/null || echo 0)
                
#                 total_proteins=$((total_proteins + sample_proteins))
#                 total_hmm_domains=$((total_hmm_domains + hmm_count))
#                 total_cdd_domains=$((total_cdd_domains + cdd_count))
                
#                 echo "SUCCESS: ${sample} - Proteins: ${sample_proteins}, HMM domains: ${hmm_count}, CDD domains: ${cdd_count}"
#             fi
#         done
        
#         echo ""
#         echo "## Overall Statistics:"
#         echo "---------------------"
#         echo "Total samples: ${total_samples}"
#         echo "Successful samples: ${successful_samples}"
#         echo "Failed samples: ${error_samples}"
#         if [[ ${total_samples} -gt 0 ]]; then
#             echo "Success rate: $(( (successful_samples * 100) / total_samples ))%"
#         fi
#         echo ""
#         echo "Total proteins extracted: ${total_proteins}"
#         echo "Total HMM domains found: ${total_hmm_domains}"
#         echo "Total CDD WRKY domains found: ${total_cdd_domains}"
        
#         if [[ ${error_samples} -gt 0 ]]; then
#             echo ""
#             echo "## Failed samples:"
#             cat "${ERROR_SAMPLE_FILE}" 2>/dev/null
#         fi
        
#     } > "${SUMMARY_FILE}"
    
#     log_info "Summary saved to: ${SUMMARY_FILE}"
# }

# # Main execution function
# main() {
#     log_info "Starting simplified WRKY analysis pipeline with gffread..."
    
#     # Get all samples
#     samples=($(get_samples))
#     log_info "Found ${#samples[@]} samples to process"
#     log_info "Samples: ${samples[*]}"
    
#     # Process samples
#     local processed=0
#     local successful=0
#     local failed=0
    
#     for sample in "${samples[@]}"; do
#         processed=$((processed + 1))
#         log_info ""
#         log_info "Processing sample ${processed}/${#samples[@]}: ${sample}"
        
#         if process_sample "${sample}"; then
#             successful=$((successful + 1))
#             log_info "Successfully completed ${sample}"
#         else
#             failed=$((failed + 1))
#             log_warn "Failed to process ${sample}"
#         fi
#     done
    
#     # Generate summary
#     generate_summary
    
#     log_info ""
#     log_info "WRKY analysis pipeline completed!"
#     log_info "Processing summary: ${processed} total, ${successful} successful, ${failed} failed"
#     log_info "Check results in: ${WORK_DIR}"
#     log_info "Summary: ${SUMMARY_FILE}"
#     log_info "Individual results: ${WORK_DIR}/04.results/"
#     log_info "Analysis finished at: $(date)"
# }

# # Check dependencies
# check_dependencies() {
#     log_info "Checking dependencies..."
    
#     local missing_deps=()
    
#     command -v gffread >/dev/null 2>&1 || missing_deps+=("gffread")
#     command -v hmmsearch >/dev/null 2>&1 || missing_deps+=("hmmsearch")
#     command -v rpsblast >/dev/null 2>&1 || missing_deps+=("rpsblast")
#     command -v python3 >/dev/null 2>&1 || missing_deps+=("python3")
    
#     if [[ ${#missing_deps[@]} -gt 0 ]]; then
#         log_error "Missing dependencies: ${missing_deps[*]}"
#         exit 1
#     fi
    
#     # Check Python packages
#     if ! python3 -c "import Bio" 2>/dev/null; then
#         log_error "Missing Python package: biopython"
#         exit 1
#     fi
    
#     log_info "All dependencies found"
# }

# # Run the pipeline
# if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
#     # Ensure log directory exists and is writable
#     mkdir -p "$(dirname "${LOG_FILE}")"
    
#     # Test log file writing
#     if ! echo "Log test" > "${LOG_FILE}" 2>/dev/null; then
#         echo "Cannot write to log file: ${LOG_FILE}"
#         exit 1
#     fi
    
#     log_info "Starting simplified WRKY pipeline with gffread..."
    
#     check_dependencies
#     main "$@"
# fi



#!/bin/bash

# WRKY Analysis Pipeline - Serial Version with Classification
# Author: Generated for WRKY transcription factor analysis
# Date: $(date)

# Note: Removed 'set -e' to allow script to continue even when individual samples fail

# 🎯 Configuration
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/98.wrky/08.苔藓类"
FASTA_DIR="/share/org/YZWL/yzwl_lixg/project/98.wrky/08.苔藓类/01.data/fasta"
GFF_DIR="/share/org/YZWL/yzwl_lixg/project/98.wrky/08.苔藓类/01.data/gff/fixed"
HMM_DB="/share/org/YZWL/yzwl_lixg/project/98.wrky/01.data/PF03106.hmm"
CDD_DB="/share/org/YZWL/yzwl_lixg/database/ncbicdd/db/ncbicdd"
WRKY_CDD_ID="460808"
EVALUE_THRESHOLD="1e-4"

# 📁 Create directories
mkdir -p "${WORK_DIR}"/{01.data/longest_pep,02.hmmsearch/{raw,clean,domain_seq,classification},03.ncbicdd/{raw,clean,domain_seq,classification},04.direct/classification,05.domain_seq,06.final_classification}

# 📝 Log files
LOG_FILE="${WORK_DIR}/pipeline.log"
ERROR_SAMPLE_FILE="${WORK_DIR}/01.data/error_sample.txt"
SUMMARY_FILE="${WORK_DIR}/analysis_summary.txt"

# 🕐 时间戳和日志函数
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

log_debug() {
    local message="$1"
    local timestamp=$(get_timestamp)
    echo "${timestamp} - DEBUG - ${message}" | tee -a "${LOG_FILE}"
}

# 🚀 Initialize log
log_info "🧬 WRKY Analysis Pipeline Started - Fixed Version"
log_info "📂 Working directory: ${WORK_DIR}"
log_info "🔧 Script PID: $$"
log_info "🗂️ Log file: ${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

# 🧹 Clean error file
> "${ERROR_SAMPLE_FILE}"

# 📋 Function to extract sample names
get_samples() {
    # Look for .fa files in the fasta directory and extract sample names
    find "${FASTA_DIR}" -name "*.fa" | sed 's/.*\///; s/\.fa$//' | sort | uniq
}

# 🔧 Function to extract longest transcript
extract_longest_transcript() {
    local sample="$1"
    local fa_file="${FASTA_DIR}/${sample}.fa"
    local gff_file="${GFF_DIR}/${sample}.gff"
    local output="${WORK_DIR}/01.data/longest_pep/${sample}.longest.pep.fa"
    
    log_info "📊 Extracting longest transcript for ${sample}..."
    
    # Check if files exist
    if [[ ! -f "${fa_file}" ]]; then
        log_error "❌ FASTA file not found: ${fa_file}"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
    
    if [[ ! -f "${gff_file}" ]]; then
        log_error "❌ GFF file not found: ${gff_file}"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
    
    # Extract longest transcript (files are already uncompressed)
    if biopytools longest-mrna -g "${fa_file}" -f "${gff_file}" -o "${output}" 2>&1; then
        if [[ -s "${output}" ]]; then
            local seq_count=$(grep -c "^>" "${output}" 2>/dev/null || echo 0)
            log_info "✅ Successfully extracted longest transcript (${seq_count} sequences)"
            return 0
        else
            log_error "❌ Output file is empty"
            echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
            return 1
        fi
    else
        log_error "❌ Failed to extract longest transcript"
        echo "${sample}" >> "${ERROR_SAMPLE_FILE}"
        return 1
    fi
}

# 🔍 Function to run HMM search (ORIGINAL WORKING VERSION + MINIMAL DEBUG)
run_hmm_search() {
    local sample="$1"
    local input="${WORK_DIR}/01.data/longest_pep/${sample}.longest.pep.fa"
    local output="${WORK_DIR}/02.hmmsearch/raw/${sample}.hmm.txt"
    
    log_info "🔍 Running HMM search for ${sample}..."
    
    # Add minimal debugging - check files exist
    if [[ ! -f "${input}" ]]; then
        log_error "❌ Input file missing: ${input}"
        return 1
    fi
    if [[ ! -f "${HMM_DB}" ]]; then
        log_error "❌ HMM database missing: ${HMM_DB}"
        return 1
    fi
    
    # Show the command being run (this should appear in log)
    log_info "🔧 Command: hmmsearch --domtblout ${output} --cut_tc --cpu 80 ${HMM_DB} ${input}"
    
    # Original working logic with error capture
    local error_log="${WORK_DIR}/02.hmmsearch/raw/${sample}.hmm.err"
    if hmmsearch --domtblout "${output}" --cut_tc --cpu 80 "${HMM_DB}" "${input}" > /dev/null 2>"${error_log}"; then
        if [[ -s "${output}" ]]; then
            log_info "✅ HMM search completed successfully"
            return 0
        else
            log_warn "⚠️ No HMM results found (empty output)"
            return 1
        fi
    else
        local exit_code=$?
        log_error "❌ HMM search failed (exit code: ${exit_code})"
        if [[ -s "${error_log}" ]]; then
            log_error "❌ Error details:"
            head -5 "${error_log}" | while read line; do
                log_error "   ${line}"
            done
        fi
        return 1
    fi
}

# 🧹 Function to clean HMM results
clean_hmm_results() {
    local sample="$1"
    local input="${WORK_DIR}/02.hmmsearch/raw/${sample}.hmm.txt"
    local output="${WORK_DIR}/02.hmmsearch/clean/${sample}.hmm.clean.txt"
    
    log_info "🧹 Cleaning HMM results for ${sample}..."
    
    # Extract useful information: protein_id, ali_from, ali_to, evalue
    grep -v "^#" "${input}" | awk '{
        if (NF >= 22) {
            protein_id = $1
            evalue = $7
            ali_from = $18
            ali_to = $19
            print protein_id "\t" ali_from "\t" ali_to "\t" evalue
        }
    }' > "${output}"
    
    local count=$(wc -l < "${output}")
    log_info "📊 Found ${count} HMM domains"
}

# 🧬 Function to extract HMM domain sequences
extract_hmm_domains() {
    local sample="$1"
    local clean_file="${WORK_DIR}/02.hmmsearch/clean/${sample}.hmm.clean.txt"
    local pep_file="${WORK_DIR}/01.data/longest_pep/${sample}.longest.pep.fa"
    local output="${WORK_DIR}/02.hmmsearch/domain_seq/${sample}.hmm.domain.fa"
    
    log_info "🧬 Extracting HMM domain sequences for ${sample}..."
    
    local python_result=$(python3 << EOF 2>&1
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    # Read protein sequences
    pep_seqs = {}
    for record in SeqIO.parse("${pep_file}", "fasta"):
        pep_seqs[record.id] = str(record.seq)

    # Extract domains
    domain_seqs = []
    with open("${clean_file}", "r") as f:
        for line in f:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    protein_id, start, end, evalue = parts[:4]
                    start, end = int(start) - 1, int(end)  # Convert to 0-based
                    
                    if protein_id in pep_seqs:
                        domain_seq = pep_seqs[protein_id][start:end]
                        domain_id = f"{protein_id}_hmm_{start+1}_{end}"
                        domain_record = SeqRecord(Seq(domain_seq), id=domain_id, 
                                                description=f"HMM_domain evalue={evalue}")
                        domain_seqs.append(domain_record)

    # Write domain sequences
    if domain_seqs:
        SeqIO.write(domain_seqs, "${output}", "fasta")
        print(f"SUCCESS: Extracted {len(domain_seqs)} HMM domain sequences")
    else:
        with open("${output}", "w") as f:
            pass  # Create empty file
        print("WARNING: No HMM domain sequences found")
        
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)
EOF
)
    
    if [[ $? -eq 0 ]]; then
        log_info "${python_result}"
    else
        log_error "${python_result}"
    fi
}

# 🔬 Function for WRKY subfamily classification
classify_wrky_subfamilies() {
    local sample="$1"
    local input="${WORK_DIR}/02.hmmsearch/domain_seq/${sample}.hmm.domain.fa"
    local output="${WORK_DIR}/02.hmmsearch/classification/${sample}.wrky_classification.txt"
    
    log_info "🔬 Classifying WRKY subfamilies for ${sample}..."
    
    # Check if input file exists and has content
    if [[ ! -f "${input}" ]] || [[ ! -s "${input}" ]]; then
        log_warn "⚠️ No HMM domain sequences available for classification"
        echo -e "Gene_ID\tGroup\tSubgroup\tConfidence\tNote" > "${output}"
        echo -e "No_domains\tN/A\tN/A\tN/A\tNo HMM domain sequences found" >> "${output}"
        return 1
    fi
    
    local python_result=$(python3 << EOF 2>&1
import re
from typing import Dict, List, Tuple
from collections import defaultdict

def quick_wrky_classifier(fasta_file: str) -> Dict:
    """
    快速WRKY分类器 - 保守但准确
    """
    # 读取FASTA文件
    sequences = read_fasta_file(fasta_file)
    
    # 按基因ID分组
    gene_groups = group_sequences_by_gene(sequences)
    
    results = {}
    for gene_id, domains in gene_groups.items():
        classification = classify_gene_quick(domains)
        # 处理Group I的特殊情况
        if classification['group'] == 'I' and 'domains' in classification:
            # 为每个域创建单独的记录
            for i, domain_info in enumerate(classification['domains']):
                domain_gene_id = f"{gene_id}_domain_{i+1}"
                results[domain_gene_id] = {
                    'group': 'I',
                    'subgroup': domain_info['subgroup'],
                    'confidence': domain_info['confidence'],
                    'reason': classification.get('reason', '-')
                }
        else:
            results[gene_id] = classification
    
    return results

def read_fasta_file(file_path: str) -> List[Tuple[str, str]]:
    """读取FASTA文件"""
    sequences = []
    try:
        with open(file_path, 'r') as f:
            header = None
            sequence = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if header is not None:
                        sequences.append((header, ''.join(sequence)))
                    header = line
                    sequence = []
                else:
                    sequence.append(line)
            # 添加最后一条序列
            if header is not None:
                sequences.append((header, ''.join(sequence)))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
    return sequences

def group_sequences_by_gene(sequences: List[Tuple[str, str]]) -> Dict:
    """按基因ID分组序列"""
    gene_groups = defaultdict(list)
    
    for header, seq in sequences:
        gene_id = extract_gene_id(header)
        position = extract_position(header)
        gene_groups[gene_id].append({
            'header': header,
            'sequence': seq,
            'position': position
        })
    
    # 按位置排序每个基因的结构域
    for gene_id in gene_groups:
        gene_groups[gene_id].sort(key=lambda x: x['position'])
    
    return dict(gene_groups)

def extract_gene_id(header: str) -> str:
    """从header提取基因ID"""
    return header.split('_hmm_')[0].replace('>', '')

def extract_position(header: str) -> int:
    """从header提取位置信息"""
    match = re.search(r'_hmm_(\d+)_', header)
    return int(match.group(1)) if match else 0

def classify_gene_quick(domains: List[Dict]) -> Dict:
    """
    快速基因分类 - 基于明确特征
    """
    domain_count = len(domains)
    
    if domain_count == 2:
        # Group I: 两个结构域，返回两个分类结果
        return classify_group_i(domains)
    elif domain_count == 1:
        # Group II/III: 单个结构域
        return classify_single_domain(domains[0])
    else:
        return {
            'group': 'other',
            'subgroup': '-',
            'reason': f'unusual_domain_count_{domain_count}',
            'confidence': 'low'
        }

def classify_group_i(domains: List[Dict]) -> Dict:
    """
    Group I 分类 - 返回两个域的分类结果
    """
    n_term = domains[0]['sequence']
    c_term = domains[1]['sequence']
    
    # 检查关键特征
    n_wrky = extract_wrky_start(n_term)
    c_wrky = extract_wrky_start(c_term)
    
    # 简单明确的规则
    if n_wrky == 'N' and c_wrky == 'R':
        # 进一步验证锌指
        n_zinc = check_c2h2_zinc_finger(n_term)
        c_zinc = check_c2h2_zinc_finger(c_term)
        
        if n_zinc and c_zinc:
            return {
                'group': 'I',
                'domains': [
                    {
                        'subgroup': 'I_NT',
                        'confidence': 'high',
                        'position': domains[0]['position']
                    },
                    {
                        'subgroup': 'I_CT',
                        'confidence': 'high',
                        'position': domains[1]['position']
                    }
                ],
                'reason': '-'
            }
    
    return {
        'group': 'other',
        'subgroup': '-',
        'reason': f'group_i_unclear_n={n_wrky}_c={c_wrky}',
        'confidence': 'low'
    }

def classify_single_domain(domain: Dict) -> Dict:
    """
    单结构域分类 - 基于明确特征
    """
    seq = domain['sequence']
    wrky_start = extract_wrky_start(seq)
    
    # 检查锌指类型
    has_c2h2 = check_c2h2_zinc_finger(seq)
    has_c2hc = check_c2hc_zinc_finger(seq)
    
    if not (has_c2h2 or has_c2hc):
        return {
            'group': 'other',
            'subgroup': '-',
            'reason': 'no_clear_zinc_finger',
            'confidence': 'low'
        }
    
    # Group III: C2HC锌指
    if has_c2hc:
        if wrky_start == 'S':
            return {'group': 'III', 'subgroup': 'IIIa', 'confidence': 'high', 'reason': '-'}
        elif wrky_start == 'Q':
            return {'group': 'III', 'subgroup': 'IIIb', 'confidence': 'high', 'reason': '-'}
        else:
            return {
                'group': 'other', 
                'subgroup': '-',
                'reason': f'group_iii_unclear_start={wrky_start}',
                'confidence': 'low'
            }
    
    # Group II: C2H2锌指
    if has_c2h2:
        return classify_group_ii(seq, wrky_start)
    
    return {
        'group': 'other',
        'subgroup': '-',
        'reason': 'zinc_finger_ambiguous',
        'confidence': 'low'
    }

def classify_group_ii(seq: str, wrky_start: str) -> Dict:
    """
    Group II 详细分类 - 基于明确模式
    """
    # IIa/IIb: 特殊基序识别
    if wrky_start == 'Q':
        if 'CPVKKKV' in seq:
            return {'group': 'II', 'subgroup': 'IIa', 'confidence': 'high', 'reason': '-'}
        elif 'CPVRKQV' in seq:
            return {'group': 'II', 'subgroup': 'IIb', 'confidence': 'high', 'reason': '-'}
        else:
            return {
                'group': 'other',
                'subgroup': '-',
                'reason': 'group_ii_q_no_special_motif',
                'confidence': 'low'
            }
    
    # IIc: R开头 + 特定特征
    elif wrky_start == 'R':
        if 'RSYYRC' in seq:  # vs RSYYKC in Group I
            return {'group': 'II', 'subgroup': 'IIc', 'confidence': 'medium', 'reason': '-'}
        else:
            return {
                'group': 'other',
                'subgroup': '-',
                'reason': 'group_ii_r_unclear',
                'confidence': 'low'
            }
    
    # IId: S开头
    elif wrky_start == 'S':
        if check_x5_spacing(seq):
            return {'group': 'II', 'subgroup': 'IId', 'confidence': 'medium', 'reason': '-'}
        else:
            return {
                'group': 'other',
                'subgroup': '-',
                'reason': 'group_ii_s_unclear_spacing',
                'confidence': 'low'
            }
    
    # IIe: A开头
    elif wrky_start == 'A':
        if check_x5_spacing(seq):
            return {'group': 'II', 'subgroup': 'IIe', 'confidence': 'medium', 'reason': '-'}
        else:
            return {
                'group': 'other',
                'subgroup': '-',
                'reason': 'group_ii_a_unclear_spacing', 
                'confidence': 'low'
            }
    
    else:
        return {
            'group': 'other',
            'subgroup': '-',
            'reason': f'group_ii_unknown_start={wrky_start}',
            'confidence': 'low'
        }

# 辅助函数
def extract_wrky_start(sequence: str) -> str:
    """提取WRKY前的起始氨基酸"""
    match = re.search(r'([NRQSAW])WRKYGQK', sequence)
    return match.group(1) if match else 'unknown'

def check_c2h2_zinc_finger(sequence: str) -> bool:
    """检查C2H2锌指模式"""
    pattern = r'C.{4,5}C.{20,25}H.H'
    return bool(re.search(pattern, sequence))

def check_c2hc_zinc_finger(sequence: str) -> bool:
    """检查C2HC锌指模式"""
    pattern = r'C.{6,8}C.{20,30}H.C'
    return bool(re.search(pattern, sequence))

def check_x5_spacing(sequence: str) -> bool:
    """检查X5间隔模式 (IId/IIe特征)"""
    pattern = r'C.{5}C.{20,25}H'
    return bool(re.search(pattern, sequence))

# 主处理函数
def process_hmm_results(file_path: str, output_path: str) -> None:
    """
    处理hmmsearch结果并输出分类
    """
    try:
        results = quick_wrky_classifier(file_path)
        
        # 统计信息
        total_genes = len(results)
        group_counts = defaultdict(int)
        confidence_counts = defaultdict(int)
        
        # 写入结果文件
        with open(output_path, 'w') as f:
            f.write("Gene_ID\tGroup\tSubgroup\tConfidence\tNote\n")
            
            for gene_id, classification in results.items():
                group = classification['group']
                subgroup = classification.get('subgroup', '-')
                confidence = classification['confidence']
                reason = classification.get('reason', '-')
                
                f.write(f"{gene_id}\t{group}\t{subgroup}\t{confidence}\t{reason}\n")
                
                # 统计
                if group in ['I', 'II', 'III']:
                    group_counts[subgroup] += 1
                else:
                    group_counts['other'] += 1
                confidence_counts[confidence] += 1
        
        # 输出统计信息
        print(f"SUCCESS: Classification completed: {total_genes} genes processed")
        print("Group distribution:")
        for group, count in sorted(group_counts.items()):
            print(f"  {group}: {count}")
        print("Confidence distribution:")
        for conf, count in sorted(confidence_counts.items()):
            print(f"  {conf}: {count}")
            
    except Exception as e:
        print(f"ERROR: WRKY classification failed: {e}")

# 执行分类
process_hmm_results("${input}", "${output}")
EOF
)
    log_info "${python_result}"
}

# 🔍 Function to run CDD search (ORIGINAL WORKING VERSION + MINIMAL DEBUG + PERIOD HANDLING)
run_cdd_search() {
    local sample="$1"
    local input="${WORK_DIR}/01.data/longest_pep/${sample}.longest.pep.fa"
    local output="${WORK_DIR}/03.ncbicdd/raw/${sample}.cdd.txt"
    local temp_input="${WORK_DIR}/03.ncbicdd/raw/${sample}.temp.fa"
    
    log_info "🔍 Running CDD search for ${sample}..."
    
    # Add minimal debugging - check files exist
    if [[ ! -f "${input}" ]]; then
        log_error "❌ Input file missing: ${input}"
        return 1
    fi
    
    # Create a cleaned temporary file - replace periods with X (unknown amino acid)
    log_info "🔧 Creating temporary file with periods replaced by X..."
    
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
        log_error "❌ Failed to create cleaned temporary file"
        log_error "${clean_result}"
        return 1
    fi
    
    log_info "${clean_result}"
    
    # Show the command being run (using temporary file)
    log_info "🔧 Command: rpsblast -query ${temp_input} -outfmt 6 -evalue ${EVALUE_THRESHOLD} -db ${CDD_DB} -out ${output} -num_threads 80"
    
    # Run CDD search with cleaned temporary file
    local error_log="${WORK_DIR}/03.ncbicdd/raw/${sample}.cdd.err"
    if rpsblast -query "${temp_input}" -outfmt 6 -evalue "${EVALUE_THRESHOLD}" \
                -db "${CDD_DB}" -out "${output}" -num_threads 80 > /dev/null 2>"${error_log}"; then
        # Clean up temporary file
        rm -f "${temp_input}"
        
        if [[ -s "${output}" ]]; then
            log_info "✅ CDD search completed successfully"
            return 0
        else
            log_warn "⚠️ No CDD results found (empty output)"
            return 1
        fi
    else
        local exit_code=$?
        # Clean up temporary file even on failure
        rm -f "${temp_input}"
        
        log_error "❌ CDD search failed (exit code: ${exit_code})"
        if [[ -s "${error_log}" ]]; then
            log_error "❌ Error details:"
            head -5 "${error_log}" | while read line; do
                log_error "   ${line}"
            done
        fi
        return 1
    fi
}

# 🧹 Function to clean CDD results
clean_cdd_results() {
    local sample="$1"
    local input="${WORK_DIR}/03.ncbicdd/raw/${sample}.cdd.txt"
    local output="${WORK_DIR}/03.ncbicdd/clean/${sample}.cdd.clean.txt"
    
    log_info "🧹 Cleaning CDD results for ${sample}..."
    
    # Filter for WRKY CDD ID (note: CDD results have "CDD:" prefix)
    # Extract: protein_id, query_start, query_end, evalue
    awk -v cdd_id="CDD:${WRKY_CDD_ID}" '$2 == cdd_id {
        print $1 "\t" $7 "\t" $8 "\t" $11
    }' "${input}" > "${output}"
    
    local count=$(wc -l < "${output}")
    log_info "📊 Found ${count} CDD domains (CDD:${WRKY_CDD_ID})"
}

# 🧬 Function to extract CDD domain sequences
extract_cdd_domains() {
    local sample="$1"
    local clean_file="${WORK_DIR}/03.ncbicdd/clean/${sample}.cdd.clean.txt"
    local pep_file="${WORK_DIR}/01.data/longest_pep/${sample}.longest.pep.fa"
    local output="${WORK_DIR}/03.ncbicdd/domain_seq/${sample}.ncbicdd.domain.fa"
    
    log_info "🧬 Extracting CDD domain sequences for ${sample}..."
    
    local python_result=$(python3 << EOF 2>&1
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    # Read protein sequences
    pep_seqs = {}
    for record in SeqIO.parse("${pep_file}", "fasta"):
        pep_seqs[record.id] = str(record.seq)

    # Extract domains
    domain_seqs = []
    with open("${clean_file}", "r") as f:
        for line in f:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    protein_id, start, end, evalue = parts[:4]
                    start, end = int(start) - 1, int(end)  # Convert to 0-based
                    
                    if protein_id in pep_seqs:
                        domain_seq = pep_seqs[protein_id][start:end]
                        domain_id = f"{protein_id}_cdd_{start+1}_{end}"
                        domain_record = SeqRecord(Seq(domain_seq), id=domain_id, 
                                                description=f"CDD_domain evalue={evalue}")
                        domain_seqs.append(domain_record)

    # Write domain sequences
    if domain_seqs:
        SeqIO.write(domain_seqs, "${output}", "fasta")
        print(f"SUCCESS: Extracted {len(domain_seqs)} CDD domain sequences")
    else:
        with open("${output}", "w") as f:
            pass  # Create empty file
        print("WARNING: No CDD domain sequences found")
        
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)
EOF
)
    
    if [[ $? -eq 0 ]]; then
        log_info "${python_result}"
    else
        log_error "${python_result}"
    fi
}

# 🔄 Function to integrate results
integrate_results() {
    local sample="$1"
    local hmm_file="${WORK_DIR}/02.hmmsearch/clean/${sample}.hmm.clean.txt"
    local cdd_file="${WORK_DIR}/03.ncbicdd/clean/${sample}.cdd.clean.txt"
    local hmm_seq="${WORK_DIR}/02.hmmsearch/domain_seq/${sample}.hmm.domain.fa"
    local cdd_seq="${WORK_DIR}/03.ncbicdd/domain_seq/${sample}.ncbicdd.domain.fa"
    
    log_info "🔄 Integrating results for ${sample}..."
    
    local python_result=$(python3 << EOF 2>&1
import os
from collections import defaultdict
import shutil

try:
    # Read HMM results
    hmm_domains = {}
    if os.path.exists("${hmm_file}") and os.path.getsize("${hmm_file}") > 0:
        with open("${hmm_file}", "r") as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        protein_id, start, end, evalue = parts[:4]
                        hmm_domains[protein_id] = hmm_domains.get(protein_id, [])
                        hmm_domains[protein_id].append((int(start), int(end)))

    # Read CDD results
    cdd_domains = {}
    if os.path.exists("${cdd_file}") and os.path.getsize("${cdd_file}") > 0:
        with open("${cdd_file}", "r") as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        protein_id, start, end, evalue = parts[:4]
                        cdd_domains[protein_id] = cdd_domains.get(protein_id, [])
                        cdd_domains[protein_id].append((int(start), int(end)))

    # Compare and integrate
    final_output = "${WORK_DIR}/05.domain_seq/${sample}.final.domain.fa"

    if os.path.exists("${hmm_seq}") and os.path.exists("${cdd_seq}"):
        # Combine both files
        with open(final_output, "w") as out:
            with open("${hmm_seq}", "r") as f:
                out.write(f.read())
            with open("${cdd_seq}", "r") as f:
                out.write(f.read())
        print("SUCCESS: Combined HMM and CDD results")
    elif os.path.exists("${hmm_seq}"):
        shutil.copy("${hmm_seq}", final_output)
        print("SUCCESS: Used HMM results only")
    elif os.path.exists("${cdd_seq}"):
        shutil.copy("${cdd_seq}", final_output)
        print("SUCCESS: Used CDD results only")
    else:
        with open(final_output, "w") as f:
            pass  # Create empty file
        print("WARNING: No domain results available")
        
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)
EOF
)
    
    if [[ $? -eq 0 ]]; then
        log_info "${python_result}"
    else
        log_error "${python_result}"
    fi
}

# 🔤 Function to search WRKY directly in sequences
search_wrky_direct() {
    local sample="$1"
    local input="${WORK_DIR}/01.data/longest_pep/${sample}.longest.pep.fa"
    local output="${WORK_DIR}/04.direct/${sample}.wrky_direct.fa"
    
    log_info "🔤 Searching WRKY directly in sequences for ${sample}..."
    
    local python_result=$(python3 << EOF 2>&1
from Bio import SeqIO
import re

try:
    wrky_sequences = []
    for record in SeqIO.parse("${input}", "fasta"):
        # Search for WRKY in the protein sequence
        if re.search(r'WRKY', str(record.seq), re.IGNORECASE):
            wrky_sequences.append(record)

    if wrky_sequences:
        SeqIO.write(wrky_sequences, "${output}", "fasta")
        print(f"SUCCESS: Found {len(wrky_sequences)} sequences containing WRKY")
    else:
        with open("${output}", "w") as f:
            pass  # Create empty file
        print("WARNING: No sequences containing WRKY found")
        
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)
EOF
)
    
    if [[ $? -eq 0 ]]; then
        log_info "${python_result}"
        return 0
    else
        log_error "${python_result}"
        return 1
    fi
}

# 🔬 Function to generate final WRKY classification summary
generate_final_classification() {
    local sample="$1"
    log_info "🔬 Generating final WRKY classification for ${sample}..."
    
    local hmm_class="${WORK_DIR}/02.hmmsearch/classification/${sample}.wrky_classification.txt"
    local cdd_seq="${WORK_DIR}/03.ncbicdd/domain_seq/${sample}.ncbicdd.domain.fa"
    local direct_seq="${WORK_DIR}/04.direct/${sample}.wrky_direct.fa"
    local final_class="${WORK_DIR}/06.final_classification/${sample}.final_wrky_summary.txt"
    
    local python_result=$(python3 << EOF 2>&1
import os
from collections import defaultdict
import re

def extract_gene_id_from_header(header):
    """从fasta header中提取基因ID"""
    gene_id = header.replace('>', '').split()[0]
    gene_id = re.split(r'_(?:hmm|cdd|domain)_', gene_id)[0]
    return gene_id

def read_fasta_gene_ids(file_path):
    """读取FASTA文件并返回基因ID集合"""
    gene_ids = set()
    if not os.path.exists(file_path) or not os.path.getsize(file_path):
        return gene_ids
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    gene_id = extract_gene_id_from_header(line.strip())
                    gene_ids.add(gene_id)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    
    return gene_ids

try:
    # 读取HMM分类结果
    hmm_classified = {}
    if os.path.exists("${hmm_class}") and os.path.getsize("${hmm_class}") > 0:
        with open("${hmm_class}", 'r') as f:
            next(f)  # 跳过header
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 5 and parts[0] != 'No_domains':
                        gene_id, group, subgroup, confidence, note = parts[:5]
                        # 处理Group I的特殊情况
                        base_gene_id = re.split(r'_domain_', gene_id)[0]
                        hmm_classified[base_gene_id] = {
                            'group': group,
                            'subgroup': subgroup,
                            'confidence': confidence,
                            'method': 'HMM'
                        }
    
    # 读取CDD序列基因ID
    cdd_genes = read_fasta_gene_ids("${cdd_seq}")
    
    # 读取Direct搜索结果基因ID
    direct_genes = read_fasta_gene_ids("${direct_seq}")
    
    # 获取所有基因ID
    hmm_genes = set(hmm_classified.keys())
    all_genes = hmm_genes | cdd_genes | direct_genes
    
    # 生成最终汇总
    stats = {
        'total_wrky_genes': len(all_genes),
        'hmm_classified': len(hmm_genes),
        'cdd_found': len(cdd_genes),
        'direct_found': len(direct_genes),
        'both_hmm_direct': len(hmm_genes & direct_genes),
        'both_cdd_direct': len(cdd_genes & direct_genes),
        'both_hmm_cdd': len(hmm_genes & cdd_genes),
        'all_three': len(hmm_genes & cdd_genes & direct_genes),
        'direct_only': len(direct_genes - hmm_genes - cdd_genes)
    }
    
    # 分组统计
    group_stats = defaultdict(int)
    for gene_data in hmm_classified.values():
        if gene_data['group'] in ['I', 'II', 'III']:
            group_stats[gene_data['subgroup']] += 1
        else:
            group_stats['other'] += 1
    
    # 写入最终分类文件
    with open("${final_class}", 'w') as f:
        f.write("# WRKY Gene Classification Summary for ${sample}\n")
        f.write("# Generated by WRKY Analysis Pipeline\n")
        f.write("#\n")
        f.write("Gene_ID\tGroup\tSubgroup\tConfidence\tMethod\tDirect_Search\n")
        
        for gene_id in sorted(all_genes):
            if gene_id in hmm_classified:
                gene_data = hmm_classified[gene_id]
                direct_found = "Yes" if gene_id in direct_genes else "No"
                subgroup = gene_data['subgroup']
                
                f.write(f"{gene_id}\t{gene_data['group']}\t{subgroup}\t"
                       f"{gene_data['confidence']}\t{gene_data['method']}\t{direct_found}\n")
            else:
                f.write(f"{gene_id}\tNA\tNA\tNA\tDirect_only\tYes\n")
        
        f.write("\n# Statistics Summary:\n")
        f.write(f"# Total WRKY genes found: {stats['total_wrky_genes']}\n")
        f.write(f"# HMM classified: {stats['hmm_classified']}\n")
        f.write(f"# CDD domain found: {stats['cdd_found']}\n")
        f.write(f"# Direct search found: {stats['direct_found']}\n")
        f.write(f"# HMM + Direct: {stats['both_hmm_direct']}\n")
        f.write(f"# CDD + Direct: {stats['both_cdd_direct']}\n")
        f.write(f"# HMM + CDD: {stats['both_hmm_cdd']}\n")
        f.write(f"# All three methods: {stats['all_three']}\n")
        f.write(f"# Direct only: {stats['direct_only']}\n")
        f.write("\n# Group Distribution:\n")
        for group, count in sorted(group_stats.items()):
            f.write(f"# {group}: {count}\n")
    
    # 输出统计信息
    print(f"SUCCESS: Final classification summary:")
    print(f"  Total WRKY genes: {stats['total_wrky_genes']}")
    print(f"  HMM classified: {stats['hmm_classified']}")
    print(f"  CDD domain found: {stats['cdd_found']}")
    print(f"  Direct search found: {stats['direct_found']}")
    print(f"  HMM + Direct: {stats['both_hmm_direct']}")
    print(f"  CDD + Direct: {stats['both_cdd_direct']}")
    print(f"  HMM + CDD: {stats['both_hmm_cdd']}")
    print(f"  All three methods: {stats['all_three']}")
    print(f"  Direct only: {stats['direct_only']}")
    if group_stats:
        print(f"  Group distribution: {dict(group_stats)}")

except Exception as e:
    print(f"ERROR: Final classification failed: {e}")
    sys.exit(1)
EOF
)
    
    if [[ $? -eq 0 ]]; then
        log_info "${python_result}"
        return 0
    else
        log_error "${python_result}"
        return 1
    fi
}

# 📊 Function to process single sample
process_sample() {
    local sample="$1"
    log_info ""
    log_info "🧬 Processing sample: ${sample}"
    log_info "========================================"
    
    # Step 1: Extract longest transcript - this is critical, if it fails, skip the sample
    if ! extract_longest_transcript "${sample}"; then
        log_error "❌ Skipping ${sample} due to extraction error"
        return 1
    fi
    
    # Step 2: HMM search and classification - continue even if this fails
    local hmm_success=false
    if run_hmm_search "${sample}"; then
        clean_hmm_results "${sample}" && hmm_success=true
        extract_hmm_domains "${sample}"
        classify_wrky_subfamilies "${sample}"
        log_info "✅ HMM analysis completed for ${sample}"
    else
        log_warn "⚠️ HMM analysis failed for ${sample}, continuing with other methods"
    fi
    
    # Step 3: CDD search - continue even if this fails  
    local cdd_success=false
    if run_cdd_search "${sample}"; then
        clean_cdd_results "${sample}" && cdd_success=true
        extract_cdd_domains "${sample}"
        log_info "✅ CDD analysis completed for ${sample}"
    else
        log_warn "⚠️ CDD analysis failed for ${sample}, continuing with other methods"
    fi
    
    # Step 4: Integrate results - only if we have some results
    if [[ "$hmm_success" == true ]] || [[ "$cdd_success" == true ]]; then
        integrate_results "${sample}" || log_warn "⚠️ Result integration had issues for ${sample}"
    else
        log_warn "⚠️ No HMM or CDD results to integrate for ${sample}"
    fi
    
    # Step 5: Direct WRKY search - this should always work if we have the protein file
    if search_wrky_direct "${sample}"; then
        log_info "✅ Direct WRKY search completed for ${sample}"
    else
        log_warn "⚠️ Direct WRKY search failed for ${sample}"
    fi
    
    # Step 6: Generate final classification summary
    if generate_final_classification "${sample}"; then
        log_info "✅ Final classification generated for ${sample}"
    else
        log_warn "⚠️ Final classification failed for ${sample}"
    fi
    
    log_info "✅ Completed processing ${sample} (some steps may have failed)"
    log_info "========================================"
    return 0  # Always return success if we got past the initial extraction step
}

# 在脚本开头添加安全的算术运算函数
safe_add() {
    local var1="$1"
    local var2="$2"
    
    # 清理变量，移除空格和非数字字符
    var1=$(echo "$var1" | tr -d ' ' | grep -o '[0-9]*' | head -1)
    var2=$(echo "$var2" | tr -d ' ' | grep -o '[0-9]*' | head -1)
    
    # 如果为空则设为0
    var1=${var1:-0}
    var2=${var2:-0}
    
    # 验证是否为数字
    if ! [[ "$var1" =~ ^[0-9]+$ ]]; then
        var1=0
    fi
    if ! [[ "$var2" =~ ^[0-9]+$ ]]; then
        var2=0
    fi
    
    echo $((var1 + var2))
}

# 📈 Function to generate summary
generate_summary() {
    log_info ""
    log_info "📈 Generating analysis summary..."
    
    {
        echo "🧬 WRKY Analysis Summary"
        echo "======================="
        echo "📅 Analysis completed at: $(date)"
        echo ""
        
        echo "📊 Sample Processing Results:"
        echo "----------------------------"
        
        local total_samples=0
        local successful_samples=0
        local error_samples=0
        
        if [[ -f "${ERROR_SAMPLE_FILE}" ]]; then
            error_samples=$(wc -l < "${ERROR_SAMPLE_FILE}" 2>/dev/null || echo 0)
        fi
        
        for sample in $(get_samples); do
            total_samples=$((total_samples + 1))
            if grep -q "${sample}" "${ERROR_SAMPLE_FILE}" 2>/dev/null; then
                echo "❌ ${sample}: Failed in transcript extraction"
            else
                successful_samples=$((successful_samples + 1))
                
                # Count results for each step
                local hmm_count=0
                local cdd_count=0
                local direct_count=0
                local final_count=0
                local classified_count=0
                
                [[ -f "${WORK_DIR}/02.hmmsearch/clean/${sample}.hmm.clean.txt" ]] && \
                    hmm_count=$(wc -l < "${WORK_DIR}/02.hmmsearch/clean/${sample}.hmm.clean.txt" 2>/dev/null || echo 0)
                
                [[ -f "${WORK_DIR}/03.ncbicdd/clean/${sample}.cdd.clean.txt" ]] && \
                    cdd_count=$(wc -l < "${WORK_DIR}/03.ncbicdd/clean/${sample}.cdd.clean.txt" 2>/dev/null || echo 0)
                
                [[ -f "${WORK_DIR}/04.direct/${sample}.wrky_direct.fa" ]] && \
                    direct_count=$(grep -c "^>" "${WORK_DIR}/04.direct/${sample}.wrky_direct.fa" 2>/dev/null || echo 0)
                
                [[ -f "${WORK_DIR}/05.domain_seq/${sample}.final.domain.fa" ]] && \
                    final_count=$(grep -c "^>" "${WORK_DIR}/05.domain_seq/${sample}.final.domain.fa" 2>/dev/null || echo 0)
                
                [[ -f "${WORK_DIR}/06.final_classification/${sample}.final_wrky_summary.txt" ]] && \
                    classified_count=$(grep -v "^#" "${WORK_DIR}/06.final_classification/${sample}.final_wrky_summary.txt" 2>/dev/null | grep -v "^Gene_ID" | wc -l 2>/dev/null || echo 0)
                
                echo "✅ ${sample}: HMM=${hmm_count:-0}, CDD=${cdd_count:-0}, Direct=${direct_count:-0}, Final=${final_count:-0}, Classified=${classified_count:-0}"
            fi
        done
        
        echo ""
        echo "📈 Overall Statistics:"
        echo "---------------------"
        echo "Total samples: ${total_samples:-0}"
        echo "Successful samples: ${successful_samples:-0}"
        echo "Failed samples: ${error_samples:-0}"
        if [[ ${total_samples:-0} -gt 0 ]]; then
            echo "Success rate: $(( (successful_samples * 100) / total_samples ))%"
        else
            echo "Success rate: 0%"
        fi
        
        echo ""
        echo "📊 WRKY Classification Summary:"
        echo "------------------------------"
        
        # 统计所有样品的分类结果
        local total_wrky_genes=0
        local total_group_i=0
        local total_group_ii=0
        local total_group_iii=0
        local total_other=0
        
        for sample in $(get_samples); do
            if ! grep -q "${sample}" "${ERROR_SAMPLE_FILE}" 2>/dev/null; then
                local class_file="${WORK_DIR}/06.final_classification/${sample}.final_wrky_summary.txt"
                if [[ -f "${class_file}" ]]; then
                    local sample_total=$(grep -v "^#" "${class_file}" 2>/dev/null | grep -v "^Gene_ID" | wc -l 2>/dev/null || echo 0)
                    total_wrky_genes=$((total_wrky_genes + ${sample_total:-0}))
                    
                    # 统计各组 - 使用安全的算术运算
                    local group_i=$(grep -v "^#" "${class_file}" 2>/dev/null | grep -v "^Gene_ID" | grep -c "^[^\t]*\tI\t" 2>/dev/null || echo 0)
                    local group_ii=$(grep -v "^#" "${class_file}" 2>/dev/null | grep -v "^Gene_ID" | grep -c "^[^\t]*\tII\t" 2>/dev/null || echo 0)
                    local group_iii=$(grep -v "^#" "${class_file}" 2>/dev/null | grep -v "^Gene_ID" | grep -c "^[^\t]*\tIII\t" 2>/dev/null || echo 0)
                    local other=$(grep -v "^#" "${class_file}" 2>/dev/null | grep -v "^Gene_ID" | grep -c "^[^\t]*\t\(other\|NA\)\t" 2>/dev/null || echo 0)
                    
                    # 使用安全的算术运算
                    total_group_i=$(safe_add "$total_group_i" "$group_i")
                    total_group_ii=$(safe_add "$total_group_ii" "$group_ii")
                    total_group_iii=$(safe_add "$total_group_iii" "$group_iii")
                    total_other=$(safe_add "$total_other" "$other")
                fi
            fi
        done
        
        echo "Total WRKY genes identified: ${total_wrky_genes:-0}"
        echo "Group I: ${total_group_i:-0}"
        echo "Group II: ${total_group_ii:-0}"
        echo "Group III: ${total_group_iii:-0}"
        echo "Other/Unclassified: ${total_other:-0}"
        
        if [[ ${error_samples:-0} -gt 0 ]]; then
            echo ""
            echo "❌ Failed samples:"
            cat "${ERROR_SAMPLE_FILE}" 2>/dev/null || echo "Error reading error file"
        fi
        
    } > "${SUMMARY_FILE}"
    
    log_info "📊 Summary saved to: ${SUMMARY_FILE}"
}

# 🚀 Main execution
main() {
    log_info "🧬 Starting WRKY analysis pipeline (Serial Mode with Classification)..."
    
    # Get all samples
    samples=($(get_samples))
    log_info "📋 Found ${#samples[@]} samples to process"
    log_info "📝 Samples: ${samples[*]}"
    
    # Process samples one by one
    local processed=0
    local successful=0
    local failed=0
    
    for sample in "${samples[@]}"; do
        processed=$((processed + 1))
        log_info ""
        log_info "🚀 Processing sample ${processed}/${#samples[@]}: ${sample}"
        
        # Process sample and handle errors gracefully
        if process_sample "${sample}"; then
            successful=$((successful + 1))
            log_info "✅ Successfully completed ${sample}"
        else
            failed=$((failed + 1))
            log_warn "⚠️ Failed to process ${sample}, continuing with next sample..."
        fi
    done
    
    # Generate summary
    generate_summary
    
    log_info ""
    log_info "🎉 WRKY analysis pipeline completed!"
    log_info "📊 Processing summary: ${processed} total, ${successful} successful, ${failed} failed"
    log_info "📊 Check results in: ${WORK_DIR}"
    log_info "📈 Summary: ${SUMMARY_FILE}"
    log_info "🔬 WRKY classifications: ${WORK_DIR}/06.final_classification/"
    log_info "⏰ Analysis finished at: $(date)"
}

# Check dependencies
check_dependencies() {
    log_info "🔍 Checking dependencies..."
    
    local missing_deps=()
    
    command -v biopytools >/dev/null 2>&1 || missing_deps+=("biopytools")
    command -v hmmsearch >/dev/null 2>&1 || missing_deps+=("hmmsearch")
    command -v rpsblast >/dev/null 2>&1 || missing_deps+=("rpsblast")
    command -v python3 >/dev/null 2>&1 || missing_deps+=("python3")
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        log_error "❌ Missing dependencies: ${missing_deps[*]}"
        exit 1
    fi
    
    # Check Python packages
    if ! python3 -c "import Bio" 2>/dev/null; then
        log_error "❌ Missing Python package: biopython"
        exit 1
    fi
    
    log_info "✅ All dependencies found"
}

# Run the pipeline
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Ensure log directory exists and is writable
    mkdir -p "$(dirname "${LOG_FILE}")"
    
    # Test log file writing
    if ! echo "🧪 Log test" > "${LOG_FILE}" 2>/dev/null; then
        echo "❌ Cannot write to log file: ${LOG_FILE}"
        echo "❌ Current directory: $(pwd)"
        echo "❌ User: $(whoami)"
        echo "❌ Permissions: $(ls -la "$(dirname "${LOG_FILE}")" 2>/dev/null || echo 'Directory not accessible')"
        exit 1
    fi
    
    log_info "🚀 Starting WRKY pipeline with classification in serial mode..."
    
    check_dependencies
    main "$@"
fi
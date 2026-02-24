#!/bin/bash
set -e # If any command fails, exit the script immediately

# --- 1. User Configuration ---
VCF_FILE="variation.filtered.snp.vcf.gz"
WINDOW_SIZE=1000000
OUTPUT_PREFIX="snp_density"

# --- Script Area - No modifications needed below ---

echo "=== Step 1: Create chromosome name map ==="
# This map defines which chromosomes to plot and what their final names should be.
cat << EOF > chr.map
NC_081805.1  Chr1
NC_081806.1  Chr2
NC_081807.1  Chr3
NC_081808.1  Chr4
NC_081809.1  Chr5
NC_081810.1  Chr6
NC_081811.1  Chr7
NC_081812.1  Chr8
NC_081813.1  Chr9
NC_081814.1 Chr10
NC_081815.1 Chr11
NC_081816.1 Chr12
NC_081817.1 Chr13
NC_081818.1 Chr14
NC_081819.1 Chr15
NC_081820.1 Chr16
NC_081821.1 Chr17
NC_081822.1 Chr18
NC_081823.1 Chr19
EOF
echo "File 'chr.map' created successfully."

echo -e "\n=== Step 2: Generate chromosome length file FROM VCF HEADER ==="
# This is a more robust method. It reads your VCF header to get lengths.
# It then uses chr.map to rename the chromosomes and filter for ONLY the ones we want to plot.
gunzip -c ${VCF_FILE} | grep '^##contig' | sed -e 's/.*<ID=//' -e 's/,length=/\t/' -e 's/>.*//' | \
awk 'FNR==NR{map[$1]=$2; next} {if ($1 in map) print map[$1] "\t" $2}' chr.map - > genome.len
echo "File 'genome.len' created successfully. Check its content:"
cat genome.len

echo -e "\n=== Step 3: Extract SNP positions from VCF and convert to BED format ==="
# This command extracts SNP positions and renames chromosomes using chr.map.
# If this step produces an empty file, it means the chromosome names in the VCF *data lines*
# do not match the first column of `chr.map`.
gunzip -c ${VCF_FILE} | grep -v '^#' | awk -v OFS='\t' '{print $1, $2-1, $2}' | \
awk 'FNR==NR{map[$1]=$2; next} {if ($1 in map) {$1=map[$1]; print}}' chr.map - > ${OUTPUT_PREFIX}.positions.bed

# Verification Step: Check if the BED file was created correctly.
if [ ! -s "${OUTPUT_PREFIX}.positions.bed" ]; then
    echo -e "\n\033[0;31m[ERROR]\033[0m The file '${OUTPUT_PREFIX}.positions.bed' is empty."
    echo "This means the chromosome names in your VCF data lines (e.g., 'NC_081805.1') do not match what's in 'chr.map'."
    echo "Please run this command to see the actual names in your data and correct 'chr.map' accordingly:"
    echo "zless ${VCF_FILE} | grep -v '^#' | cut -f1 | sort | uniq"
    exit 1
fi
echo "File '${OUTPUT_PREFIX}.positions.bed' created successfully."

echo -e "\n=== Step 4: Calculate SNP density in specified windows ==="
python -m jcvi.formats.bed bins ${OUTPUT_PREFIX}.positions.bed ${WINDOW_SIZE} > ${OUTPUT_PREFIX}.density.bed
echo "File '${OUTPUT_PREFIX}.density.bed' created successfully."

echo -e "\n=== Step 5: Create JCVI plot layout file ==="
cat << EOF > ${OUTPUT_PREFIX}.layout
# y, x, rotation, height, width
.canvas
8, 5, 0, 1200, 600

# seqid, length (read from the file we generated)
.seqids
$(cat genome.len)

# trackname, color, min, max, height, file
.tracks
heatmap, viridis, 0, db, 80, ${OUTPUT_PREFIX}.density.bed

# options
.config
proportional=no
spacing=20
EOF
echo "File '${OUTPUT_PREFIX}.layout' created successfully."

echo -e "\n=== Step 6: Use JCVI to draw the chromosome density plot ==="
python -m jcvi.graphics.karyotype ${OUTPUT_PREFIX}.layout
echo -e "\n🎉 Plotting complete! Output files are 'karyotype.pdf' and 'karyotype.png'."
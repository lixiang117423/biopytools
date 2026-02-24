#!/bin/bash
# ~/.local/bin 软链接恢复脚本

echo "🔗 开始恢复 ~/.local/bin 软链接..."
echo "=================================="

# 确保目标目录存在
mkdir -p ~/.local/bin

success_count=0
error_count=0

# 恢复软链接: jellyfish
if [ ! -e "$HOME/.local/bin/jellyfish" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/K-mer/bin/jellyfish" "$HOME/.local/bin/jellyfish" 2>/dev/null; then
        echo "✅ 创建软链接: jellyfish"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: jellyfish"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: jellyfish"
fi

# 恢复软链接: hicTransform
if [ ! -e "$HOME/.local/bin/hicTransform" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicTransform" "$HOME/.local/bin/hicTransform" 2>/dev/null; then
        echo "✅ 创建软链接: hicTransform"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicTransform"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicTransform"
fi

# 恢复软链接: xclip
if [ ! -e "$HOME/.local/bin/xclip" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/BioinfTools/bin/xclip" "$HOME/.local/bin/xclip" 2>/dev/null; then
        echo "✅ 创建软链接: xclip"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: xclip"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: xclip"
fi

# 恢复软链接: bowtie2
if [ ! -e "$HOME/.local/bin/bowtie2" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/bowtie2" "$HOME/.local/bin/bowtie2" 2>/dev/null; then
        echo "✅ 创建软链接: bowtie2"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bowtie2"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bowtie2"
fi

# 恢复软链接: R
if [ ! -e "$HOME/.local/bin/R" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/R_v.4.5.1/bin/R" "$HOME/.local/bin/R" 2>/dev/null; then
        echo "✅ 创建软链接: R"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: R"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: R"
fi

# 恢复软链接: nucmer
if [ ! -e "$HOME/.local/bin/nucmer" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mummer_v.4.0.1/bin/nucmer" "$HOME/.local/bin/nucmer" 2>/dev/null; then
        echo "✅ 创建软链接: nucmer"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: nucmer"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: nucmer"
fi

# 恢复软链接: hicCompareMatrices
if [ ! -e "$HOME/.local/bin/hicCompareMatrices" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicCompareMatrices" "$HOME/.local/bin/hicCompareMatrices" 2>/dev/null; then
        echo "✅ 创建软链接: hicCompareMatrices"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicCompareMatrices"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicCompareMatrices"
fi

# 恢复软链接: run_annovar
if [ ! -e "$HOME/.local/bin/run_annovar" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_annovar" "$HOME/.local/bin/run_annovar" 2>/dev/null; then
        echo "✅ 创建软链接: run_annovar"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_annovar"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_annovar"
fi

# 恢复软链接: mutmap
if [ ! -e "$HOME/.local/bin/mutmap" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mutmap/bin/mutmap" "$HOME/.local/bin/mutmap" 2>/dev/null; then
        echo "✅ 创建软链接: mutmap"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mutmap"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mutmap"
fi

# 恢复软链接: assembly-stats
if [ ! -e "$HOME/.local/bin/assembly-stats" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/assembly_stats_v.1.0.1/bin/assembly-stats" "$HOME/.local/bin/assembly-stats" 2>/dev/null; then
        echo "✅ 创建软链接: assembly-stats"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: assembly-stats"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: assembly-stats"
fi

# 恢复软链接: edgeturbo
if [ ! -e "$HOME/.local/bin/edgeturbo" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/edgeturbo/edgeturbo-client/edgeturbo" "$HOME/.local/bin/edgeturbo" 2>/dev/null; then
        echo "✅ 创建软链接: edgeturbo"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: edgeturbo"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: edgeturbo"
fi

# 恢复软链接: run_parabricks
if [ ! -e "$HOME/.local/bin/run_parabricks" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_parabricks" "$HOME/.local/bin/run_parabricks" 2>/dev/null; then
        echo "✅ 创建软链接: run_parabricks"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_parabricks"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_parabricks"
fi

# 恢复软链接: prefetch
if [ ! -e "$HOME/.local/bin/prefetch" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/sratoolkit_v.2.5.7/bin/prefetch" "$HOME/.local/bin/prefetch" 2>/dev/null; then
        echo "✅ 创建软链接: prefetch"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: prefetch"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: prefetch"
fi

# 恢复软链接: ascp
if [ ! -e "$HOME/.local/bin/ascp" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/aspera_v.3.9.6/bin/ascp" "$HOME/.local/bin/ascp" 2>/dev/null; then
        echo "✅ 创建软链接: ascp"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ascp"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ascp"
fi

# 恢复软链接: taxonkit
if [ ! -e "$HOME/.local/bin/taxonkit" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/taxonkit_v.0.20.0/bin/taxonkit" "$HOME/.local/bin/taxonkit" 2>/dev/null; then
        echo "✅ 创建软链接: taxonkit"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: taxonkit"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: taxonkit"
fi

# 恢复软链接: caster-pair
if [ ! -e "$HOME/.local/bin/caster-pair" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/caster_v.1.23/bin/caster-pair" "$HOME/.local/bin/caster-pair" 2>/dev/null; then
        echo "✅ 创建软链接: caster-pair"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: caster-pair"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: caster-pair"
fi

# 恢复软链接: ragtag.py
if [ ! -e "$HOME/.local/bin/ragtag.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag.py" "$HOME/.local/bin/ragtag.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag.py"
fi

# 恢复软链接: fastANI
if [ ! -e "$HOME/.local/bin/fastANI" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/fastani_v.1.34/bin/fastANI" "$HOME/.local/bin/fastANI" 2>/dev/null; then
        echo "✅ 创建软链接: fastANI"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fastANI"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fastANI"
fi

# 恢复软链接: ragtag_asmstats.py
if [ ! -e "$HOME/.local/bin/ragtag_asmstats.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_asmstats.py" "$HOME/.local/bin/ragtag_asmstats.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_asmstats.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_asmstats.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_asmstats.py"
fi

# 恢复软链接: parse_gene_info
if [ ! -e "$HOME/.local/bin/parse_gene_info" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/parse_gene_info" "$HOME/.local/bin/parse_gene_info" 2>/dev/null; then
        echo "✅ 创建软链接: parse_gene_info"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: parse_gene_info"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: parse_gene_info"
fi

# 恢复软链接: hifiasm
if [ ! -e "$HOME/.local/bin/hifiasm" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hifiasm_v.0.25.0/bin/hifiasm" "$HOME/.local/bin/hifiasm" 2>/dev/null; then
        echo "✅ 创建软链接: hifiasm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hifiasm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hifiasm"
fi

# 恢复软链接: multiqc
if [ ! -e "$HOME/.local/bin/multiqc" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/BioinfTools/bin/multiqc" "$HOME/.local/bin/multiqc" 2>/dev/null; then
        echo "✅ 创建软链接: multiqc"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: multiqc"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: multiqc"
fi

# 恢复软链接: hmmscan
if [ ! -e "$HOME/.local/bin/hmmscan" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hmmer_v.3.4/bin/hmmscan" "$HOME/.local/bin/hmmscan" 2>/dev/null; then
        echo "✅ 创建软链接: hmmscan"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hmmscan"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hmmscan"
fi

# 恢复软链接: juicer_tools
if [ ! -e "$HOME/.local/bin/juicer_tools" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/juicer_v.1.6/bin/juicer_tools" "$HOME/.local/bin/juicer_tools" 2>/dev/null; then
        echo "✅ 创建软链接: juicer_tools"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: juicer_tools"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: juicer_tools"
fi

# 恢复软链接: plink
if [ ! -e "$HOME/.local/bin/plink" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Population_genetics/bin/plink" "$HOME/.local/bin/plink" 2>/dev/null; then
        echo "✅ 创建软链接: plink"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: plink"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: plink"
fi

# 恢复软链接: run_kmer_pav
if [ ! -e "$HOME/.local/bin/run_kmer_pav" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_kmer_pav" "$HOME/.local/bin/run_kmer_pav" 2>/dev/null; then
        echo "✅ 创建软链接: run_kmer_pav"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_kmer_pav"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_kmer_pav"
fi

# 恢复软链接: hicBuildMatrix
if [ ! -e "$HOME/.local/bin/hicBuildMatrix" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicBuildMatrix" "$HOME/.local/bin/hicBuildMatrix" 2>/dev/null; then
        echo "✅ 创建软链接: hicBuildMatrix"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicBuildMatrix"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicBuildMatrix"
fi

# 恢复软链接: juicer.sh
if [ ! -e "$HOME/.local/bin/juicer.sh" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/juicer/scripts/juicer.sh" "$HOME/.local/bin/juicer.sh" 2>/dev/null; then
        echo "✅ 创建软链接: juicer.sh"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: juicer.sh"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: juicer.sh"
fi

# 恢复软链接: ragtag_stats.py
if [ ! -e "$HOME/.local/bin/ragtag_stats.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_stats.py" "$HOME/.local/bin/ragtag_stats.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_stats.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_stats.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_stats.py"
fi

# 恢复软链接: RepeatClassifier
if [ ! -e "$HOME/.local/bin/RepeatClassifier" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatClassifier" "$HOME/.local/bin/RepeatClassifier" 2>/dev/null; then
        echo "✅ 创建软链接: RepeatClassifier"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: RepeatClassifier"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: RepeatClassifier"
fi

# 恢复软链接: ragtag_splitasm.py
if [ ! -e "$HOME/.local/bin/ragtag_splitasm.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_splitasm.py" "$HOME/.local/bin/ragtag_splitasm.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_splitasm.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_splitasm.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_splitasm.py"
fi

# 恢复软链接: hicCompartmentalization
if [ ! -e "$HOME/.local/bin/hicCompartmentalization" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicCompartmentalization" "$HOME/.local/bin/hicCompartmentalization" 2>/dev/null; then
        echo "✅ 创建软链接: hicCompartmentalization"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicCompartmentalization"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicCompartmentalization"
fi

# 恢复软链接: orthofinder
if [ ! -e "$HOME/.local/bin/orthofinder" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Orthofinder_v.3.0.1b1/bin/orthofinder" "$HOME/.local/bin/orthofinder" 2>/dev/null; then
        echo "✅ 创建软链接: orthofinder"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: orthofinder"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: orthofinder"
fi

# 恢复软链接: gffread
if [ ! -e "$HOME/.local/bin/gffread" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/gffread" "$HOME/.local/bin/gffread" 2>/dev/null; then
        echo "✅ 创建软链接: gffread"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: gffread"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: gffread"
fi

# 恢复软链接: kmindex
if [ ! -e "$HOME/.local/bin/kmindex" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmindex_v.0.6.0/bin/kmindex" "$HOME/.local/bin/kmindex" 2>/dev/null; then
        echo "✅ 创建软链接: kmindex"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: kmindex"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: kmindex"
fi

# 恢复软链接: kmercountexact.sh
if [ ! -e "$HOME/.local/bin/kmercountexact.sh" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bbtools_v.37.62/bin/kmercountexact.sh" "$HOME/.local/bin/kmercountexact.sh" 2>/dev/null; then
        echo "✅ 创建软链接: kmercountexact.sh"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: kmercountexact.sh"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: kmercountexact.sh"
fi

# 恢复软链接: wgsim
if [ ! -e "$HOME/.local/bin/wgsim" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/GATK_v.4.6.2.0/bin/wgsim" "$HOME/.local/bin/wgsim" 2>/dev/null; then
        echo "✅ 创建软链接: wgsim"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: wgsim"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: wgsim"
fi

# 恢复软链接: meryl
if [ ! -e "$HOME/.local/bin/meryl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/merqury_v.1.3/bin/meryl" "$HOME/.local/bin/meryl" 2>/dev/null; then
        echo "✅ 创建软链接: meryl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: meryl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: meryl"
fi

# 恢复软链接: interproscan-5.jar
if [ ! -e "$HOME/.local/bin/interproscan-5.jar" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan-5.jar" "$HOME/.local/bin/interproscan-5.jar" 2>/dev/null; then
        echo "✅ 创建软链接: interproscan-5.jar"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: interproscan-5.jar"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: interproscan-5.jar"
fi

# 恢复软链接: buildInFont.pl
if [ ! -e "$HOME/.local/bin/buildInFont.pl" ]; then
    if ln -s "buildInFont.pl" "$HOME/.local/bin/buildInFont.pl" 2>/dev/null; then
        echo "✅ 创建软链接: buildInFont.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: buildInFont.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: buildInFont.pl"
fi

# 恢复软链接: coding_change.pl
if [ ! -e "$HOME/.local/bin/coding_change.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/annovar/annovar/coding_change.pl" "$HOME/.local/bin/coding_change.pl" 2>/dev/null; then
        echo "✅ 创建软链接: coding_change.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: coding_change.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: coding_change.pl"
fi

# 恢复软链接: run_methylation_analysis
if [ ! -e "$HOME/.local/bin/run_methylation_analysis" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_methylation_analysis" "$HOME/.local/bin/run_methylation_analysis" 2>/dev/null; then
        echo "✅ 创建软链接: run_methylation_analysis"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_methylation_analysis"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_methylation_analysis"
fi

# 恢复软链接: unikmer
if [ ! -e "$HOME/.local/bin/unikmer" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/unikmer/unikmer" "$HOME/.local/bin/unikmer" 2>/dev/null; then
        echo "✅ 创建软链接: unikmer"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: unikmer"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: unikmer"
fi

# 恢复软链接: blastn
if [ ! -e "$HOME/.local/bin/blastn" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Blast_v.2.16.0/bin/blastn" "$HOME/.local/bin/blastn" 2>/dev/null; then
        echo "✅ 创建软链接: blastn"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: blastn"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: blastn"
fi

# 恢复软链接: KaKs_Calculator
if [ ! -e "$HOME/.local/bin/KaKs_Calculator" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kakscalculator2_v.2.0.1/bin/KaKs_Calculator" "$HOME/.local/bin/KaKs_Calculator" 2>/dev/null; then
        echo "✅ 创建软链接: KaKs_Calculator"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: KaKs_Calculator"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: KaKs_Calculator"
fi

# 恢复软链接: gatk
if [ ! -e "$HOME/.local/bin/gatk" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/GATK_v.4.6.2.0/bin/gatk" "$HOME/.local/bin/gatk" 2>/dev/null; then
        echo "✅ 创建软链接: gatk"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: gatk"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: gatk"
fi

# 恢复软链接: FontSize.pm
if [ ! -e "$HOME/.local/bin/FontSize.pm" ]; then
    if ln -s "FontSize.pm" "$HOME/.local/bin/FontSize.pm" 2>/dev/null; then
        echo "✅ 创建软链接: FontSize.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: FontSize.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: FontSize.pm"
fi

# 恢复软链接: hicNormalize
if [ ! -e "$HOME/.local/bin/hicNormalize" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicNormalize" "$HOME/.local/bin/hicNormalize" 2>/dev/null; then
        echo "✅ 创建软链接: hicNormalize"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicNormalize"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicNormalize"
fi

# 恢复软链接: run_metawrap_pipeline
if [ ! -e "$HOME/.local/bin/run_metawrap_pipeline" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_metawrap_pipeline" "$HOME/.local/bin/run_metawrap_pipeline" 2>/dev/null; then
        echo "✅ 创建软链接: run_metawrap_pipeline"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_metawrap_pipeline"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_metawrap_pipeline"
fi

# 恢复软链接: kmc_tools
if [ ! -e "$HOME/.local/bin/kmc_tools" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmc_v.3.2.4/bin/kmc_tools" "$HOME/.local/bin/kmc_tools" 2>/dev/null; then
        echo "✅ 创建软链接: kmc_tools"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: kmc_tools"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: kmc_tools"
fi

# 恢复软链接: purge_dups
if [ ! -e "$HOME/.local/bin/purge_dups" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/purge_dups_v.1.2.6/bin/purge_dups" "$HOME/.local/bin/purge_dups" 2>/dev/null; then
        echo "✅ 创建软链接: purge_dups"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: purge_dups"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: purge_dups"
fi

# 恢复软链接: hicPlotDistVsCounts
if [ ! -e "$HOME/.local/bin/hicPlotDistVsCounts" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicPlotDistVsCounts" "$HOME/.local/bin/hicPlotDistVsCounts" 2>/dev/null; then
        echo "✅ 创建软链接: hicPlotDistVsCounts"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicPlotDistVsCounts"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicPlotDistVsCounts"
fi

# 恢复软链接: SVG.pm
if [ ! -e "$HOME/.local/bin/SVG.pm" ]; then
    if ln -s "SVG.pm" "$HOME/.local/bin/SVG.pm" 2>/dev/null; then
        echo "✅ 创建软链接: SVG.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: SVG.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: SVG.pm"
fi

# 恢复软链接: checkm
if [ ! -e "$HOME/.local/bin/checkm" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/metaWRAP_v.1.2/bin/checkm" "$HOME/.local/bin/checkm" 2>/dev/null; then
        echo "✅ 创建软链接: checkm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: checkm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: checkm"
fi

# 恢复软链接: RepeatModeler
if [ ! -e "$HOME/.local/bin/RepeatModeler" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatModeler" "$HOME/.local/bin/RepeatModeler" 2>/dev/null; then
        echo "✅ 创建软链接: RepeatModeler"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: RepeatModeler"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: RepeatModeler"
fi

# 恢复软链接: trimal
if [ ! -e "$HOME/.local/bin/trimal" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/trimal_v.1.5.0/bin/trimal" "$HOME/.local/bin/trimal" 2>/dev/null; then
        echo "✅ 创建软链接: trimal"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: trimal"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: trimal"
fi

# 恢复软链接: deduplicate_bismark
if [ ! -e "$HOME/.local/bin/deduplicate_bismark" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/deduplicate_bismark" "$HOME/.local/bin/deduplicate_bismark" 2>/dev/null; then
        echo "✅ 创建软链接: deduplicate_bismark"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: deduplicate_bismark"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: deduplicate_bismark"
fi

# 恢复软链接: run_genome_syn
if [ ! -e "$HOME/.local/bin/run_genome_syn" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_genome_syn" "$HOME/.local/bin/run_genome_syn" 2>/dev/null; then
        echo "✅ 创建软链接: run_genome_syn"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_genome_syn"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_genome_syn"
fi

# 恢复软链接: run_kmer_count
if [ ! -e "$HOME/.local/bin/run_kmer_count" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_kmer_count" "$HOME/.local/bin/run_kmer_count" 2>/dev/null; then
        echo "✅ 创建软链接: run_kmer_count"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_kmer_count"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_kmer_count"
fi

# 恢复软链接: run_mtehylation_pipeline
if [ ! -e "$HOME/.local/bin/run_mtehylation_pipeline" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_mtehylation_pipeline" "$HOME/.local/bin/run_mtehylation_pipeline" 2>/dev/null; then
        echo "✅ 创建软链接: run_mtehylation_pipeline"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_mtehylation_pipeline"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_mtehylation_pipeline"
fi

# 恢复软链接: EDTA_raw.pl
if [ ! -e "$HOME/.local/bin/EDTA_raw.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/EDTA_raw.pl" "$HOME/.local/bin/EDTA_raw.pl" 2>/dev/null; then
        echo "✅ 创建软链接: EDTA_raw.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: EDTA_raw.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: EDTA_raw.pl"
fi

# 恢复软链接: hicSumMatrices
if [ ! -e "$HOME/.local/bin/hicSumMatrices" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicSumMatrices" "$HOME/.local/bin/hicSumMatrices" 2>/dev/null; then
        echo "✅ 创建软链接: hicSumMatrices"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicSumMatrices"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicSumMatrices"
fi

# 恢复软链接: run_genomescope
if [ ! -e "$HOME/.local/bin/run_genomescope" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_genomescope" "$HOME/.local/bin/run_genomescope" 2>/dev/null; then
        echo "✅ 创建软链接: run_genomescope"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_genomescope"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_genomescope"
fi

# 恢复软链接: bowtie2-build
if [ ! -e "$HOME/.local/bin/bowtie2-build" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/bowtie2-build" "$HOME/.local/bin/bowtie2-build" 2>/dev/null; then
        echo "✅ 创建软链接: bowtie2-build"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bowtie2-build"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bowtie2-build"
fi

# 恢复软链接: salmon
if [ ! -e "$HOME/.local/bin/salmon" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/salmon_v.1.10.3/bin/salmon" "$HOME/.local/bin/salmon" 2>/dev/null; then
        echo "✅ 创建软链接: salmon"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: salmon"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: salmon"
fi

# 恢复软链接: hicMergeDomains
if [ ! -e "$HOME/.local/bin/hicMergeDomains" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicMergeDomains" "$HOME/.local/bin/hicMergeDomains" 2>/dev/null; then
        echo "✅ 创建软链接: hicMergeDomains"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicMergeDomains"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicMergeDomains"
fi

# 恢复软链接: plotsr
if [ ! -e "$HOME/.local/bin/plotsr" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Syri_v.1.7.1/bin/plotsr" "$HOME/.local/bin/plotsr" 2>/dev/null; then
        echo "✅ 创建软链接: plotsr"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: plotsr"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: plotsr"
fi

# 恢复软链接: TransDecoder.Predict
if [ ! -e "$HOME/.local/bin/TransDecoder.Predict" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/transdecoder_v.5.5.0/bin/TransDecoder.Predict" "$HOME/.local/bin/TransDecoder.Predict" 2>/dev/null; then
        echo "✅ 创建软链接: TransDecoder.Predict"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: TransDecoder.Predict"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: TransDecoder.Predict"
fi

# 恢复软链接: FastTree
if [ ! -e "$HOME/.local/bin/FastTree" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Orthofinder_v.3.0.1b1/bin/FastTree" "$HOME/.local/bin/FastTree" 2>/dev/null; then
        echo "✅ 创建软链接: FastTree"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: FastTree"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: FastTree"
fi

# 恢复软链接: fasttree
if [ ! -e "$HOME/.local/bin/fasttree" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Orthofinder_v.3.0.1b1/bin/fasttree" "$HOME/.local/bin/fasttree" 2>/dev/null; then
        echo "✅ 创建软链接: fasttree"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fasttree"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fasttree"
fi

# 恢复软链接: hicHyperoptDetectLoops
if [ ! -e "$HOME/.local/bin/hicHyperoptDetectLoops" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicHyperoptDetectLoops" "$HOME/.local/bin/hicHyperoptDetectLoops" 2>/dev/null; then
        echo "✅ 创建软链接: hicHyperoptDetectLoops"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicHyperoptDetectLoops"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicHyperoptDetectLoops"
fi

# 恢复软链接: run_kmer_extractor
if [ ! -e "$HOME/.local/bin/run_kmer_extractor" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_kmer_extractor" "$HOME/.local/bin/run_kmer_extractor" 2>/dev/null; then
        echo "✅ 创建软链接: run_kmer_extractor"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_kmer_extractor"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_kmer_extractor"
fi

# 恢复软链接: filterBam
if [ ! -e "$HOME/.local/bin/filterBam" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bamtools/bin/filterBam" "$HOME/.local/bin/filterBam" 2>/dev/null; then
        echo "✅ 创建软链接: filterBam"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: filterBam"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: filterBam"
fi

# 恢复软链接: hicFindTADs
if [ ! -e "$HOME/.local/bin/hicFindTADs" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicFindTADs" "$HOME/.local/bin/hicFindTADs" 2>/dev/null; then
        echo "✅ 创建软链接: hicFindTADs"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicFindTADs"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicFindTADs"
fi

# 恢复软链接: esummary
if [ ! -e "$HOME/.local/bin/esummary" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/entrez-direct_v.24.0/bin/esummary" "$HOME/.local/bin/esummary" 2>/dev/null; then
        echo "✅ 创建软链接: esummary"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: esummary"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: esummary"
fi

# 恢复软链接: hicDetectLoops
if [ ! -e "$HOME/.local/bin/hicDetectLoops" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicDetectLoops" "$HOME/.local/bin/hicDetectLoops" 2>/dev/null; then
        echo "✅ 创建软链接: hicDetectLoops"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicDetectLoops"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicDetectLoops"
fi

# 恢复软链接: yahs
if [ ! -e "$HOME/.local/bin/yahs" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/yahs_v.1.2.2/bin/yahs" "$HOME/.local/bin/yahs" 2>/dev/null; then
        echo "✅ 创建软链接: yahs"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: yahs"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: yahs"
fi

# 恢复软链接: hicPlotMatrix
if [ ! -e "$HOME/.local/bin/hicPlotMatrix" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicPlotMatrix" "$HOME/.local/bin/hicPlotMatrix" 2>/dev/null; then
        echo "✅ 创建软链接: hicPlotMatrix"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicPlotMatrix"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicPlotMatrix"
fi

# 恢复软链接: parent.pm
if [ ! -e "$HOME/.local/bin/parent.pm" ]; then
    if ln -s "parent.pm" "$HOME/.local/bin/parent.pm" 2>/dev/null; then
        echo "✅ 创建软链接: parent.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: parent.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: parent.pm"
fi

# 恢复软链接: faketime
if [ ! -e "$HOME/.local/bin/faketime" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/faketime/bin/faketime" "$HOME/.local/bin/faketime" 2>/dev/null; then
        echo "✅ 创建软链接: faketime"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: faketime"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: faketime"
fi

# 恢复软链接: Trinity_gene_splice_modeler.py
if [ ! -e "$HOME/.local/bin/Trinity_gene_splice_modeler.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/trinity_v.2.15.2/bin/Trinity_gene_splice_modeler.py" "$HOME/.local/bin/Trinity_gene_splice_modeler.py" 2>/dev/null; then
        echo "✅ 创建软链接: Trinity_gene_splice_modeler.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Trinity_gene_splice_modeler.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Trinity_gene_splice_modeler.py"
fi

# 恢复软链接: TrinityStats.pl
if [ ! -e "$HOME/.local/bin/TrinityStats.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/trinity_v.2.15.2/bin/TrinityStats.pl" "$HOME/.local/bin/TrinityStats.pl" 2>/dev/null; then
        echo "✅ 创建软链接: TrinityStats.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: TrinityStats.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: TrinityStats.pl"
fi

# 恢复软链接: mafft
if [ ! -e "$HOME/.local/bin/mafft" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mafft_v.7.525/bin/mafft" "$HOME/.local/bin/mafft" 2>/dev/null; then
        echo "✅ 创建软链接: mafft"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mafft"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mafft"
fi

# 恢复软链接: ragtag_agp2fa.py
if [ ! -e "$HOME/.local/bin/ragtag_agp2fa.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_agp2fa.py" "$HOME/.local/bin/ragtag_agp2fa.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_agp2fa.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_agp2fa.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_agp2fa.py"
fi

# 恢复软链接: Symmex
if [ ! -e "$HOME/.local/bin/Symmex" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/Symmex" "$HOME/.local/bin/Symmex" 2>/dev/null; then
        echo "✅ 创建软链接: Symmex"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Symmex"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Symmex"
fi

# 恢复软链接: clean_fasta
if [ ! -e "$HOME/.local/bin/clean_fasta" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Python_v.3.13.5/bin/clean_fasta" "$HOME/.local/bin/clean_fasta" 2>/dev/null; then
        echo "✅ 创建软链接: clean_fasta"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: clean_fasta"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: clean_fasta"
fi

# 恢复软链接: tidk
if [ ! -e "$HOME/.local/bin/tidk" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/tidk_v.0.2.65/bin/tidk" "$HOME/.local/bin/tidk" 2>/dev/null; then
        echo "✅ 创建软链接: tidk"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: tidk"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: tidk"
fi

# 恢复软链接: mcmctree
if [ ! -e "$HOME/.local/bin/mcmctree" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/paml_v.4.10.9/bin/mcmctree" "$HOME/.local/bin/mcmctree" 2>/dev/null; then
        echo "✅ 创建软链接: mcmctree"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mcmctree"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mcmctree"
fi

# 恢复软链接: ragtag_agpcheck.py
if [ ! -e "$HOME/.local/bin/ragtag_agpcheck.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_agpcheck.py" "$HOME/.local/bin/ragtag_agpcheck.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_agpcheck.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_agpcheck.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_agpcheck.py"
fi

# 恢复软链接: glnexus_cli
if [ ! -e "$HOME/.local/bin/glnexus_cli" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/glnexus_v.1.4.1/bin/glnexus_cli" "$HOME/.local/bin/glnexus_cli" 2>/dev/null; then
        echo "✅ 创建软链接: glnexus_cli"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: glnexus_cli"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: glnexus_cli"
fi

# 恢复软链接: ragtag_paf2delta.py
if [ ! -e "$HOME/.local/bin/ragtag_paf2delta.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_paf2delta.py" "$HOME/.local/bin/ragtag_paf2delta.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_paf2delta.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_paf2delta.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_paf2delta.py"
fi

# 恢复软链接: fithic
if [ ! -e "$HOME/.local/bin/fithic" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/fithic-v.2.0.8/bin/fithic" "$HOME/.local/bin/fithic" 2>/dev/null; then
        echo "✅ 创建软链接: fithic"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fithic"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fithic"
fi

# 恢复软链接: hicAdjustMatrix
if [ ! -e "$HOME/.local/bin/hicAdjustMatrix" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicAdjustMatrix" "$HOME/.local/bin/hicAdjustMatrix" 2>/dev/null; then
        echo "✅ 创建软链接: hicAdjustMatrix"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicAdjustMatrix"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicAdjustMatrix"
fi

# 恢复软链接: run_repeat_masker
if [ ! -e "$HOME/.local/bin/run_repeat_masker" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_repeat_masker" "$HOME/.local/bin/run_repeat_masker" 2>/dev/null; then
        echo "✅ 创建软链接: run_repeat_masker"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_repeat_masker"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_repeat_masker"
fi

# 恢复软链接: LTR_HARVEST_parallel
if [ ! -e "$HOME/.local/bin/LTR_HARVEST_parallel" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/ltr_harvest_parallel_v.1.2/bin/LTR_HARVEST_parallel" "$HOME/.local/bin/LTR_HARVEST_parallel" 2>/dev/null; then
        echo "✅ 创建软链接: LTR_HARVEST_parallel"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: LTR_HARVEST_parallel"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: LTR_HARVEST_parallel"
fi

# 恢复软链接: esearch
if [ ! -e "$HOME/.local/bin/esearch" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/entrez-direct_v.24.0/bin/esearch" "$HOME/.local/bin/esearch" 2>/dev/null; then
        echo "✅ 创建软链接: esearch"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: esearch"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: esearch"
fi

# 恢复软链接: minimap2
if [ ! -e "$HOME/.local/bin/minimap2" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Genome_dedup/bin/minimap2" "$HOME/.local/bin/minimap2" 2>/dev/null; then
        echo "✅ 创建软链接: minimap2"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: minimap2"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: minimap2"
fi

# 恢复软链接: qtlplot
if [ ! -e "$HOME/.local/bin/qtlplot" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/qtlseq/bin/qtlplot" "$HOME/.local/bin/qtlplot" 2>/dev/null; then
        echo "✅ 创建软链接: qtlplot"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: qtlplot"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: qtlplot"
fi

# 恢复软链接: run_kmer_analysis
if [ ! -e "$HOME/.local/bin/run_kmer_analysis" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_kmer_analysis" "$HOME/.local/bin/run_kmer_analysis" 2>/dev/null; then
        echo "✅ 创建软链接: run_kmer_analysis"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_kmer_analysis"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_kmer_analysis"
fi

# 恢复软链接: vcftools
if [ ! -e "$HOME/.local/bin/vcftools" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Population_genetics/bin/vcftools" "$HOME/.local/bin/vcftools" 2>/dev/null; then
        echo "✅ 创建软链接: vcftools"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: vcftools"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: vcftools"
fi

# 恢复软链接: run_vcf_pca
if [ ! -e "$HOME/.local/bin/run_vcf_pca" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_vcf_pca" "$HOME/.local/bin/run_vcf_pca" 2>/dev/null; then
        echo "✅ 创建软链接: run_vcf_pca"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_vcf_pca"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_vcf_pca"
fi

# 恢复软链接: RepeatMasker
if [ ! -e "$HOME/.local/bin/RepeatMasker" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/repeat_identiy/bin/RepeatMasker" "$HOME/.local/bin/RepeatMasker" 2>/dev/null; then
        echo "✅ 创建软链接: RepeatMasker"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: RepeatMasker"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: RepeatMasker"
fi

# 恢复软链接: run_gtx
if [ ! -e "$HOME/.local/bin/run_gtx" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_gtx" "$HOME/.local/bin/run_gtx" 2>/dev/null; then
        echo "✅ 创建软链接: run_gtx"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_gtx"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_gtx"
fi

# 恢复软链接: show-coords
if [ ! -e "$HOME/.local/bin/show-coords" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mummer_v.3.23/bin/show-coords" "$HOME/.local/bin/show-coords" 2>/dev/null; then
        echo "✅ 创建软链接: show-coords"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: show-coords"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: show-coords"
fi

# 恢复软链接: bismark
if [ ! -e "$HOME/.local/bin/bismark" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/bismark" "$HOME/.local/bin/bismark" 2>/dev/null; then
        echo "✅ 创建软链接: bismark"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bismark"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bismark"
fi

# 恢复软链接: trf
if [ ! -e "$HOME/.local/bin/trf" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/trf" "$HOME/.local/bin/trf" 2>/dev/null; then
        echo "✅ 创建软链接: trf"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: trf"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: trf"
fi

# 恢复软链接: hicCorrectMatrix
if [ ! -e "$HOME/.local/bin/hicCorrectMatrix" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicCorrectMatrix" "$HOME/.local/bin/hicCorrectMatrix" 2>/dev/null; then
        echo "✅ 创建软链接: hicCorrectMatrix"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicCorrectMatrix"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicCorrectMatrix"
fi

# 恢复软链接: variants_reduction.pl
if [ ! -e "$HOME/.local/bin/variants_reduction.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/annovar/annovar/variants_reduction.pl" "$HOME/.local/bin/variants_reduction.pl" 2>/dev/null; then
        echo "✅ 创建软链接: variants_reduction.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: variants_reduction.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: variants_reduction.pl"
fi

# 恢复软链接: checkm2
if [ ! -e "$HOME/.local/bin/checkm2" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/checkm_v.1.1.0/bin/checkm2" "$HOME/.local/bin/checkm2" 2>/dev/null; then
        echo "✅ 创建软链接: checkm2"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: checkm2"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: checkm2"
fi

# 恢复软链接: hicConvertFormat
if [ ! -e "$HOME/.local/bin/hicConvertFormat" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicConvertFormat" "$HOME/.local/bin/hicConvertFormat" 2>/dev/null; then
        echo "✅ 创建软链接: hicConvertFormat"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicConvertFormat"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicConvertFormat"
fi

# 恢复软链接: jcvi
if [ ! -e "$HOME/.local/bin/jcvi" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/JCVI_v.1.5.6/bin/jcvi" "$HOME/.local/bin/jcvi" 2>/dev/null; then
        echo "✅ 创建软链接: jcvi"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: jcvi"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: jcvi"
fi

# 恢复软链接: makeprofiledb
if [ ! -e "$HOME/.local/bin/makeprofiledb" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Blast_v.2.16.0/bin/makeprofiledb" "$HOME/.local/bin/makeprofiledb" 2>/dev/null; then
        echo "✅ 创建软链接: makeprofiledb"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: makeprofiledb"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: makeprofiledb"
fi

# 恢复软链接: mmseqs
if [ ! -e "$HOME/.local/bin/mmseqs" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mmseqs2_v.16.747c6/bin/mmseqs" "$HOME/.local/bin/mmseqs" 2>/dev/null; then
        echo "✅ 创建软链接: mmseqs"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mmseqs"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mmseqs"
fi

# 恢复软链接: run_transcriptome_prediction
if [ ! -e "$HOME/.local/bin/run_transcriptome_prediction" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_transcriptome_prediction" "$HOME/.local/bin/run_transcriptome_prediction" 2>/dev/null; then
        echo "✅ 创建软链接: run_transcriptome_prediction"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_transcriptome_prediction"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_transcriptome_prediction"
fi

# 恢复软链接: caster-site
if [ ! -e "$HOME/.local/bin/caster-site" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/caster_v.1.23/bin/caster-site" "$HOME/.local/bin/caster-site" 2>/dev/null; then
        echo "✅ 创建软链接: caster-site"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: caster-site"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: caster-site"
fi

# 恢复软链接: sniffles
if [ ! -e "$HOME/.local/bin/sniffles" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/sniffles_v.2.6.3/bin/sniffles" "$HOME/.local/bin/sniffles" 2>/dev/null; then
        echo "✅ 创建软链接: sniffles"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: sniffles"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: sniffles"
fi

# 恢复软链接: extract_splice_sites.py
if [ ! -e "$HOME/.local/bin/extract_splice_sites.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/extract_splice_sites.py" "$HOME/.local/bin/extract_splice_sites.py" 2>/dev/null; then
        echo "✅ 创建软链接: extract_splice_sites.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: extract_splice_sites.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: extract_splice_sites.py"
fi

# 恢复软链接: convert2annovar.pl
if [ ! -e "$HOME/.local/bin/convert2annovar.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/annovar/annovar/convert2annovar.pl" "$HOME/.local/bin/convert2annovar.pl" 2>/dev/null; then
        echo "✅ 创建软链接: convert2annovar.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: convert2annovar.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: convert2annovar.pl"
fi

# 恢复软链接: python
if [ ! -e "$HOME/.local/bin/python" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/python3" "$HOME/.local/bin/python" 2>/dev/null; then
        echo "✅ 创建软链接: python"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: python"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: python"
fi

# 恢复软链接: ragtag_patch.py
if [ ! -e "$HOME/.local/bin/ragtag_patch.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_patch.py" "$HOME/.local/bin/ragtag_patch.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_patch.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_patch.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_patch.py"
fi

# 恢复软链接: hisat2-build
if [ ! -e "$HOME/.local/bin/hisat2-build" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/hisat2-build" "$HOME/.local/bin/hisat2-build" 2>/dev/null; then
        echo "✅ 创建软链接: hisat2-build"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hisat2-build"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hisat2-build"
fi

# 恢复软链接: extract_exons.py
if [ ! -e "$HOME/.local/bin/extract_exons.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/extract_exons.py" "$HOME/.local/bin/extract_exons.py" 2>/dev/null; then
        echo "✅ 创建软链接: extract_exons.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: extract_exons.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: extract_exons.py"
fi

# 恢复软链接: hicTrainTADClassifier
if [ ! -e "$HOME/.local/bin/hicTrainTADClassifier" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicTrainTADClassifier" "$HOME/.local/bin/hicTrainTADClassifier" 2>/dev/null; then
        echo "✅ 创建软链接: hicTrainTADClassifier"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicTrainTADClassifier"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicTrainTADClassifier"
fi

# 恢复软链接: iqtree
if [ ! -e "$HOME/.local/bin/iqtree" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/iqtree_v.3.0.1/bin/iqtree" "$HOME/.local/bin/iqtree" 2>/dev/null; then
        echo "✅ 创建软链接: iqtree"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: iqtree"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: iqtree"
fi

# 恢复软链接: hisat2
if [ ! -e "$HOME/.local/bin/hisat2" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/hisat2" "$HOME/.local/bin/hisat2" 2>/dev/null; then
        echo "✅ 创建软链接: hisat2"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hisat2"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hisat2"
fi

# 恢复软链接: merqury.sh
if [ ! -e "$HOME/.local/bin/merqury.sh" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/merqury_v.1.3/bin/merqury.sh" "$HOME/.local/bin/merqury.sh" 2>/dev/null; then
        echo "✅ 创建软链接: merqury.sh"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: merqury.sh"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: merqury.sh"
fi

# 恢复软链接: MCScanX2Link.pl
if [ ! -e "$HOME/.local/bin/MCScanX2Link.pl" ]; then
    if ln -s "MCScanX2Link.pl" "$HOME/.local/bin/MCScanX2Link.pl" 2>/dev/null; then
        echo "✅ 创建软链接: MCScanX2Link.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: MCScanX2Link.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: MCScanX2Link.pl"
fi

# 恢复软链接: nhmmscan
if [ ! -e "$HOME/.local/bin/nhmmscan" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hmmer_v.3.4/bin/nhmmscan" "$HOME/.local/bin/nhmmscan" 2>/dev/null; then
        echo "✅ 创建软链接: nhmmscan"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: nhmmscan"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: nhmmscan"
fi

# 恢复软链接: NGenomeSyn
if [ ! -e "$HOME/.local/bin/NGenomeSyn" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/NGenomeSyn/NGenomeSyn-1.43/bin/NGenomeSyn" "$HOME/.local/bin/NGenomeSyn" 2>/dev/null; then
        echo "✅ 创建软链接: NGenomeSyn"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: NGenomeSyn"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: NGenomeSyn"
fi

# 恢复软链接: kraken2
if [ ! -e "$HOME/.local/bin/kraken2" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kraken_v.2.17/bin/kraken2" "$HOME/.local/bin/kraken2" 2>/dev/null; then
        echo "✅ 创建软链接: kraken2"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: kraken2"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: kraken2"
fi

# 恢复软链接: BuildDatabase
if [ ! -e "$HOME/.local/bin/BuildDatabase" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/BuildDatabase" "$HOME/.local/bin/BuildDatabase" 2>/dev/null; then
        echo "✅ 创建软链接: BuildDatabase"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: BuildDatabase"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: BuildDatabase"
fi

# 恢复软链接: bedtools
if [ ! -e "$HOME/.local/bin/bedtools" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Population_genetics/bin/bedtools" "$HOME/.local/bin/bedtools" 2>/dev/null; then
        echo "✅ 创建软链接: bedtools"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bedtools"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bedtools"
fi

# 恢复软链接: blastx
if [ ! -e "$HOME/.local/bin/blastx" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Blast_v.2.16.0/bin/blastx" "$HOME/.local/bin/blastx" 2>/dev/null; then
        echo "✅ 创建软链接: blastx"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: blastx"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: blastx"
fi

# 恢复软链接: signalp6
if [ ! -e "$HOME/.local/bin/signalp6" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/signalp6/bin/signalp6" "$HOME/.local/bin/signalp6" 2>/dev/null; then
        echo "✅ 创建软链接: signalp6"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: signalp6"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: signalp6"
fi

# 恢复软链接: hicPlotAverageRegions
if [ ! -e "$HOME/.local/bin/hicPlotAverageRegions" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicPlotAverageRegions" "$HOME/.local/bin/hicPlotAverageRegions" 2>/dev/null; then
        echo "✅ 创建软链接: hicPlotAverageRegions"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicPlotAverageRegions"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicPlotAverageRegions"
fi

# 恢复软链接: Element.pm
if [ ! -e "$HOME/.local/bin/Element.pm" ]; then
    if ln -s "Element.pm" "$HOME/.local/bin/Element.pm" 2>/dev/null; then
        echo "✅ 创建软链接: Element.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Element.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Element.pm"
fi

# 恢复软链接: run_rnaseq
if [ ! -e "$HOME/.local/bin/run_rnaseq" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_rnaseq" "$HOME/.local/bin/run_rnaseq" 2>/dev/null; then
        echo "✅ 创建软链接: run_rnaseq"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_rnaseq"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_rnaseq"
fi

# 恢复软链接: racon
if [ ! -e "$HOME/.local/bin/racon" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/racon_v.1.5.0/bin/racon" "$HOME/.local/bin/racon" 2>/dev/null; then
        echo "✅ 创建软链接: racon"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: racon"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: racon"
fi

# 恢复软链接: run_bam_depth
if [ ! -e "$HOME/.local/bin/run_bam_depth" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_bam_depth" "$HOME/.local/bin/run_bam_depth" 2>/dev/null; then
        echo "✅ 创建软链接: run_bam_depth"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_bam_depth"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_bam_depth"
fi

# 恢复软链接: LTR_retriever
if [ ! -e "$HOME/.local/bin/LTR_retriever" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/ltr_retriever_v.3.0.1/bin/LTR_retriever" "$HOME/.local/bin/LTR_retriever" 2>/dev/null; then
        echo "✅ 创建软链接: LTR_retriever"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: LTR_retriever"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: LTR_retriever"
fi

# 恢复软链接: XML.pm
if [ ! -e "$HOME/.local/bin/XML.pm" ]; then
    if ln -s "XML.pm" "$HOME/.local/bin/XML.pm" 2>/dev/null; then
        echo "✅ 创建软链接: XML.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: XML.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: XML.pm"
fi

# 恢复软链接: DOM.pm
if [ ! -e "$HOME/.local/bin/DOM.pm" ]; then
    if ln -s "DOM.pm" "$HOME/.local/bin/DOM.pm" 2>/dev/null; then
        echo "✅ 创建软链接: DOM.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: DOM.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: DOM.pm"
fi

# 恢复软链接: xtract
if [ ! -e "$HOME/.local/bin/xtract" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/entrez-direct_v.24.0/bin/xtract" "$HOME/.local/bin/xtract" 2>/dev/null; then
        echo "✅ 创建软链接: xtract"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: xtract"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: xtract"
fi

# 恢复软链接: run_admixture
if [ ! -e "$HOME/.local/bin/run_admixture" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_admixture" "$HOME/.local/bin/run_admixture" 2>/dev/null; then
        echo "✅ 创建软链接: run_admixture"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_admixture"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_admixture"
fi

# 恢复软链接: clustalo
if [ ! -e "$HOME/.local/bin/clustalo" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/clustalo_v.1.2.4/bin/clustalo" "$HOME/.local/bin/clustalo" 2>/dev/null; then
        echo "✅ 创建软链接: clustalo"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: clustalo"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: clustalo"
fi

# 恢复软链接: bismark_genome_preparation
if [ ! -e "$HOME/.local/bin/bismark_genome_preparation" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/bismark_genome_preparation" "$HOME/.local/bin/bismark_genome_preparation" 2>/dev/null; then
        echo "✅ 创建软链接: bismark_genome_preparation"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bismark_genome_preparation"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bismark_genome_preparation"
fi

# 恢复软链接: run_vcf_extractor
if [ ! -e "$HOME/.local/bin/run_vcf_extractor" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_vcf_extractor" "$HOME/.local/bin/run_vcf_extractor" 2>/dev/null; then
        echo "✅ 创建软链接: run_vcf_extractor"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_vcf_extractor"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_vcf_extractor"
fi

# 恢复软链接: seqtk
if [ ! -e "$HOME/.local/bin/seqtk" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/seqtk/bin/seqtk" "$HOME/.local/bin/seqtk" 2>/dev/null; then
        echo "✅ 创建软链接: seqtk"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: seqtk"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: seqtk"
fi

# 恢复软链接: run_gtx_joint
if [ ! -e "$HOME/.local/bin/run_gtx_joint" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_gtx_joint" "$HOME/.local/bin/run_gtx_joint" 2>/dev/null; then
        echo "✅ 创建软链接: run_gtx_joint"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_gtx_joint"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_gtx_joint"
fi

# 恢复软链接: bam2fastq
if [ ! -e "$HOME/.local/bin/bam2fastq" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/pbbam_v.2.4.0/bin/bam2fastq" "$HOME/.local/bin/bam2fastq" 2>/dev/null; then
        echo "✅ 创建软链接: bam2fastq"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bam2fastq"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bam2fastq"
fi

# 恢复软链接: bwa
if [ ! -e "$HOME/.local/bin/bwa" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Population_genetics/bin/bwa" "$HOME/.local/bin/bwa" 2>/dev/null; then
        echo "✅ 创建软链接: bwa"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bwa"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bwa"
fi

# 恢复软链接: gtdbtk
if [ ! -e "$HOME/.local/bin/gtdbtk" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/GTDB-Tk_v.2.4.1/bin/gtdbtk" "$HOME/.local/bin/gtdbtk" 2>/dev/null; then
        echo "✅ 创建软链接: gtdbtk"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: gtdbtk"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: gtdbtk"
fi

# 恢复软链接: RepeatProteinMask
if [ ! -e "$HOME/.local/bin/RepeatProteinMask" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/repeat_identiy/bin/RepeatProteinMask" "$HOME/.local/bin/RepeatProteinMask" 2>/dev/null; then
        echo "✅ 创建软链接: RepeatProteinMask"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: RepeatProteinMask"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: RepeatProteinMask"
fi

# 恢复软链接: gemma
if [ ! -e "$HOME/.local/bin/gemma" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/gemma_v.0.98.5/bin/gemma" "$HOME/.local/bin/gemma" 2>/dev/null; then
        echo "✅ 创建软链接: gemma"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: gemma"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: gemma"
fi

# 恢复软链接: TransDecoder.LongOrfs
if [ ! -e "$HOME/.local/bin/TransDecoder.LongOrfs" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/transdecoder_v.5.5.0/bin/TransDecoder.LongOrfs" "$HOME/.local/bin/TransDecoder.LongOrfs" 2>/dev/null; then
        echo "✅ 创建软链接: TransDecoder.LongOrfs"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: TransDecoder.LongOrfs"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: TransDecoder.LongOrfs"
fi

# 恢复软链接: pasa
if [ ! -e "$HOME/.local/bin/pasa" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/pasa_v.2.5.3/bin/pasa" "$HOME/.local/bin/pasa" 2>/dev/null; then
        echo "✅ 创建软链接: pasa"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: pasa"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: pasa"
fi

# 恢复软链接: hicAggregateContacts
if [ ! -e "$HOME/.local/bin/hicAggregateContacts" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicAggregateContacts" "$HOME/.local/bin/hicAggregateContacts" 2>/dev/null; then
        echo "✅ 创建软链接: hicAggregateContacts"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicAggregateContacts"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicAggregateContacts"
fi

# 恢复软链接: busco
if [ ! -e "$HOME/.local/bin/busco" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/BUSCO_v.6.0.0/bin/busco" "$HOME/.local/bin/busco" 2>/dev/null; then
        echo "✅ 创建软链接: busco"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: busco"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: busco"
fi

# 恢复软链接: run_hifiasm
if [ ! -e "$HOME/.local/bin/run_hifiasm" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_hifiasm" "$HOME/.local/bin/run_hifiasm" 2>/dev/null; then
        echo "✅ 创建软链接: run_hifiasm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_hifiasm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_hifiasm"
fi

# 恢复软链接: LDBlockShow
if [ ! -e "$HOME/.local/bin/LDBlockShow" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/LDBlockShow/LDBlockShow-1.41/bin/LDBlockShow" "$HOME/.local/bin/LDBlockShow" 2>/dev/null; then
        echo "✅ 创建软链接: LDBlockShow"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: LDBlockShow"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: LDBlockShow"
fi

# 恢复软链接: ragtag_rename.py
if [ ! -e "$HOME/.local/bin/ragtag_rename.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_rename.py" "$HOME/.local/bin/ragtag_rename.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_rename.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_rename.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_rename.py"
fi

# 恢复软链接: ragtag_update_gff.py
if [ ! -e "$HOME/.local/bin/ragtag_update_gff.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_update_gff.py" "$HOME/.local/bin/ragtag_update_gff.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_update_gff.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_update_gff.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_update_gff.py"
fi

# 恢复软链接: aria2c
if [ ! -e "$HOME/.local/bin/aria2c" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/aria2/bin/aria2c" "$HOME/.local/bin/aria2c" 2>/dev/null; then
        echo "✅ 创建软链接: aria2c"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: aria2c"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: aria2c"
fi

# 恢复软链接: sra-stat
if [ ! -e "$HOME/.local/bin/sra-stat" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/sratoolkit_v.2.5.7/bin/sra-stat" "$HOME/.local/bin/sra-stat" 2>/dev/null; then
        echo "✅ 创建软链接: sra-stat"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: sra-stat"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: sra-stat"
fi

# 恢复软链接: HelitronScanner
if [ ! -e "$HOME/.local/bin/HelitronScanner" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/HelitronScanner" "$HOME/.local/bin/HelitronScanner" 2>/dev/null; then
        echo "✅ 创建软链接: HelitronScanner"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: HelitronScanner"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: HelitronScanner"
fi

# 恢复软链接: 3d-dna
if [ ! -e "$HOME/.local/bin/3d-dna" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/3d-dna_v.201008/bin/3d-dna" "$HOME/.local/bin/3d-dna" 2>/dev/null; then
        echo "✅ 创建软链接: 3d-dna"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: 3d-dna"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: 3d-dna"
fi

# 恢复软链接: run_vcf_njtree
if [ ! -e "$HOME/.local/bin/run_vcf_njtree" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_vcf_njtree" "$HOME/.local/bin/run_vcf_njtree" 2>/dev/null; then
        echo "✅ 创建软链接: run_vcf_njtree"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_vcf_njtree"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_vcf_njtree"
fi

# 恢复软链接: filter_bam
if [ ! -e "$HOME/.local/bin/filter_bam" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/filter_bam" "$HOME/.local/bin/filter_bam" 2>/dev/null; then
        echo "✅ 创建软链接: filter_bam"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: filter_bam"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: filter_bam"
fi

# 恢复软链接: hicCorrelate
if [ ! -e "$HOME/.local/bin/hicCorrelate" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicCorrelate" "$HOME/.local/bin/hicCorrelate" 2>/dev/null; then
        echo "✅ 创建软链接: hicCorrelate"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicCorrelate"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicCorrelate"
fi

# 恢复软链接: blastp
if [ ! -e "$HOME/.local/bin/blastp" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Blast_v.2.16.0/bin/blastp" "$HOME/.local/bin/blastp" 2>/dev/null; then
        echo "✅ 创建软链接: blastp"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: blastp"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: blastp"
fi

# 恢复软链接: Pipeliner.pm
if [ ! -e "$HOME/.local/bin/Pipeliner.pm" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/evidencemodeler/bin/PerlLib/Pipeliner.pm" "$HOME/.local/bin/Pipeliner.pm" 2>/dev/null; then
        echo "✅ 创建软链接: Pipeliner.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Pipeliner.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Pipeliner.pm"
fi

# 恢复软链接: run_vcf_filter
if [ ! -e "$HOME/.local/bin/run_vcf_filter" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_vcf_filter" "$HOME/.local/bin/run_vcf_filter" 2>/dev/null; then
        echo "✅ 创建软链接: run_vcf_filter"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_vcf_filter"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_vcf_filter"
fi

# 恢复软链接: hicDifferentialTAD
if [ ! -e "$HOME/.local/bin/hicDifferentialTAD" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicDifferentialTAD" "$HOME/.local/bin/hicDifferentialTAD" 2>/dev/null; then
        echo "✅ 创建软链接: hicDifferentialTAD"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicDifferentialTAD"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicDifferentialTAD"
fi

# 恢复软链接: bcftools
if [ ! -e "$HOME/.local/bin/bcftools" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bcftools_v.1.22/bin/bcftools" "$HOME/.local/bin/bcftools" 2>/dev/null; then
        echo "✅ 创建软链接: bcftools"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bcftools"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bcftools"
fi

# 恢复软链接: mumemto
if [ ! -e "$HOME/.local/bin/mumemto" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Mumemto_v.1.3.0/bin/mumemto" "$HOME/.local/bin/mumemto" 2>/dev/null; then
        echo "✅ 创建软链接: mumemto"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mumemto"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mumemto"
fi

# 恢复软链接: hicQuickQC
if [ ! -e "$HOME/.local/bin/hicQuickQC" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicQuickQC" "$HOME/.local/bin/hicQuickQC" 2>/dev/null; then
        echo "✅ 创建软链接: hicQuickQC"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicQuickQC"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicQuickQC"
fi

# 恢复软链接: run_blast
if [ ! -e "$HOME/.local/bin/run_blast" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_blast" "$HOME/.local/bin/run_blast" 2>/dev/null; then
        echo "✅ 创建软链接: run_blast"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_blast"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_blast"
fi

# 恢复软链接: generate_font_embed_defs.pl
if [ ! -e "$HOME/.local/bin/generate_font_embed_defs.pl" ]; then
    if ln -s "generate_font_embed_defs.pl" "$HOME/.local/bin/generate_font_embed_defs.pl" 2>/dev/null; then
        echo "✅ 创建软链接: generate_font_embed_defs.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: generate_font_embed_defs.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: generate_font_embed_defs.pl"
fi

# 恢复软链接: ragtag_correct.py
if [ ! -e "$HOME/.local/bin/ragtag_correct.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_correct.py" "$HOME/.local/bin/ragtag_correct.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_correct.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_correct.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_correct.py"
fi

# 恢复软链接: hicHyperoptDetectLoopsHiCCUPS
if [ ! -e "$HOME/.local/bin/hicHyperoptDetectLoopsHiCCUPS" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicHyperoptDetectLoopsHiCCUPS" "$HOME/.local/bin/hicHyperoptDetectLoopsHiCCUPS" 2>/dev/null; then
        echo "✅ 创建软链接: hicHyperoptDetectLoopsHiCCUPS"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicHyperoptDetectLoopsHiCCUPS"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicHyperoptDetectLoopsHiCCUPS"
fi

# 恢复软链接: tabix
if [ ! -e "$HOME/.local/bin/tabix" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/GATK_v.4.6.2.0/bin/tabix" "$HOME/.local/bin/tabix" 2>/dev/null; then
        echo "✅ 创建软链接: tabix"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: tabix"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: tabix"
fi

# 恢复软链接: Switch.pm
if [ ! -e "$HOME/.local/bin/Switch.pm" ]; then
    if ln -s "Switch.pm" "$HOME/.local/bin/Switch.pm" 2>/dev/null; then
        echo "✅ 创建软链接: Switch.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Switch.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Switch.pm"
fi

# 恢复软链接: TEsorter
if [ ! -e "$HOME/.local/bin/TEsorter" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/tesorter_v.1.4.7/bin/TEsorter" "$HOME/.local/bin/TEsorter" 2>/dev/null; then
        echo "✅ 创建软链接: TEsorter"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: TEsorter"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: TEsorter"
fi

# 恢复软链接: hicPCA
if [ ! -e "$HOME/.local/bin/hicPCA" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicPCA" "$HOME/.local/bin/hicPCA" 2>/dev/null; then
        echo "✅ 创建软链接: hicPCA"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicPCA"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicPCA"
fi

# 恢复软链接: smc++
if [ ! -e "$HOME/.local/bin/smc++" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/smcpp_v.1.15.2/bin/smc++" "$HOME/.local/bin/smc++" 2>/dev/null; then
        echo "✅ 创建软链接: smc++"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: smc++"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: smc++"
fi

# 恢复软链接: epost
if [ ! -e "$HOME/.local/bin/epost" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/entrez-direct_v.24.0/bin/epost" "$HOME/.local/bin/epost" 2>/dev/null; then
        echo "✅ 创建软链接: epost"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: epost"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: epost"
fi

# 恢复软链接: interproscan.properties
if [ ! -e "$HOME/.local/bin/interproscan.properties" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.properties" "$HOME/.local/bin/interproscan.properties" 2>/dev/null; then
        echo "✅ 创建软链接: interproscan.properties"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: interproscan.properties"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: interproscan.properties"
fi

# 恢复软链接: kmc_dump
if [ ! -e "$HOME/.local/bin/kmc_dump" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmc_v.3.2.4/bin/kmc_dump" "$HOME/.local/bin/kmc_dump" 2>/dev/null; then
        echo "✅ 创建软链接: kmc_dump"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: kmc_dump"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: kmc_dump"
fi

# 恢复软链接: run_bam_stats
if [ ! -e "$HOME/.local/bin/run_bam_stats" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_bam_stats" "$HOME/.local/bin/run_bam_stats" 2>/dev/null; then
        echo "✅ 创建软链接: run_bam_stats"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_bam_stats"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_bam_stats"
fi

# 恢复软链接: FastK
if [ ! -e "$HOME/.local/bin/FastK" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/FastK" "$HOME/.local/bin/FastK" 2>/dev/null; then
        echo "✅ 创建软链接: FastK"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: FastK"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: FastK"
fi

# 恢复软链接: augustus
if [ ! -e "$HOME/.local/bin/augustus" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Augustus_v.3.5.0/bin/augustus" "$HOME/.local/bin/augustus" 2>/dev/null; then
        echo "✅ 创建软链接: augustus"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: augustus"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: augustus"
fi

# 恢复软链接: samtools
if [ ! -e "$HOME/.local/bin/samtools" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools" "$HOME/.local/bin/samtools" 2>/dev/null; then
        echo "✅ 创建软链接: samtools"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: samtools"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: samtools"
fi

# 恢复软链接: Extension.pm
if [ ! -e "$HOME/.local/bin/Extension.pm" ]; then
    if ln -s "Extension.pm" "$HOME/.local/bin/Extension.pm" 2>/dev/null; then
        echo "✅ 创建软链接: Extension.pm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Extension.pm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Extension.pm"
fi

# 恢复软链接: syri
if [ ! -e "$HOME/.local/bin/syri" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Syri_v.1.7.1/bin/syri" "$HOME/.local/bin/syri" 2>/dev/null; then
        echo "✅ 创建软链接: syri"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: syri"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: syri"
fi

# 恢复软链接: mdust
if [ ! -e "$HOME/.local/bin/mdust" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/mdust" "$HOME/.local/bin/mdust" 2>/dev/null; then
        echo "✅ 创建软链接: mdust"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mdust"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mdust"
fi

# 恢复软链接: kmc
if [ ! -e "$HOME/.local/bin/kmc" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmc_v.3.2.4/bin/kmc" "$HOME/.local/bin/kmc" 2>/dev/null; then
        echo "✅ 创建软链接: kmc"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: kmc"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: kmc"
fi

# 恢复软链接: stringtie
if [ ! -e "$HOME/.local/bin/stringtie" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/stringtie" "$HOME/.local/bin/stringtie" 2>/dev/null; then
        echo "✅ 创建软链接: stringtie"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: stringtie"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: stringtie"
fi

# 恢复软链接: annotate_variation.pl
if [ ! -e "$HOME/.local/bin/annotate_variation.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/annovar/annovar/annotate_variation.pl" "$HOME/.local/bin/annotate_variation.pl" 2>/dev/null; then
        echo "✅ 创建软链接: annotate_variation.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: annotate_variation.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: annotate_variation.pl"
fi

# 恢复软链接: efetch
if [ ! -e "$HOME/.local/bin/efetch" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/entrez-direct_v.24.0/bin/efetch" "$HOME/.local/bin/efetch" 2>/dev/null; then
        echo "✅ 创建软链接: efetch"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: efetch"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: efetch"
fi

# 恢复软链接: diamond
if [ ! -e "$HOME/.local/bin/diamond" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/diamond_v.2.1.13/bin/diamond" "$HOME/.local/bin/diamond" 2>/dev/null; then
        echo "✅ 创建软链接: diamond"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: diamond"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: diamond"
fi

# 恢复软链接: matlock
if [ ! -e "$HOME/.local/bin/matlock" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/juicer_v.1.6/bin/matlock" "$HOME/.local/bin/matlock" 2>/dev/null; then
        echo "✅ 创建软链接: matlock"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: matlock"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: matlock"
fi

# 恢复软链接: hic2cool
if [ ! -e "$HOME/.local/bin/hic2cool" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hic2cool" "$HOME/.local/bin/hic2cool" 2>/dev/null; then
        echo "✅ 创建软链接: hic2cool"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hic2cool"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hic2cool"
fi

# 恢复软链接: hicInterIntraTAD
if [ ! -e "$HOME/.local/bin/hicInterIntraTAD" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicInterIntraTAD" "$HOME/.local/bin/hicInterIntraTAD" 2>/dev/null; then
        echo "✅ 创建软链接: hicInterIntraTAD"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicInterIntraTAD"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicInterIntraTAD"
fi

# 恢复软链接: smudgeplot
if [ ! -e "$HOME/.local/bin/smudgeplot" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/smudgeplot" "$HOME/.local/bin/smudgeplot" 2>/dev/null; then
        echo "✅ 创建软链接: smudgeplot"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: smudgeplot"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: smudgeplot"
fi

# 恢复软链接: cd-hit-est
if [ ! -e "$HOME/.local/bin/cd-hit-est" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/cd-hit-est" "$HOME/.local/bin/cd-hit-est" 2>/dev/null; then
        echo "✅ 创建软链接: cd-hit-est"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: cd-hit-est"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: cd-hit-est"
fi

# 恢复软链接: parallel
if [ ! -e "$HOME/.local/bin/parallel" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/juicer_v.1.6/bin/parallel" "$HOME/.local/bin/parallel" 2>/dev/null; then
        echo "✅ 创建软链接: parallel"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: parallel"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: parallel"
fi

# 恢复软链接: Trinity
if [ ! -e "$HOME/.local/bin/Trinity" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/trinity_v.2.15.2/bin/Trinity" "$HOME/.local/bin/Trinity" 2>/dev/null; then
        echo "✅ 创建软链接: Trinity"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Trinity"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Trinity"
fi

# 恢复软链接: run_augustus_multi_rnaseq
if [ ! -e "$HOME/.local/bin/run_augustus_multi_rnaseq" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_augustus_multi_rnaseq" "$HOME/.local/bin/run_augustus_multi_rnaseq" 2>/dev/null; then
        echo "✅ 创建软链接: run_augustus_multi_rnaseq"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_augustus_multi_rnaseq"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_augustus_multi_rnaseq"
fi

# 恢复软链接: mutplot
if [ ! -e "$HOME/.local/bin/mutplot" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mutmap/bin/mutplot" "$HOME/.local/bin/mutplot" 2>/dev/null; then
        echo "✅ 创建软链接: mutplot"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mutplot"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mutplot"
fi

# 恢复软链接: kat
if [ ! -e "$HOME/.local/bin/kat" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kat_v.2.4.2/bin/kat" "$HOME/.local/bin/kat" 2>/dev/null; then
        echo "✅ 创建软链接: kat"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: kat"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: kat"
fi

# 恢复软链接: tblastn
if [ ! -e "$HOME/.local/bin/tblastn" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Blast_v.2.16.0/bin/tblastn" "$HOME/.local/bin/tblastn" 2>/dev/null; then
        echo "✅ 创建软链接: tblastn"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: tblastn"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: tblastn"
fi

# 恢复软链接: TIR-Learner
if [ ! -e "$HOME/.local/bin/TIR-Learner" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/TIR-Learner" "$HOME/.local/bin/TIR-Learner" 2>/dev/null; then
        echo "✅ 创建软链接: TIR-Learner"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: TIR-Learner"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: TIR-Learner"
fi

# 恢复软链接: ragtag_scaffold.py
if [ ! -e "$HOME/.local/bin/ragtag_scaffold.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_scaffold.py" "$HOME/.local/bin/ragtag_scaffold.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_scaffold.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_scaffold.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_scaffold.py"
fi

# 恢复软链接: parallel-fastq-dump
if [ ! -e "$HOME/.local/bin/parallel-fastq-dump" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/sratoolkit_v.2.5.7/bin/parallel-fastq-dump" "$HOME/.local/bin/parallel-fastq-dump" 2>/dev/null; then
        echo "✅ 创建软链接: parallel-fastq-dump"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: parallel-fastq-dump"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: parallel-fastq-dump"
fi

# 恢复软链接: run_kaks_calculator
if [ ! -e "$HOME/.local/bin/run_kaks_calculator" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_kaks_calculator" "$HOME/.local/bin/run_kaks_calculator" 2>/dev/null; then
        echo "✅ 创建软链接: run_kaks_calculator"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_kaks_calculator"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_kaks_calculator"
fi

# 恢复软链接: run_bismark
if [ ! -e "$HOME/.local/bin/run_bismark" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_bismark" "$HOME/.local/bin/run_bismark" 2>/dev/null; then
        echo "✅ 创建软链接: run_bismark"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_bismark"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_bismark"
fi

# 恢复软链接: metawrap
if [ ! -e "$HOME/.local/bin/metawrap" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/metaWRAP/bin/metawrap" "$HOME/.local/bin/metawrap" 2>/dev/null; then
        echo "✅ 创建软链接: metawrap"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: metawrap"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: metawrap"
fi

# 恢复软链接: Fastrm
if [ ! -e "$HOME/.local/bin/Fastrm" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/smudgeplot/bin/Fastrm" "$HOME/.local/bin/Fastrm" 2>/dev/null; then
        echo "✅ 创建软链接: Fastrm"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Fastrm"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Fastrm"
fi

# 恢复软链接: seqkit
if [ ! -e "$HOME/.local/bin/seqkit" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/BioinfTools/bin/seqkit" "$HOME/.local/bin/seqkit" 2>/dev/null; then
        echo "✅ 创建软链接: seqkit"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: seqkit"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: seqkit"
fi

# 恢复软链接: retrieve_seq_from_fasta.pl
if [ ! -e "$HOME/.local/bin/retrieve_seq_from_fasta.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/annovar/annovar/retrieve_seq_from_fasta.pl" "$HOME/.local/bin/retrieve_seq_from_fasta.pl" 2>/dev/null; then
        echo "✅ 创建软链接: retrieve_seq_from_fasta.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: retrieve_seq_from_fasta.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: retrieve_seq_from_fasta.pl"
fi

# 恢复软链接: ragtag_create_links.py
if [ ! -e "$HOME/.local/bin/ragtag_create_links.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_create_links.py" "$HOME/.local/bin/ragtag_create_links.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_create_links.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_create_links.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_create_links.py"
fi

# 恢复软链接: canu
if [ ! -e "$HOME/.local/bin/canu" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/canu_v.2.3/bin/canu" "$HOME/.local/bin/canu" 2>/dev/null; then
        echo "✅ 创建软链接: canu"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: canu"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: canu"
fi

# 恢复软链接: makeblastdb
if [ ! -e "$HOME/.local/bin/makeblastdb" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Blast_v.2.16.0/bin/makeblastdb" "$HOME/.local/bin/makeblastdb" 2>/dev/null; then
        echo "✅ 创建软链接: makeblastdb"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: makeblastdb"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: makeblastdb"
fi

# 恢复软链接: run_minimap2
if [ ! -e "$HOME/.local/bin/run_minimap2" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_minimap2" "$HOME/.local/bin/run_minimap2" 2>/dev/null; then
        echo "✅ 创建软链接: run_minimap2"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_minimap2"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_minimap2"
fi

# 恢复软链接: hicTADClassifier
if [ ! -e "$HOME/.local/bin/hicTADClassifier" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicTADClassifier" "$HOME/.local/bin/hicTADClassifier" 2>/dev/null; then
        echo "✅ 创建软链接: hicTADClassifier"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicTADClassifier"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicTADClassifier"
fi

# 恢复软链接: EDTA.pl
if [ ! -e "$HOME/.local/bin/EDTA.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/EDTA.pl" "$HOME/.local/bin/EDTA.pl" 2>/dev/null; then
        echo "✅ 创建软链接: EDTA.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: EDTA.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: EDTA.pl"
fi

# 恢复软链接: Launch_PASA_pipeline.pl
if [ ! -e "$HOME/.local/bin/Launch_PASA_pipeline.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/pasa_v.2.5.3/bin/Launch_PASA_pipeline.pl" "$HOME/.local/bin/Launch_PASA_pipeline.pl" 2>/dev/null; then
        echo "✅ 创建软链接: Launch_PASA_pipeline.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Launch_PASA_pipeline.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Launch_PASA_pipeline.pl"
fi

# 恢复软链接: cafe5
if [ ! -e "$HOME/.local/bin/cafe5" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/cafe_v.5.1.0/bin/cafe5" "$HOME/.local/bin/cafe5" 2>/dev/null; then
        echo "✅ 创建软链接: cafe5"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: cafe5"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: cafe5"
fi

# 恢复软链接: java
if [ ! -e "$HOME/.local/bin/java" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/GATK_v.4.6.2.0/bin/java" "$HOME/.local/bin/java" 2>/dev/null; then
        echo "✅ 创建软链接: java"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: java"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: java"
fi

# 恢复软链接: parse_sequence_vcf
if [ ! -e "$HOME/.local/bin/parse_sequence_vcf" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/parse_sequence_vcf" "$HOME/.local/bin/parse_sequence_vcf" 2>/dev/null; then
        echo "✅ 创建软链接: parse_sequence_vcf"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: parse_sequence_vcf"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: parse_sequence_vcf"
fi

# 恢复软链接: eza
if [ ! -e "$HOME/.local/bin/eza" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/eza/eza" "$HOME/.local/bin/eza" 2>/dev/null; then
        echo "✅ 创建软链接: eza"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: eza"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: eza"
fi

# 恢复软链接: gff3ToGenePred
if [ ! -e "$HOME/.local/bin/gff3ToGenePred" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/BioinfTools/bin/gff3ToGenePred" "$HOME/.local/bin/gff3ToGenePred" 2>/dev/null; then
        echo "✅ 创建软链接: gff3ToGenePred"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: gff3ToGenePred"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: gff3ToGenePred"
fi

# 恢复软链接: ragtag_break_query.py
if [ ! -e "$HOME/.local/bin/ragtag_break_query.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_break_query.py" "$HOME/.local/bin/ragtag_break_query.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_break_query.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_break_query.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_break_query.py"
fi

# 恢复软链接: ragtag_merge.py
if [ ! -e "$HOME/.local/bin/ragtag_merge.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_merge.py" "$HOME/.local/bin/ragtag_merge.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_merge.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_merge.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_merge.py"
fi

# 恢复软链接: hicQC
if [ ! -e "$HOME/.local/bin/hicQC" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicQC" "$HOME/.local/bin/hicQC" 2>/dev/null; then
        echo "✅ 创建软链接: hicQC"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicQC"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicQC"
fi

# 恢复软链接: run_split_fasta_id
if [ ! -e "$HOME/.local/bin/run_split_fasta_id" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_split_fasta_id" "$HOME/.local/bin/run_split_fasta_id" 2>/dev/null; then
        echo "✅ 创建软链接: run_split_fasta_id"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_split_fasta_id"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_split_fasta_id"
fi

# 恢复软链接: bismark2report
if [ ! -e "$HOME/.local/bin/bismark2report" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/bismark2report" "$HOME/.local/bin/bismark2report" 2>/dev/null; then
        echo "✅ 创建软链接: bismark2report"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bismark2report"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bismark2report"
fi

# 恢复软链接: Rscript
if [ ! -e "$HOME/.local/bin/Rscript" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/R_v.4.5.1/bin/Rscript" "$HOME/.local/bin/Rscript" 2>/dev/null; then
        echo "✅ 创建软链接: Rscript"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: Rscript"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: Rscript"
fi

# 恢复软链接: bismark2summary
if [ ! -e "$HOME/.local/bin/bismark2summary" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/bismark2summary" "$HOME/.local/bin/bismark2summary" 2>/dev/null; then
        echo "✅ 创建软链接: bismark2summary"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bismark2summary"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bismark2summary"
fi

# 恢复软链接: exa
if [ ! -e "$HOME/.local/bin/exa" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/exa/bin/exa" "$HOME/.local/bin/exa" 2>/dev/null; then
        echo "✅ 创建软链接: exa"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: exa"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: exa"
fi

# 恢复软链接: hicPlotViewpoint
if [ ! -e "$HOME/.local/bin/hicPlotViewpoint" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicPlotViewpoint" "$HOME/.local/bin/hicPlotViewpoint" 2>/dev/null; then
        echo "✅ 创建软链接: hicPlotViewpoint"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicPlotViewpoint"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicPlotViewpoint"
fi

# 恢复软链接: fastme
if [ ! -e "$HOME/.local/bin/fastme" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/fastme_v.2.1.6.3/bin/fastme" "$HOME/.local/bin/fastme" 2>/dev/null; then
        echo "✅ 创建软链接: fastme"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fastme"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fastme"
fi

# 恢复软链接: hicInfo
if [ ! -e "$HOME/.local/bin/hicInfo" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicInfo" "$HOME/.local/bin/hicInfo" 2>/dev/null; then
        echo "✅ 创建软链接: hicInfo"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicInfo"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicInfo"
fi

# 恢复软链接: hicPlotSVL
if [ ! -e "$HOME/.local/bin/hicPlotSVL" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicPlotSVL" "$HOME/.local/bin/hicPlotSVL" 2>/dev/null; then
        echo "✅ 创建软链接: hicPlotSVL"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicPlotSVL"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicPlotSVL"
fi

# 恢复软链接: EDTA_processK.pl
if [ ! -e "$HOME/.local/bin/EDTA_processK.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/EDTA/bin/EDTA_processK.pl" "$HOME/.local/bin/EDTA_processK.pl" 2>/dev/null; then
        echo "✅ 创建软链接: EDTA_processK.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: EDTA_processK.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: EDTA_processK.pl"
fi

# 恢复软链接: hicFindRestSite
if [ ! -e "$HOME/.local/bin/hicFindRestSite" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicFindRestSite" "$HOME/.local/bin/hicFindRestSite" 2>/dev/null; then
        echo "✅ 创建软链接: hicFindRestSite"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicFindRestSite"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicFindRestSite"
fi

# 恢复软链接: fasterq-dump
if [ ! -e "$HOME/.local/bin/fasterq-dump" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/sratoolkit_v.2.5.7/bin/fasterq-dump" "$HOME/.local/bin/fasterq-dump" 2>/dev/null; then
        echo "✅ 创建软链接: fasterq-dump"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fasterq-dump"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fasterq-dump"
fi

# 恢复软链接: mcl
if [ ! -e "$HOME/.local/bin/mcl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Orthofinder_v.3.0.1b1/bin/mcl" "$HOME/.local/bin/mcl" 2>/dev/null; then
        echo "✅ 创建软链接: mcl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mcl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mcl"
fi

# 恢复软链接: delta-filter
if [ ! -e "$HOME/.local/bin/delta-filter" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mummer_v.3.23/bin/delta-filter" "$HOME/.local/bin/delta-filter" 2>/dev/null; then
        echo "✅ 创建软链接: delta-filter"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: delta-filter"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: delta-filter"
fi

# 恢复软链接: show-snps
if [ ! -e "$HOME/.local/bin/show-snps" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mummer_v.4.0.1/bin/show-snps" "$HOME/.local/bin/show-snps" 2>/dev/null; then
        echo "✅ 创建软链接: show-snps"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: show-snps"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: show-snps"
fi

# 恢复软链接: admixture
if [ ! -e "$HOME/.local/bin/admixture" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Population_genetics/bin/admixture" "$HOME/.local/bin/admixture" 2>/dev/null; then
        echo "✅ 创建软链接: admixture"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: admixture"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: admixture"
fi

# 恢复软链接: prodigal
if [ ! -e "$HOME/.local/bin/prodigal" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/checkm_v.1.1.0/bin/prodigal" "$HOME/.local/bin/prodigal" 2>/dev/null; then
        echo "✅ 创建软链接: prodigal"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: prodigal"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: prodigal"
fi

# 恢复软链接: VCF2Dis
if [ ! -e "$HOME/.local/bin/VCF2Dis" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/vcf2dis/VCF2Dis-1.54/bin/VCF2Dis" "$HOME/.local/bin/VCF2Dis" 2>/dev/null; then
        echo "✅ 创建软链接: VCF2Dis"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: VCF2Dis"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: VCF2Dis"
fi

# 恢复软链接: hicPlotTADs
if [ ! -e "$HOME/.local/bin/hicPlotTADs" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicPlotTADs" "$HOME/.local/bin/hicPlotTADs" 2>/dev/null; then
        echo "✅ 创建软链接: hicPlotTADs"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicPlotTADs"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicPlotTADs"
fi

# 恢复软链接: biopytools
if [ ! -e "$HOME/.local/bin/biopytools" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/biopytools" "$HOME/.local/bin/biopytools" 2>/dev/null; then
        echo "✅ 创建软链接: biopytools"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: biopytools"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: biopytools"
fi

# 恢复软链接: baseml
if [ ! -e "$HOME/.local/bin/baseml" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/paml_v.4.10.9/bin/baseml" "$HOME/.local/bin/baseml" 2>/dev/null; then
        echo "✅ 创建软链接: baseml"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: baseml"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: baseml"
fi

# 恢复软链接: genomescope2
if [ ! -e "$HOME/.local/bin/genomescope2" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/genomescope_v.2.0.1/bin/genomescope2" "$HOME/.local/bin/genomescope2" 2>/dev/null; then
        echo "✅ 创建软链接: genomescope2"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: genomescope2"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: genomescope2"
fi

# 恢复软链接: quast
if [ ! -e "$HOME/.local/bin/quast" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/quast_v.5.3.0/bin/quast" "$HOME/.local/bin/quast" 2>/dev/null; then
        echo "✅ 创建软链接: quast"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: quast"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: quast"
fi

# 恢复软链接: nw_reroot
if [ ! -e "$HOME/.local/bin/nw_reroot" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/newick_utils_v.1.6/bin/nw_reroot" "$HOME/.local/bin/nw_reroot" 2>/dev/null; then
        echo "✅ 创建软链接: nw_reroot"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: nw_reroot"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: nw_reroot"
fi

# 恢复软链接: fastp
if [ ! -e "$HOME/.local/bin/fastp" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/fastp" "$HOME/.local/bin/fastp" 2>/dev/null; then
        echo "✅ 创建软链接: fastp"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fastp"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fastp"
fi

# 恢复软链接: hmmsearch
if [ ! -e "$HOME/.local/bin/hmmsearch" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch" "$HOME/.local/bin/hmmsearch" 2>/dev/null; then
        echo "✅ 创建软链接: hmmsearch"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hmmsearch"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hmmsearch"
fi

# 恢复软链接: hicBuildMatrixMicroC
if [ ! -e "$HOME/.local/bin/hicBuildMatrixMicroC" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicBuildMatrixMicroC" "$HOME/.local/bin/hicBuildMatrixMicroC" 2>/dev/null; then
        echo "✅ 创建软链接: hicBuildMatrixMicroC"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicBuildMatrixMicroC"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicBuildMatrixMicroC"
fi

# 恢复软链接: HiC-Pro
if [ ! -e "$HOME/.local/bin/HiC-Pro" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/HiC-Pro_3.1.0/bin/HiC-Pro" "$HOME/.local/bin/HiC-Pro" 2>/dev/null; then
        echo "✅ 创建软链接: HiC-Pro"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: HiC-Pro"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: HiC-Pro"
fi

# 恢复软链接: samblaster
if [ ! -e "$HOME/.local/bin/samblaster" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/samblaster" "$HOME/.local/bin/samblaster" 2>/dev/null; then
        echo "✅ 创建软链接: samblaster"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: samblaster"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: samblaster"
fi

# 恢复软链接: hicCreateThresholdFile
if [ ! -e "$HOME/.local/bin/hicCreateThresholdFile" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicCreateThresholdFile" "$HOME/.local/bin/hicCreateThresholdFile" 2>/dev/null; then
        echo "✅ 创建软链接: hicCreateThresholdFile"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicCreateThresholdFile"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicCreateThresholdFile"
fi

# 恢复软链接: hicexplorer
if [ ! -e "$HOME/.local/bin/hicexplorer" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicexplorer" "$HOME/.local/bin/hicexplorer" 2>/dev/null; then
        echo "✅ 创建软链接: hicexplorer"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicexplorer"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicexplorer"
fi

# 恢复软链接: bismark_methylation_extractor
if [ ! -e "$HOME/.local/bin/bismark_methylation_extractor" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bismark_v.0.24.2/bin/bismark_methylation_extractor" "$HOME/.local/bin/bismark_methylation_extractor" 2>/dev/null; then
        echo "✅ 创建软链接: bismark_methylation_extractor"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bismark_methylation_extractor"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bismark_methylation_extractor"
fi

# 恢复软链接: parse_longest_mrna
if [ ! -e "$HOME/.local/bin/parse_longest_mrna" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/parse_longest_mrna" "$HOME/.local/bin/parse_longest_mrna" 2>/dev/null; then
        echo "✅ 创建软链接: parse_longest_mrna"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: parse_longest_mrna"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: parse_longest_mrna"
fi

# 恢复软链接: generate_font_width_hash.pl
if [ ! -e "$HOME/.local/bin/generate_font_width_hash.pl" ]; then
    if ln -s "generate_font_width_hash.pl" "$HOME/.local/bin/generate_font_width_hash.pl" 2>/dev/null; then
        echo "✅ 创建软链接: generate_font_width_hash.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: generate_font_width_hash.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: generate_font_width_hash.pl"
fi

# 恢复软链接: run_vcf_sample_hete
if [ ! -e "$HOME/.local/bin/run_vcf_sample_hete" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_vcf_sample_hete" "$HOME/.local/bin/run_vcf_sample_hete" 2>/dev/null; then
        echo "✅ 创建软链接: run_vcf_sample_hete"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_vcf_sample_hete"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_vcf_sample_hete"
fi

# 恢复软链接: get_novogene
if [ ! -e "$HOME/.local/bin/get_novogene" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/novogene/linuxnd/get_novogene" "$HOME/.local/bin/get_novogene" 2>/dev/null; then
        echo "✅ 创建软链接: get_novogene"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: get_novogene"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: get_novogene"
fi

# 恢复软链接: rpsblast
if [ ! -e "$HOME/.local/bin/rpsblast" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Blast_v.2.16.0/bin/rpsblast" "$HOME/.local/bin/rpsblast" 2>/dev/null; then
        echo "✅ 创建软链接: rpsblast"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: rpsblast"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: rpsblast"
fi

# 恢复软链接: bbmap.sh
if [ ! -e "$HOME/.local/bin/bbmap.sh" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bbtools_v.37.62/bin/bbmap.sh" "$HOME/.local/bin/bbmap.sh" 2>/dev/null; then
        echo "✅ 创建软链接: bbmap.sh"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bbmap.sh"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bbmap.sh"
fi

# 恢复软链接: interproscan.sh
if [ ! -e "$HOME/.local/bin/interproscan.sh" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh" "$HOME/.local/bin/interproscan.sh" 2>/dev/null; then
        echo "✅ 创建软链接: interproscan.sh"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: interproscan.sh"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: interproscan.sh"
fi

# 恢复软链接: resistify
if [ ! -e "$HOME/.local/bin/resistify" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/resistify_v.1.3.0/bin/resistify" "$HOME/.local/bin/resistify" 2>/dev/null; then
        echo "✅ 创建软链接: resistify"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: resistify"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: resistify"
fi

# 恢复软链接: bam2hints
if [ ! -e "$HOME/.local/bin/bam2hints" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/bamtools/bin/bam2hints" "$HOME/.local/bin/bam2hints" 2>/dev/null; then
        echo "✅ 创建软链接: bam2hints"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bam2hints"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bam2hints"
fi

# 恢复软链接: haphic
if [ ! -e "$HOME/.local/bin/haphic" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/haphic" "$HOME/.local/bin/haphic" 2>/dev/null; then
        echo "✅ 创建软链接: haphic"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: haphic"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: haphic"
fi

# 恢复软链接: in2csv
if [ ! -e "$HOME/.local/bin/in2csv" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/csvkit/bin/in2csv" "$HOME/.local/bin/in2csv" 2>/dev/null; then
        echo "✅ 创建软链接: in2csv"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: in2csv"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: in2csv"
fi

# 恢复软链接: fastqc
if [ ! -e "$HOME/.local/bin/fastqc" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/fastqc_v.0.12.1/bin/fastqc" "$HOME/.local/bin/fastqc" 2>/dev/null; then
        echo "✅ 创建软链接: fastqc"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fastqc"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fastqc"
fi

# 恢复软链接: GetTwoGenomeSyn.pl
if [ ! -e "$HOME/.local/bin/GetTwoGenomeSyn.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/NGenomeSyn/NGenomeSyn-1.43/bin/GetTwoGenomeSyn.pl" "$HOME/.local/bin/GetTwoGenomeSyn.pl" 2>/dev/null; then
        echo "✅ 创建软链接: GetTwoGenomeSyn.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: GetTwoGenomeSyn.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: GetTwoGenomeSyn.pl"
fi

# 恢复软链接: axel
if [ ! -e "$HOME/.local/bin/axel" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/axel_v.2.17.13/bin/axel" "$HOME/.local/bin/axel" 2>/dev/null; then
        echo "✅ 创建软链接: axel"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: axel"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: axel"
fi

# 恢复软链接: run_genome_collinearity
if [ ! -e "$HOME/.local/bin/run_genome_collinearity" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_genome_collinearity" "$HOME/.local/bin/run_genome_collinearity" 2>/dev/null; then
        echo "✅ 创建软链接: run_genome_collinearity"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_genome_collinearity"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_genome_collinearity"
fi

# 恢复软链接: run_vcf_genotype
if [ ! -e "$HOME/.local/bin/run_vcf_genotype" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_vcf_genotype" "$HOME/.local/bin/run_vcf_genotype" 2>/dev/null; then
        echo "✅ 创建软链接: run_vcf_genotype"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_vcf_genotype"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_vcf_genotype"
fi

# 恢复软链接: PopLDdecay
if [ ! -e "$HOME/.local/bin/PopLDdecay" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/poplddecay_v.3.43/bin/PopLDdecay" "$HOME/.local/bin/PopLDdecay" 2>/dev/null; then
        echo "✅ 创建软链接: PopLDdecay"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: PopLDdecay"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: PopLDdecay"
fi

# 恢复软链接: parse_sample_hete
if [ ! -e "$HOME/.local/bin/parse_sample_hete" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/parse_sample_hete" "$HOME/.local/bin/parse_sample_hete" 2>/dev/null; then
        echo "✅ 创建软链接: parse_sample_hete"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: parse_sample_hete"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: parse_sample_hete"
fi

# 恢复软链接: metagraph
if [ ! -e "$HOME/.local/bin/metagraph" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/metagraph/bin/metagraph" "$HOME/.local/bin/metagraph" 2>/dev/null; then
        echo "✅ 创建软链接: metagraph"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: metagraph"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: metagraph"
fi

# 恢复软链接: pbsv
if [ ! -e "$HOME/.local/bin/pbsv" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/pbsv_v.2.11.0/bin/pbsv" "$HOME/.local/bin/pbsv" 2>/dev/null; then
        echo "✅ 创建软链接: pbsv"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: pbsv"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: pbsv"
fi

# 恢复软链接: ragtag_delta2paf.py
if [ ! -e "$HOME/.local/bin/ragtag_delta2paf.py" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RagTag_v2.10./bin/ragtag_delta2paf.py" "$HOME/.local/bin/ragtag_delta2paf.py" 2>/dev/null; then
        echo "✅ 创建软链接: ragtag_delta2paf.py"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: ragtag_delta2paf.py"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: ragtag_delta2paf.py"
fi

# 恢复软链接: bracken
if [ ! -e "$HOME/.local/bin/bracken" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/kraken_v.2.17/bin/bracken" "$HOME/.local/bin/bracken" 2>/dev/null; then
        echo "✅ 创建软链接: bracken"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bracken"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bracken"
fi

# 恢复软链接: tiberius
if [ ! -e "$HOME/.local/bin/tiberius" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Tiberius_v.1.1.1/bin/tiberius" "$HOME/.local/bin/tiberius" 2>/dev/null; then
        echo "✅ 创建软链接: tiberius"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: tiberius"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: tiberius"
fi

# 恢复软链接: run_ena_downloader
if [ ! -e "$HOME/.local/bin/run_ena_downloader" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_ena_downloader" "$HOME/.local/bin/run_ena_downloader" 2>/dev/null; then
        echo "✅ 创建软链接: run_ena_downloader"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_ena_downloader"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_ena_downloader"
fi

# 恢复软链接: julia
if [ ! -e "$HOME/.local/bin/julia" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/julia_v.1.12.2/bin/julia" "$HOME/.local/bin/julia" 2>/dev/null; then
        echo "✅ 创建软链接: julia"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: julia"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: julia"
fi

# 恢复软链接: EVidenceModeler
if [ ! -e "$HOME/.local/bin/EVidenceModeler" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/evidencemodeler/bin/EVidenceModeler" "$HOME/.local/bin/EVidenceModeler" 2>/dev/null; then
        echo "✅ 创建软链接: EVidenceModeler"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: EVidenceModeler"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: EVidenceModeler"
fi

# 恢复软链接: hicAverageRegions
if [ ! -e "$HOME/.local/bin/hicAverageRegions" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicAverageRegions" "$HOME/.local/bin/hicAverageRegions" 2>/dev/null; then
        echo "✅ 创建软链接: hicAverageRegions"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicAverageRegions"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicAverageRegions"
fi

# 恢复软链接: hicValidateLocations
if [ ! -e "$HOME/.local/bin/hicValidateLocations" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicValidateLocations" "$HOME/.local/bin/hicValidateLocations" 2>/dev/null; then
        echo "✅ 创建软链接: hicValidateLocations"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicValidateLocations"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicValidateLocations"
fi

# 恢复软链接: bgzip
if [ ! -e "$HOME/.local/bin/bgzip" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/RNA_Seq/bin/bgzip" "$HOME/.local/bin/bgzip" 2>/dev/null; then
        echo "✅ 创建软链接: bgzip"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: bgzip"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: bgzip"
fi

# 恢复软链接: run_fastp
if [ ! -e "$HOME/.local/bin/run_fastp" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_fastp" "$HOME/.local/bin/run_fastp" 2>/dev/null; then
        echo "✅ 创建软链接: run_fastp"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_fastp"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_fastp"
fi

# 恢复软链接: hicMergeLoops
if [ ! -e "$HOME/.local/bin/hicMergeLoops" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicMergeLoops" "$HOME/.local/bin/hicMergeLoops" 2>/dev/null; then
        echo "✅ 创建软链接: hicMergeLoops"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicMergeLoops"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicMergeLoops"
fi

# 恢复软链接: table_annovar.pl
if [ ! -e "$HOME/.local/bin/table_annovar.pl" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/annovar/annovar/table_annovar.pl" "$HOME/.local/bin/table_annovar.pl" 2>/dev/null; then
        echo "✅ 创建软链接: table_annovar.pl"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: table_annovar.pl"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: table_annovar.pl"
fi

# 恢复软链接: hicMergeMatrixBins
if [ ! -e "$HOME/.local/bin/hicMergeMatrixBins" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/hicexplorer/bin/hicMergeMatrixBins" "$HOME/.local/bin/hicMergeMatrixBins" 2>/dev/null; then
        echo "✅ 创建软链接: hicMergeMatrixBins"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: hicMergeMatrixBins"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: hicMergeMatrixBins"
fi

# 恢复软链接: fastq-dump
if [ ! -e "$HOME/.local/bin/fastq-dump" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/sratoolkit_v.2.5.7/bin/fastq-dump" "$HOME/.local/bin/fastq-dump" 2>/dev/null; then
        echo "✅ 创建软链接: fastq-dump"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: fastq-dump"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: fastq-dump"
fi

# 恢复软链接: run_plink_gwas
if [ ! -e "$HOME/.local/bin/run_plink_gwas" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_plink_gwas" "$HOME/.local/bin/run_plink_gwas" 2>/dev/null; then
        echo "✅ 创建软链接: run_plink_gwas"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_plink_gwas"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_plink_gwas"
fi

# 恢复软链接: mummer
if [ ! -e "$HOME/.local/bin/mummer" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/mummer_v.3.23/bin/mummer" "$HOME/.local/bin/mummer" 2>/dev/null; then
        echo "✅ 创建软链接: mummer"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: mummer"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: mummer"
fi

# 恢复软链接: pandepth
if [ ! -e "$HOME/.local/bin/pandepth" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/software/PanDepth-2.26-Linux-x86_64/pandepth" "$HOME/.local/bin/pandepth" 2>/dev/null; then
        echo "✅ 创建软链接: pandepth"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: pandepth"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: pandepth"
fi

# 恢复软链接: megahit
if [ ! -e "$HOME/.local/bin/megahit" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/metaWRAP_v.1.2/bin/megahit" "$HOME/.local/bin/megahit" 2>/dev/null; then
        echo "✅ 创建软链接: megahit"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: megahit"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: megahit"
fi

# 恢复软链接: maker
if [ ! -e "$HOME/.local/bin/maker" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/NLR_Annotation_Pipeline/bin/maker" "$HOME/.local/bin/maker" 2>/dev/null; then
        echo "✅ 创建软链接: maker"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: maker"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: maker"
fi

# 恢复软链接: run_variant_filter
if [ ! -e "$HOME/.local/bin/run_variant_filter" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/biopytools/bin/run_variant_filter" "$HOME/.local/bin/run_variant_filter" 2>/dev/null; then
        echo "✅ 创建软链接: run_variant_filter"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: run_variant_filter"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: run_variant_filter"
fi

# 恢复软链接: gffcompare
if [ ! -e "$HOME/.local/bin/gffcompare" ]; then
    if ln -s "/share/org/YZWL/yzwl_lixg/miniforge3/envs/BioinfTools/bin/gffcompare" "$HOME/.local/bin/gffcompare" 2>/dev/null; then
        echo "✅ 创建软链接: gffcompare"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: gffcompare"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: gffcompare"
fi

# 恢复普通文件: starship  
if [ ! -f "$HOME/.local/bin/starship" ]; then
    if cp "$(dirname "$0")/starship" "$HOME/.local/bin/starship" 2>/dev/null; then
        chmod +x "$HOME/.local/bin/starship" 2>/dev/null || true
        echo "✅ 复制文件: starship"
        ((success_count++))
    else
        echo "❌ 复制文件失败: starship"
        ((error_count++))
    fi
else
    echo "⚠️  文件已存在: starship"
fi

# 恢复普通文件: export_all_conda_envs.py  
if [ ! -f "$HOME/.local/bin/export_all_conda_envs.py" ]; then
    if cp "$(dirname "$0")/export_all_conda_envs.py" "$HOME/.local/bin/export_all_conda_envs.py" 2>/dev/null; then
        chmod +x "$HOME/.local/bin/export_all_conda_envs.py" 2>/dev/null || true
        echo "✅ 复制文件: export_all_conda_envs.py"
        ((success_count++))
    else
        echo "❌ 复制文件失败: export_all_conda_envs.py"
        ((error_count++))
    fi
else
    echo "⚠️  文件已存在: export_all_conda_envs.py"
fi

# 恢复普通文件: interproscan.sh.backup_20250920_101059  
if [ ! -f "$HOME/.local/bin/interproscan.sh.backup_20250920_101059" ]; then
    if cp "$(dirname "$0")/interproscan.sh.backup_20250920_101059" "$HOME/.local/bin/interproscan.sh.backup_20250920_101059" 2>/dev/null; then
        chmod +x "$HOME/.local/bin/interproscan.sh.backup_20250920_101059" 2>/dev/null || true
        echo "✅ 复制文件: interproscan.sh.backup_20250920_101059"
        ((success_count++))
    else
        echo "❌ 复制文件失败: interproscan.sh.backup_20250920_101059"
        ((error_count++))
    fi
else
    echo "⚠️  文件已存在: interproscan.sh.backup_20250920_101059"
fi

echo ""
echo "📊 恢复完成统计:"
echo "  成功: $success_count"
echo "  失败: $error_count"
echo "✅ 恢复脚本执行完毕"

# #!/usr/bin/env python3
# """
# PLINK GWAS分析运行脚本 | PLINK GWAS Analysis Runner Script
# 这是一个简化的入口脚本，用于运行PLINK GWAS分析 | Simple entry script for running PLINK GWAS analysis

# 用法 | Usage:
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative -o results
    
# 示例 | Examples:
#     # 质量性状分析 (0/1 -> 1/2转换) | Qualitative trait analysis (0/1 -> 1/2 conversion)
#     python run_plink.py -v my_data.vcf.gz -p my_pheno.txt -t qualitative -o plink_results
    
#     # 数量性状分析 (保持原值) | Quantitative trait analysis (keep original values)
#     python run_plink.py -v my_data.vcf.gz -p my_pheno.txt -t quantitative -o plink_results
    
#     # 使用所有校正方法 (默认) | Use all correction methods (default)
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method all
    
#     # 仅使用Bonferroni校正 | Only use Bonferroni correction
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method bonferroni
    
#     # 仅使用FDR校正 | Only use FDR correction
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method fdr --fdr-alpha 0.1
    
#     # 自定义校正参数 | Custom correction parameters
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method all \\
#         --bonferroni-alpha 0.01 --suggestive-threshold 5e-6 --fdr-alpha 0.05
    
#     # 自定义质控参数 | Custom QC parameters
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative -o results \\
#         --mind 0.1 --geno 0.1 --maf 0.05 --hwe 1e-5
    
#     # 调整主成分分析参数 | Adjust PCA parameters
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t quantitative -o results \\
#         --pca-components 15 --pca-use 8
    
#     # 使用多线程 | Use multiple threads
#     python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative -o results --threads 8
# """

# from biopytools.plink_gwas.main import main

# if __name__ == "__main__":
#     main()

# 20250726添加显性模型和隐性模型的配置选项 | Added options for dominant and recessive models
#!/usr/bin/env python3
"""
PLINK GWAS分析运行脚本 | PLINK GWAS Analysis Runner Script
这是一个增强版的入口脚本，支持多种遗传模型的GWAS分析 | Enhanced entry script supporting multi-model GWAS analysis

新增功能 | New Features:
    - 支持加性、显性、隐性和全模型分析 | Support additive, dominant, recessive and all-model analysis
    - 智能表型转换，根据遗传模型调整编码策略 | Intelligent phenotype conversion with model-specific coding
    - 多模型结果比较和可视化 | Multi-model result comparison and visualization
    - 增强的统计报告和建议 | Enhanced statistical reports and recommendations

用法 | Usage:
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m additive -o results
    
示例 | Examples:
    # 基本加性模型分析（默认）| Basic additive model analysis (default)
    python run_plink_gwas.py -v my_data.vcf.gz -p my_pheno.txt -t qualitative -o plink_results
    
    # 显性模型分析 | Dominant model analysis
    python run_plink_gwas.py -v my_data.vcf.gz -p my_pheno.txt -t qualitative -m dominant -o plink_results
    
    # 隐性模型分析 | Recessive model analysis
    python run_plink_gwas.py -v my_data.vcf.gz -p my_pheno.txt -t qualitative -m recessive -o plink_results
    
    # 测试所有遗传模型 | Test all genetic models
    python run_plink_gwas.py -v my_data.vcf.gz -p my_pheno.txt -t qualitative -m all -o plink_results
    
    # 数量性状分析 | Quantitative trait analysis
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t quantitative -m additive -o results
    
    # 使用所有校正方法 (默认) | Use all correction methods (default)
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m additive --correction-method all
    
    # 仅使用Bonferroni校正 | Only use Bonferroni correction
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m dominant --correction-method bonferroni
    
    # 仅使用FDR校正 | Only use FDR correction
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m all --correction-method fdr --fdr-alpha 0.1
    
    # 自定义校正参数 | Custom correction parameters
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m all --correction-method all \\
        --bonferroni-alpha 0.01 --suggestive-threshold 5e-6 --fdr-alpha 0.05
    
    # 自定义质控参数 | Custom QC parameters
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m additive -o results \\
        --mind 0.1 --geno 0.1 --maf 0.05 --hwe 1e-5
    
    # 调整主成分分析参数 | Adjust PCA parameters
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t quantitative -m additive -o results \\
        --pca-components 15 --pca-use 8
    
    # 使用多线程 | Use multiple threads
    python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m all -o results --threads 8

遗传模型说明 | Genetic Model Description:
    additive  : 加性模型，假设杂合子效应是两个纯合子的中间值 | Additive model, assumes heterozygote effect is intermediate
    dominant  : 显性模型，假设一个风险等位基因拷贝就足以产生效应 | Dominant model, one risk allele copy is sufficient
    recessive : 隐性模型，假设需要两个风险等位基因拷贝才产生效应 | Recessive model, two risk allele copies are needed
    all       : 测试所有模型并比较结果 | Test all models and compare results

表型转换说明 | Phenotype Conversion Notes:
    - 质量性状: 0(抗病) -> 1(对照), 1(感病) -> 2(病例) | Qualitative: 0(resistant) -> 1(control), 1(susceptible) -> 2(case)
    - 数量性状: 保持原始数值，缺失值转为-9 | Quantitative: keep original values, missing -> -9
    - 遗传模型差异主要体现在基因型编码而非表型编码 | Model differences mainly in genotype coding, not phenotype coding

输出文件说明 | Output Files Description:
    当选择单一模型时 | When single model selected:
    - gwas_results_[MODEL].txt: 主要关联分析结果 | Main association results
    - significant_hits_bonferroni_[MODEL].txt: Bonferroni校正显著位点 | Bonferroni significant loci
    - suggestive_hits_[MODEL].txt: 提示性关联位点 | Suggestive association loci
    - fdr_significant_hits_[MODEL].txt: FDR校正显著位点 | FDR significant loci
    - gwas_summary_report.txt: 总结报告 | Summary report
    
    当选择所有模型时 | When all models selected:
    - gwas_results_ADD.txt: 加性模型结果 | Additive model results
    - gwas_results_DOM.txt: 显性模型结果 | Dominant model results  
    - gwas_results_GENO.txt: 基因型模型结果 | Genotypic model results
    - model_comparison_report.txt: 模型比较报告 | Model comparison report
    - manhattan_plot_[model].png: 各模型Manhattan图 | Manhattan plots for each model
    - qq_plot_[model].png: 各模型QQ图 | QQ plots for each model
    - model_comparison_pvalue_distribution.png: 模型P值分布比较图 | Model P-value distribution comparison

F2代抗感病分析建议 | Recommendations for F2 Disease Resistance Analysis:
    - 推荐首先使用加性模型 | Recommend starting with additive model
    - 如果加性模型效果不佳，尝试显性或隐性模型 | If additive model performs poorly, try dominant/recessive
    - 使用all模型选项可以同时比较所有模型 | Use 'all' option to compare all models simultaneously
    - 注意样本量：每组至少50个样本获得可靠结果 | Note sample size: at least 50 samples per group for reliable results
    - 考虑群体分层和主成分校正 | Consider population stratification and PC correction
"""

from biopytools.plink_gwas.main import main

if __name__ == "__main__":
    main()
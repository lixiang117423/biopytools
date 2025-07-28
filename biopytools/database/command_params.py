# biopytools/database/command_params.py

PLINK_PARAM_DB = {
    "plink": {
        # ... (您plink_gwas模块的所有plink参数) ...
        "--bfile": "指定二进制PLINK文件集前缀",
        "--vcf": "指定VCF格式输入文件",
        "--make-bed": "生成二进制PLINK文件格式",
        # ... etc ...
    }
}

BWA_PARAM_DB = {
    "bwa": {
        "mem": "使用BWA-MEM算法进行比对",
        "-t": "设置线程数",
        "-R": "添加Read Group头信息 (格式: '@RG\\tID:group1\\tSM:sample1')",
        "index": "为参考基因组构建索引",
    }
}

SAMTOOLS_PARAM_DB = {
    "samtools": {
        "view": "查看和转换SAM/BAM/CRAM文件",
        "sort": "对BAM文件进行排序",
        "index": "为排序后的BAM文件创建索引",
        "flagstat": "统计BAM文件的比对信息",
        "-@": "设置线程数",
        "-b": "输出BAM格式",
        "-o": "指定输出文件名",
    }
}

# 可以继续为 GATK, STAR, HISAT2 等添加
# ...

# 将所有数据库合并到一个字典中，方便CommandRunner使用
ALL_COMMANDS_DB = {
    **PLINK_PARAM_DB,
    **BWA_PARAM_DB,
    **SAMTOOLS_PARAM_DB,
}
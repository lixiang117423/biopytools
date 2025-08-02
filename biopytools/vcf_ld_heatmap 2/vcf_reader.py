class VCFReader:
    """VCF文件读取器 | VCF File Reader"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def read_vcf_data(self) -> Tuple[allel.GenotypeArray, np.ndarray, np.ndarray, List[str]]:
        """读取VCF文件数据 | Read VCF file data"""
        if self.config.verbose:
            self.logger.info(f"正在读取VCF文件 | Reading VCF file: {self.config.vcf_file}")
        
        try:
            if self.config.region:
                # 解析区域参数 | Parse region parameter
                if ':' in self.config.region and '-' in self.config.region:
                    chrom, pos_range = self.config.region.split(':')
                    start, end = map(int, pos_range.split('-'))
                    callset = allel.read_vcf(self.config.vcf_file, region=f"{chrom}:{start}-{end}")
                else:
                    callset = allel.read_vcf(self.config.vcf_file, region=self.config.region)
            else:
                callset = allel.read_vcf(self.config.vcf_file)
                
            if callset is None:
                raise ValueError("无法读取VCF文件或指定区域无数据 | Cannot read VCF file or no data in specified region")
                
        except Exception as e:
            self.logger.error(f"错误: 读取VCF文件失败 | Error: Failed to read VCF file - {e}")
            raise
        
        # 提取基因型数据 | Extract genotype data
        gt = allel.GenotypeArray(callset['calldata/GT'])
        positions = callset['variants/POS']
        chrom = callset['variants/CHROM']
        
        # 如果指定了样本，进行过滤 | Filter samples if specified
        if self.config.samples:
            sample_indices = []
            all_samples = list(callset['samples'])
            for sample in self.config.samples:
                if sample in all_samples:
                    sample_indices.append(all_samples.index(sample))
                else:
                    self.logger.warning(f"警告: 样本 {sample} 不在VCF文件中 | Warning: Sample {sample} not in VCF file")
            
            if sample_indices:
                gt = gt.take(sample_indices, axis=1)
            else:
                raise ValueError("错误: 没有找到指定的样本 | Error: No specified samples found")
        
        if self.config.verbose:
            self.logger.info(f"读取到 {gt.shape[0]} 个变异位点, {gt.shape[1]} 个样本 | "
                           f"Read {gt.shape[0]} variants, {gt.shape[1]} samples")
        
        return gt, positions, chrom, callset.get('samples', [])

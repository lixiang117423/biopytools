"""
🌾 EDTA分析统计模块 | EDTA Analysis Statistics Module
"""

import os
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any
from collections import defaultdict, Counter
import re

class EDTAAnalysisProcessor:
    """EDTA分析结果处理器 | EDTA Analysis Results Processor"""
    
    def __init__(self, config, logger, directories: Dict[str, Path]):
        self.config = config
        self.logger = logger
        self.directories = directories
    
    def process_results(self, edta_results: Dict[str, Dict]) -> Dict[str, Any]:
        """处理EDTA结果 | Process EDTA results"""
        self.logger.info("📊 处理EDTA分析结果 | Processing EDTA analysis results")
        
        processed_results = {}
        
        for analysis_name, result_info in edta_results.items():
            if result_info["status"] == "success":
                self.logger.info(f"📋 处理成功结果 | Processing successful result: {analysis_name}")
                
                # 解析注释文件
                annotation_stats = self._parse_annotation_file(result_info, analysis_name)
                
                # 解析统计摘要
                summary_stats = self._parse_summary_file(result_info, analysis_name)
                
                processed_results[analysis_name] = {
                    "original_result": result_info,
                    "annotation_stats": annotation_stats,
                    "summary_stats": summary_stats,
                    "status": "processed"
                }
                
            else:
                self.logger.warning(f"⚠️ 跳过失败结果 | Skipping failed result: {analysis_name}")
                processed_results[analysis_name] = {
                    "original_result": result_info,
                    "status": "failed"
                }
        
        self.logger.info("✅ 结果处理完成 | Results processing completed")
        return processed_results
    
    def _parse_annotation_file(self, result_info: Dict, analysis_name: str) -> Dict[str, Any]:
        """解析GFF3注释文件 | Parse GFF3 annotation file"""
        self.logger.info(f"📄 解析注释文件 | Parsing annotation file: {analysis_name}")
        
        output_files = result_info.get("output_files", {})
        gff3_file = None
        
        # 查找GFF3文件
        for filename, filepath in output_files.items():
            if filename.endswith(".TEanno.gff3"):
                gff3_file = filepath
                break
        
        if not gff3_file or not os.path.exists(gff3_file):
            self.logger.warning(f"⚠️ GFF3注释文件不存在 | GFF3 annotation file not found: {analysis_name}")
            return {"error": "GFF3 file not found"}
        
        try:
            # 解析GFF3文件
            annotations = []
            with open(gff3_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 9:
                        annotations.append({
                            'seqname': parts[0],
                            'source': parts[1],
                            'feature': parts[2],
                            'start': int(parts[3]),
                            'end': int(parts[4]),
                            'score': parts[5],
                            'strand': parts[6],
                            'frame': parts[7],
                            'attributes': parts[8],
                            'length': int(parts[4]) - int(parts[3]) + 1
                        })
            
            # 统计分析
            stats = {
                "total_elements": len(annotations),
                "total_length": sum(ann['length'] for ann in annotations),
                "chromosome_distribution": Counter(ann['seqname'] for ann in annotations),
                "feature_distribution": Counter(ann['feature'] for ann in annotations),
                "strand_distribution": Counter(ann['strand'] for ann in annotations),
                "length_statistics": self._calculate_length_stats([ann['length'] for ann in annotations])
            }
            
            # 解析Classification属性
            classifications = []
            for ann in annotations:
                attrs = ann['attributes']
                if 'Classification=' in attrs:
                    match = re.search(r'Classification=([^;]+)', attrs)
                    if match:
                        classifications.append(match.group(1))
            
            stats["classification_distribution"] = Counter(classifications)
            
            self.logger.info(f"✅ 注释文件解析完成 | Annotation file parsing completed: {len(annotations)} elements")
            return stats
            
        except Exception as e:
            self.logger.error(f"❌ 注释文件解析失败 | Annotation file parsing failed: {e}")
            return {"error": str(e)}
    
    def _parse_summary_file(self, result_info: Dict, analysis_name: str) -> Dict[str, Any]:
        """解析统计摘要文件 | Parse summary file"""
        self.logger.info(f"📄 解析摘要文件 | Parsing summary file: {analysis_name}")
        
        output_files = result_info.get("output_files", {})
        summary_file = None
        
        # 查找摘要文件
        for filename, filepath in output_files.items():
            if filename.endswith(".anno.sum") or filename.endswith(".EDTA.anno.sum"):
                summary_file = filepath
                break
        
        if not summary_file or not os.path.exists(summary_file):
            self.logger.warning(f"⚠️ 摘要文件不存在 | Summary file not found: {analysis_name}")
            return {"error": "Summary file not found"}
        
        try:
            summary_data = {}
            with open(summary_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 解析摘要内容
            lines = content.split('\n')
            for line in lines:
                if line.strip() and 'Total' in line and 'bp' in line:
                    match = re.search(r'(\d+)\s*bp.*?(\d+\.?\d*)\s*%', line)
                    if match:
                        if 'Total length' in line:
                            summary_data['total_length_bp'] = int(match.group(1))
                            summary_data['genome_percentage'] = float(match.group(2))
            
            self.logger.info(f"✅ 摘要文件解析完成 | Summary file parsing completed")
            return summary_data
            
        except Exception as e:
            self.logger.error(f"❌ 摘要文件解析失败 | Summary file parsing failed: {e}")
            return {"error": str(e)}
    
    def _calculate_length_stats(self, lengths: List[int]) -> Dict[str, float]:
        """计算长度统计 | Calculate length statistics"""
        if not lengths:
            return {}
        
        sorted_lengths = sorted(lengths)
        n = len(sorted_lengths)
        
        return {
            "count": n,
            "min": min(sorted_lengths),
            "max": max(sorted_lengths),
            "mean": sum(sorted_lengths) / n,
            "median": sorted_lengths[n//2] if n % 2 == 1 else (sorted_lengths[n//2-1] + sorted_lengths[n//2]) / 2,
            "q1": sorted_lengths[n//4],
            "q3": sorted_lengths[3*n//4]
        }

# """
# 输出写入模块 | Output Writer Module
# """

# import pandas as pd
# from datetime import datetime
# from pathlib import Path
# from typing import List, Dict

# class OutputWriter:
#     """输出文件写入器 | Output File Writer"""
    
#     def __init__(self, logger):
#         self.logger = logger
    
#     def write_txt_output(self, variants: List[Dict], output_file: Path, 
#                         region_info: Dict, variant_type: str = "变异"):
#         """写入格式化的TXT文件 | Write formatted TXT file"""
#         with open(output_file, 'w', encoding='utf-8') as f:
#             f.write("=" * 120 + "\n")
#             f.write(f"🧬 基因区域{variant_type}提取报告 | Gene Region {variant_type} Extraction Report\n")
#             f.write("=" * 120 + "\n")
#             f.write(f"⏰ 生成时间 | Generation time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
#             f.write(f"🎯 基因ID | Gene ID: {region_info['gene_id']}\n")
#             f.write(f"📍 染色体 | Chromosome: {region_info['chrom']}\n")
#             f.write(f"📏 基因位置 | Gene location: {region_info['gene_start']}-{region_info['gene_end']} (链方向|Strand: {region_info['strand']})\n")
#             f.write(f"📏 基因长度 | Gene length: {region_info['gene_end'] - region_info['gene_start'] + 1} bp\n")
#             f.write(f"📐 分析区域 | Analysis region: {region_info['region_start']}-{region_info['region_end']} (±{region_info['extend']}bp)\n")
#             f.write(f"📐 区域总长度 | Total region length: {region_info['region_end'] - region_info['region_start'] + 1} bp\n")
#             f.write(f"📊 变异总数 | Total variants: {len(variants)}\n")
#             f.write("=" * 120 + "\n\n")
            
#             f.write("📖 坐标说明 | Coordinate explanation:\n")
#             f.write("  - 基因坐标: 相对于基因起始的1-based坐标 | Gene coordinate: 1-based relative to gene start\n")
#             if region_info['strand'] == '+':
#                 f.write("  - 正链基因: 基因坐标从5'端(左侧)开始计数 | Forward strand: counted from 5' end (left)\n")
#                 f.write("    * 上游(左侧): 负数, 如-1000 = 上游1000bp | Upstream(left): negative, e.g., -1000 = 1000bp upstream\n")
#                 f.write("    * 下游(右侧): 正数, 如+500 = 下游500bp | Downstream(right): positive, e.g., +500 = 500bp downstream\n")
#             else:
#                 f.write("  - 负链基因: 基因坐标从5'端(右侧)开始计数 | Reverse strand: counted from 5' end (right)\n")
#                 f.write("    * 上游(右侧): 负数, 如-1000 = 上游1000bp | Upstream(right): negative, e.g., -1000 = 1000bp upstream\n")
#                 f.write("    * 下游(左侧): 正数, 如+500 = 下游500bp | Downstream(left): positive, e.g., +500 = 500bp downstream\n")
#             f.write("=" * 120 + "\n\n")
            
#             if len(variants) == 0:
#                 f.write("ℹ️  未找到变异 | No variants found.\n")
#                 return
            
#             # 写入变异表格 | Write variant table
#             f.write(f"{'序号':<6} {'基因组位置':<12} {'基因坐标':<12} {'位置类型':<10} {'参考':<8} {'变异':<8} {'类型/功能':<20} 基因信息\n")
#             f.write("-" * 120 + "\n")
            
#             for idx, var in enumerate(variants, 1):
#                 pos = var['position']
#                 gene_coord = var['gene_coordinate']
#                 location = var['location_type']
#                 ref = var['ref']
#                 alt = var['alt']
                
#                 # 位置类型显示 | Location type display
#                 if location == 'gene':
#                     loc_display = '基因内'
#                 elif location == 'upstream':
#                     loc_display = '上游'
#                 else:
#                     loc_display = '下游'
                
#                 if 'variant_type' in var:  # exonic
#                     type_info = var['variant_type']
#                 else:  # all variants
#                     type_info = var['function']
                
#                 gene_info = var['gene_info'][:25] + '...' if len(var['gene_info']) > 25 else var['gene_info']
                
#                 f.write(f"{idx:<6} {pos:<12} {gene_coord:<12} {loc_display:<10} {ref:<8} {alt:<8} {type_info:<20} {gene_info}\n")
            
#             f.write("\n" + "=" * 120 + "\n")
#             f.write("📝 详细信息（完整行）| Detailed information (complete lines)\n")
#             f.write("=" * 120 + "\n\n")
            
#             for idx, var in enumerate(variants, 1):
#                 f.write(f"[{idx}] 基因组:{var['position']} | 基因:{var['gene_coordinate']} | {var['raw_line']}\n")
        
#         self.logger.info(f"✅ TXT文件已保存 | TXT file saved: {output_file}")
    
#     def write_excel_output(self, exonic_variants: List[Dict], all_variants: List[Dict],
#                           output_file: Path, region_info: Dict):
#         """写入Excel文件 | Write Excel file"""
#         with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
#             # Sheet 1: 基因信息摘要 | Gene information summary
#             summary_data = {
#                 '项目 | Item': [
#                     '基因ID | Gene ID',
#                     '染色体 | Chromosome',
#                     '链方向 | Strand',
#                     '基因起始位置 | Gene start',
#                     '基因终止位置 | Gene end',
#                     '基因长度(bp) | Gene length',
#                     '分析区域起始 | Region start',
#                     '分析区域终止 | Region end',
#                     '分析区域长度(bp) | Region length',
#                     '上下游扩展(bp) | Extension',
#                     '外显子变异数 | Exonic variants',
#                     '总变异数 | Total variants',
#                     '生成时间 | Generated time'
#                 ],
#                 '值 | Value': [
#                     region_info['gene_id'],
#                     region_info['chrom'],
#                     region_info['strand'],
#                     region_info['gene_start'],
#                     region_info['gene_end'],
#                     region_info['gene_end'] - region_info['gene_start'] + 1,
#                     region_info['region_start'],
#                     region_info['region_end'],
#                     region_info['region_end'] - region_info['region_start'] + 1,
#                     region_info['extend'],
#                     len(exonic_variants),
#                     len(all_variants),
#                     datetime.now().strftime('%Y-%m-%d %H:%M:%S')
#                 ]
#             }
#             df_summary = pd.DataFrame(summary_data)
#             df_summary.to_excel(writer, sheet_name='摘要信息 Summary', index=False)
            
#             # Sheet 2: 外显子变异 | Exonic variants
#             if exonic_variants:
#                 exonic_df = pd.DataFrame([{
#                     '序号 | No.': idx,
#                     '染色体 | Chr': var['chrom'],
#                     '基因组位置 | Genomic Pos': var['position'],
#                     '基因坐标 | Gene Coord': var['gene_coordinate'],
#                     '位置类型 | Location': '基因内' if var['location_type'] == 'gene' else ('上游' if var['location_type'] == 'upstream' else '下游'),
#                     '结束位置 | End': var['end_pos'],
#                     '参考碱基 | Ref': var['ref'],
#                     '变异碱基 | Alt': var['alt'],
#                     '变异类型 | Type': var['variant_type'],
#                     '基因信息 | Gene Info': var['gene_info'],
#                     '原始行号 | Line No.': var['line_number']
#                 } for idx, var in enumerate(exonic_variants, 1)])
#             else:
#                 exonic_df = pd.DataFrame(columns=['序号', '染色体', '基因组位置', '基因坐标', '位置类型', 
#                                                   '结束位置', '参考碱基', '变异碱基', '变异类型', '基因信息'])
#             exonic_df.to_excel(writer, sheet_name='外显子变异 Exonic', index=False)
            
#             # Sheet 3: 所有变异 | All variants
#             if all_variants:
#                 all_df = pd.DataFrame([{
#                     '序号 | No.': idx,
#                     '染色体 | Chr': var['chrom'],
#                     '基因组位置 | Genomic Pos': var['position'],
#                     '基因坐标 | Gene Coord': var['gene_coordinate'],
#                     '位置类型 | Location': '基因内' if var['location_type'] == 'gene' else ('上游' if var['location_type'] == 'upstream' else '下游'),
#                     '结束位置 | End': var['end_pos'],
#                     '参考碱基 | Ref': var['ref'],
#                     '变异碱基 | Alt': var['alt'],
#                     '功能区域 | Function': var['function'],
#                     '基因信息 | Gene Info': var['gene_info']
#                 } for idx, var in enumerate(all_variants, 1)])
#             else:
#                 all_df = pd.DataFrame(columns=['序号', '染色体', '基因组位置', '基因坐标', '位置类型',
#                                                '结束位置', '参考碱基', '变异碱基', '功能区域', '基因信息'])
#             all_df.to_excel(writer, sheet_name='所有变异 All Variants', index=False)
        
#         self.logger.info(f"✅ Excel文件已保存 | Excel file saved: {output_file}")


# ===== FILE: filter_annovar/writer.py =====
"""
输出写入模块 | Output Writer Module
"""

import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import List, Dict

class OutputWriter:
    """输出文件写入器 | Output File Writer"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def write_txt_output(self, variants: List[Dict], output_file: Path, 
                        region_info: Dict, variant_type: str = "变异"):
        """写入格式化的TXT文件 | Write formatted TXT file"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("=" * 120 + "\n")
            f.write(f"🧬 基因区域{variant_type}提取报告 | Gene Region {variant_type} Extraction Report\n")
            f.write("=" * 120 + "\n")
            f.write(f"⏰ 生成时间 | Generation time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"🎯 基因ID | Gene ID: {region_info['gene_id']}\n")
            f.write(f"📍 染色体 | Chromosome: {region_info['chrom']}\n")
            f.write(f"📏 基因位置 | Gene location: {region_info['gene_start']}-{region_info['gene_end']} (链方向|Strand: {region_info['strand']})\n")
            f.write(f"📏 基因长度 | Gene length: {region_info['gene_end'] - region_info['gene_start'] + 1} bp\n")
            f.write(f"📐 分析区域 | Analysis region: {region_info['region_start']}-{region_info['region_end']} (±{region_info['extend']}bp)\n")
            f.write(f"📐 区域总长度 | Total region length: {region_info['region_end'] - region_info['region_start'] + 1} bp\n")
            f.write(f"📊 变异总数 | Total variants: {len(variants)}\n")
            f.write("=" * 120 + "\n\n")
            
            f.write("📖 坐标说明 | Coordinate explanation:\n")
            f.write("  - 基因坐标: 相对于基因起始的1-based坐标 | Gene coordinate: 1-based relative to gene start\n")
            if region_info['strand'] == '+':
                f.write("  - 正链基因: 基因坐标从5'端(左侧)开始计数 | Forward strand: counted from 5' end (left)\n")
                f.write("    * 上游(左侧): 负数, 如-1000 = 上游1000bp | Upstream(left): negative, e.g., -1000 = 1000bp upstream\n")
                f.write("    * 下游(右侧): 正数, 如+500 = 下游500bp | Downstream(right): positive, e.g., +500 = 500bp downstream\n")
            else:
                f.write("  - 负链基因: 基因坐标从5'端(右侧)开始计数 | Reverse strand: counted from 5' end (right)\n")
                f.write("    * 上游(右侧): 负数, 如-1000 = 上游1000bp | Upstream(right): negative, e.g., -1000 = 1000bp upstream\n")
                f.write("    * 下游(左侧): 正数, 如+500 = 下游500bp | Downstream(left): positive, e.g., +500 = 500bp downstream\n")
            f.write("=" * 120 + "\n\n")
            
            if len(variants) == 0:
                f.write("ℹ️  未找到变异 | No variants found.\n")
                return
            
            # 写入变异表格 | Write variant table
            f.write(f"{'序号':<6} {'基因组位置':<12} {'基因坐标':<12} {'位置类型':<10} {'参考':<8} {'变异':<8} {'类型/功能':<20} 基因信息\n")
            f.write("-" * 120 + "\n")
            
            for idx, var in enumerate(variants, 1):
                pos = var['position']
                gene_coord = var['gene_coordinate']
                location = var['location_type']
                ref = var['ref']
                alt = var['alt']
                
                # 位置类型显示 | Location type display
                if location == 'gene':
                    loc_display = '基因内'
                elif location == 'upstream':
                    loc_display = '上游'
                else:
                    loc_display = '下游'
                
                # 智能获取类型信息
                if 'variant_type' in var:
                    type_info = var['variant_type']
                elif 'function' in var:
                    type_info = var['function']
                else:
                    type_info = 'unknown'
                
                gene_info = var['gene_info'][:25] + '...' if len(var['gene_info']) > 25 else var['gene_info']
                
                f.write(f"{idx:<6} {pos:<12} {gene_coord:<12} {loc_display:<10} {ref:<8} {alt:<8} {type_info:<20} {gene_info}\n")
            
            f.write("\n" + "=" * 120 + "\n")
            f.write("📝 详细信息（完整行）| Detailed information (complete lines)\n")
            f.write("=" * 120 + "\n\n")
            
            for idx, var in enumerate(variants, 1):
                f.write(f"[{idx}] 基因组:{var['position']} | 基因:{var['gene_coordinate']} | {var['raw_line']}\n")
        
        self.logger.info(f"✅ TXT文件已保存 | TXT file saved: {output_file}")
    
    def write_excel_output(self, exonic_variants: List[Dict], all_variants: List[Dict],
                          output_file: Path, region_info: Dict):
        """写入Excel文件 | Write Excel file"""
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Sheet 1: 基因信息摘要 | Gene information summary
            summary_data = {
                '项目 | Item': [
                    '基因ID | Gene ID',
                    '染色体 | Chromosome',
                    '链方向 | Strand',
                    '基因起始位置 | Gene start',
                    '基因终止位置 | Gene end',
                    '基因长度(bp) | Gene length',
                    '分析区域起始 | Region start',
                    '分析区域终止 | Region end',
                    '分析区域长度(bp) | Region length',
                    '上下游扩展(bp) | Extension',
                    '外显子变异数 | Exonic variants',
                    '总变异数 | Total variants',
                    '生成时间 | Generated time'
                ],
                '值 | Value': [
                    region_info['gene_id'],
                    region_info['chrom'],
                    region_info['strand'],
                    region_info['gene_start'],
                    region_info['gene_end'],
                    region_info['gene_end'] - region_info['gene_start'] + 1,
                    region_info['region_start'],
                    region_info['region_end'],
                    region_info['region_end'] - region_info['region_start'] + 1,
                    region_info['extend'],
                    len(exonic_variants),
                    len(all_variants),
                    datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                ]
            }
            df_summary = pd.DataFrame(summary_data)
            df_summary.to_excel(writer, sheet_name='摘要信息 Summary', index=False)
            
            # Sheet 2: 外显子变异 | Exonic variants
            if exonic_variants:
                exonic_df = pd.DataFrame([{
                    '序号 | No.': idx,
                    '染色体 | Chr': var['chrom'],
                    '基因组位置 | Genomic Pos': var['position'],
                    '基因坐标 | Gene Coord': var['gene_coordinate'],
                    '位置类型 | Location': '基因内' if var['location_type'] == 'gene' else ('上游' if var['location_type'] == 'upstream' else '下游'),
                    '结束位置 | End': var['end_pos'],
                    '参考碱基 | Ref': var['ref'],
                    '变异碱基 | Alt': var['alt'],
                    '变异类型 | Type': var.get('variant_type', 'unknown'),
                    '基因信息 | Gene Info': var['gene_info'],
                    '原始行号 | Line No.': var.get('line_number', '')
                } for idx, var in enumerate(exonic_variants, 1)])
            else:
                exonic_df = pd.DataFrame(columns=['序号', '染色体', '基因组位置', '基因坐标', '位置类型', 
                                                  '结束位置', '参考碱基', '变异碱基', '变异类型', '基因信息'])
            exonic_df.to_excel(writer, sheet_name='外显子变异 Exonic', index=False)
            
            # Sheet 3: 所有变异 | All variants
            if all_variants:
                all_df = pd.DataFrame([{
                    '序号 | No.': idx,
                    '染色体 | Chr': var['chrom'],
                    '基因组位置 | Genomic Pos': var['position'],
                    '基因坐标 | Gene Coord': var['gene_coordinate'],
                    '位置类型 | Location': '基因内' if var['location_type'] == 'gene' else ('上游' if var['location_type'] == 'upstream' else '下游'),
                    '结束位置 | End': var['end_pos'],
                    '参考碱基 | Ref': var['ref'],
                    '变异碱基 | Alt': var['alt'],
                    '功能区域 | Function': var.get('function', var.get('variant_type', 'unknown')),
                    '基因信息 | Gene Info': var['gene_info']
                } for idx, var in enumerate(all_variants, 1)])
            else:
                all_df = pd.DataFrame(columns=['序号', '染色体', '基因组位置', '基因坐标', '位置类型',
                                               '结束位置', '参考碱基', '变异碱基', '功能区域', '基因信息'])
            all_df.to_excel(writer, sheet_name='所有变异 All Variants', index=False)
        
        self.logger.info(f"✅ Excel文件已保存 | Excel file saved: {output_file}")

# ===== END FILE =====
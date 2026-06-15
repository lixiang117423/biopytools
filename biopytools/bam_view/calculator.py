"""
BAM比对可视化核心计算模块|BAM Alignment Visualization Core Calculation Module
"""

import os
import subprocess
from .config import BamViewConfig
from .utils import format_number, build_conda_command


class BamViewCalculator:
    """BAM比对可视化计算器|BAM Alignment Visualization Calculator"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def generate_visualization(self):
        """生成比对可视化|Generate alignment visualization"""
        self.logger.info("开始生成BAM比对可视化|Starting BAM alignment visualization generation")

        # 构建alignoth命令|Build alignoth command
        cmd = self._build_alignoth_command()

        # 根据输出格式执行不同处理|Execute different processing based on output format
        if self.config.output_format == 'html':
            success = self._generate_html_output(cmd)
        elif self.config.output_format == 'json':
            success = self._generate_json_output(cmd)
        elif self.config.output_format == 'svg':
            success = self._generate_svg_output(cmd)
        elif self.config.output_format == 'pdf':
            success = self._generate_pdf_output(cmd)
        else:
            self.logger.error(f"不支持的输出格式|Unsupported output format: {self.config.output_format}")
            return False

        if success:
            self.logger.info("BAM比对可视化生成完成|BAM alignment visualization generation completed")
        else:
            self.logger.error("BAM比对可视化生成失败|BAM alignment visualization generation failed")

        return success

    def _build_alignoth_command(self):
        """构建alignoth命令|Build alignoth command"""
        # 构建参数列表|Build arguments list
        args = [
            '-b', self.config.bam_file,
            '-r', self.config.reference,
            '-g', self.config.region,
            '-d', str(self.config.max_read_depth),
            '-w', str(self.config.max_width),
            '--mismatch-display-min-percent', str(self.config.mismatch_display_min_percent)
        ]

        # 高亮选项|Highlight options
        if self.config.vcf_file:
            args.extend(['-v', self.config.vcf_file])

        if self.config.bed_file:
            args.extend(['--bed', self.config.bed_file])

        if self.config.highlight_intervals:
            for interval in self.config.highlight_intervals:
                args.extend(['-h', interval])

        # 辅助标签|Aux tags
        if self.config.aux_tags:
            for tag in self.config.aux_tags:
                args.extend(['-x', tag])

        # 输出选项|Output options
        if self.config.plot_all:
            args.append('--plot-all')

        # HTML输出选项|HTML output options
        if self.config.output_format == 'html':
            args.append('--html')
            if self.config.no_embed_js:
                args.append('--no-embed-js')

        # 使用conda wrapper|Use conda wrapper
        cmd = build_conda_command(self.config.alignoth_path, args)
        return cmd

    def _generate_html_output(self, cmd: list) -> bool:
        """生成HTML输出|Generate HTML output"""
        self.logger.info("生成HTML格式输出|Generating HTML format output")

        # 生成输出文件名|Generate output filename
        bam_basename = os.path.basename(self.config.bam_file).replace('.bam', '').replace('.sam', '').replace('.cram', '')
        output_file = os.path.join(self.config.output_dir, f"{bam_basename}.html")

        # 执行命令并保存输出|Execute command and save output
        success, stdout = self.cmd_runner.run_with_output(cmd, "生成HTML可视化|Generate HTML visualization")

        if success:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(stdout)
            self.logger.info(f"HTML文件已保存|HTML file saved: {output_file}")
            return True
        else:
            return False

    def _generate_json_output(self, cmd: list) -> bool:
        """生成JSON输出|Generate JSON output"""
        self.logger.info("生成Vega-Lite JSON格式输出|Generating Vega-Lite JSON format output")

        # 生成输出文件名|Generate output filename
        bam_basename = os.path.basename(self.config.bam_file).replace('.bam', '').replace('.sam', '').replace('.cram', '')
        output_file = os.path.join(self.config.output_dir, f"{bam_basename}.vl.json")

        # 执行命令并保存输出|Execute command and save output
        success, stdout = self.cmd_runner.run_with_output(cmd, "生成JSON可视化|Generate JSON visualization")

        if success:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(stdout)
            self.logger.info(f"JSON文件已保存|JSON file saved: {output_file}")
            return True
        else:
            return False

    def _generate_svg_output(self, cmd: list) -> bool:
        """生成SVG输出|Generate SVG output"""
        self.logger.info("生成SVG格式输出|Generating SVG format output")

        # 生成输出文件名|Generate output filename
        bam_basename = os.path.basename(self.config.bam_file).replace('.bam', '').replace('.sam', '').replace('.cram', '')
        output_file = os.path.join(self.config.output_dir, f"{bam_basename}.svg")

        # 检查vl2vg是否可用|Check if vl2vg is available
        check_cmd = ["which", "vl2vg"]
        success, _ = self.cmd_runner.run_with_output(check_cmd, "检查vl2vg|Check vl2vg")

        if not success:
            self.logger.error(
                "vl2vg未找到，请先安装vega-cli|"
                "vl2vg not found, please install vega-cli first. "
                "npm install -g vega-cli vega-lite-cli"
            )
            return False

        # 执行命令，使用subprocess.Popen管道|Execute command using subprocess.Popen pipe
        try:
            self.logger.info(f"命令|Command: {' '.join(cmd)} | vl2vg > {output_file}")

            # 第一个进程：运行alignoth|First process: run alignoth
            p1 = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.config.output_dir
            )

            # 第二个进程：vl2vg|Second process: vl2vg
            p2 = subprocess.Popen(
                ["vl2vg"],
                stdin=p1.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.config.output_dir
            )

            # 关闭p1的stdout，允许p1接收SIGPIPE如果p2退出|Close p1.stdout to allow p1 to receive SIGPIPE if p2 exits
            p1.stdout.close()

            # 写入文件|Write to file
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(p2.stdout.read().decode('utf-8'))

            # 等待进程完成|Wait for processes to complete
            p1.wait()
            p2.wait()

            if p1.returncode == 0 and p2.returncode == 0:
                self.logger.info(f"SVG文件已保存|SVG file saved: {output_file}")
                return True
            else:
                self.logger.error(f"alignoth执行失败|alignoth execution failed, return code: {p1.returncode}")
                if p1.stderr:
                    self.logger.error(f"错误信息|Error: {p1.stderr.read().decode('utf-8')}")
                return False

        except Exception as e:
            self.logger.error(f"SVG生成失败|SVG generation failed: {e}")
            return False

    def _generate_pdf_output(self, cmd: list) -> bool:
        """生成PDF输出|Generate PDF output"""
        self.logger.info("生成PDF格式输出|Generating PDF format output")

        # 生成输出文件名|Generate output filename
        bam_basename = os.path.basename(self.config.bam_file).replace('.bam', '').replace('.sam', '').replace('.cram', '')
        output_file = os.path.join(self.config.output_dir, f"{bam_basename}.pdf")

        # 检查vl2vg和vg2pdf是否可用|Check if vl2vg and vg2pdf are available
        check_cmd = ["which", "vl2vg", "vg2pdf"]
        success, _ = self.cmd_runner.run_with_output(check_cmd, "检查vega-cli|Check vega-cli")

        if not success:
            self.logger.error(
                "vega-cli未找到，请先安装vega-cli|"
                "vega-cli not found, please install vega-cli first. "
                "npm install -g vega-cli vega-lite-cli"
            )
            return False

        # 执行命令，使用subprocess.Popen管道|Execute command using subprocess.Popen pipe
        try:
            self.logger.info(f"命令|Command: {' '.join(cmd)} | vl2vg | vg2pdf > {output_file}")

            # 第一个进程：运行alignoth|First process: run alignoth
            p1 = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.config.output_dir
            )

            # 第二个进程：vl2vg|Second process: vl2vg
            p2 = subprocess.Popen(
                ["vl2vg"],
                stdin=p1.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.config.output_dir
            )

            # 第三个进程：vg2pdf|Third process: vg2pdf
            p3 = subprocess.Popen(
                ["vg2pdf"],
                stdin=p2.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.config.output_dir
            )

            # 关闭p2的stdout|Close p2.stdout
            p2.stdout.close()

            # 写入文件|Write to file
            with open(output_file, 'wb') as f:
                f.write(p3.stdout.read())

            # 等待进程完成|Wait for processes to complete
            p1.wait()
            p2.wait()
            p3.wait()

            if p1.returncode == 0 and p2.returncode == 0 and p3.returncode == 0:
                self.logger.info(f"PDF文件已保存|PDF file saved: {output_file}")
                return True
            else:
                self.logger.error(f"alignoth执行失败|alignoth execution failed, return code: {p1.returncode}")
                if p1.stderr:
                    self.logger.error(f"错误信息|Error: {p1.stderr.read().decode('utf-8')}")
                return False

        except Exception as e:
            self.logger.error(f"PDF生成失败|PDF generation failed: {e}")
            return False

    def get_statistics(self):
        """获取统计信息|Get statistics"""
        # TODO: 可以添加更多统计信息，如reads数量、覆盖度等
        # TODO: Can add more statistics, such as read count, coverage, etc.
        stats = {
            'bam_file': self.config.bam_file,
            'reference': self.config.reference,
            'region': self.config.region,
            'output_format': self.config.output_format,
            'max_read_depth': self.config.max_read_depth,
            'max_width': self.config.max_width
        }
        return stats

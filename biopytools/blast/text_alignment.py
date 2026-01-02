"""
Text Alignment Generator
"""

from pathlib import Path
from typing import Dict, List

class TextAlignmentGenerator:
    """Generate text format alignment visualizations"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.output_dir = config.output_path / config.alignment_output_dir / "text"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_alignments(self, alignments_data: Dict):
        """Generate text alignment files for all samples"""
        self.logger.info("Generating text format alignment visualizations...")

        sample_files = []
        all_alignments = []

        for sample_name, sample_data in alignments_data.items():
            # Write individual sample file
            sample_file = self._write_sample_file(sample_name, sample_data)
            sample_files.append(sample_file)
            all_alignments.extend(sample_data['alignments'])

        # Write summary file
        summary_file = self._write_summary_file(alignments_data, all_alignments)

        self.logger.info(f"Text output directory: {self.output_dir}")
        return sample_files, summary_file

    def _write_sample_file(self, sample_name: str, sample_data: Dict) -> str:
        """Write alignment file for a single sample"""
        output_file = self.output_dir / f"{sample_name}_alignments.txt"
        alignments = sample_data['alignments']

        with open(output_file, 'w', encoding='utf-8') as f:
            # Header
            f.write("=" * 80 + "\n")
            f.write("BLAST Alignment Details\n")
            f.write("=" * 80 + "\n")
            f.write(f"Sample Name: {sample_name}\n")
            f.write(f"Input File: {sample_data['file_name']}\n")
            f.write(f"Alignment Count: {len(alignments)}\n")

            if alignments:
                avg_identity = sum(a['identity'] for a in alignments) / len(alignments)
                avg_coverage = sum(a['coverage'] for a in alignments) / len(alignments)
                f.write(f"Average Identity: {avg_identity:.2f}%\n")
                f.write(f"Average Coverage: {avg_coverage:.2f}%\n")

            f.write("=" * 80 + "\n\n")

            # Alignments
            for idx, alignment in enumerate(alignments, 1):
                alignment_text = self._format_alignment(idx, alignment)
                f.write(alignment_text)
                f.write("\n")

        self.logger.info(f"  Sample {sample_name}: {len(alignments)} alignments")
        return str(output_file)

    def _write_summary_file(self, alignments_data: Dict, all_alignments: List) -> str:
        """Write summary file for all samples"""
        output_file = self.output_dir / "all_samples_alignments.txt"

        with open(output_file, 'w', encoding='utf-8') as f:
            # Header
            f.write("=" * 80 + "\n")
            f.write("BLAST Alignment Summary - All Samples\n")
            f.write("=" * 80 + "\n")
            f.write(f"Total Samples: {len(alignments_data)}\n")
            f.write(f"Total Alignments: {len(all_alignments)}\n")

            if all_alignments:
                avg_identity = sum(a['identity'] for a in all_alignments) / len(all_alignments)
                avg_coverage = sum(a['coverage'] for a in all_alignments) / len(all_alignments)
                f.write(f"Average Identity: {avg_identity:.2f}%\n")
                f.write(f"Average Coverage: {avg_coverage:.2f}%\n")

            f.write("=" * 80 + "\n\n")

            # Sample details
            for sample_name, sample_data in alignments_data.items():
                f.write("\n" + "=" * 80 + "\n")
                f.write(f"Sample: {sample_name}\n")
                f.write("=" * 80 + "\n\n")

                for idx, alignment in enumerate(sample_data['alignments'], 1):
                    alignment_text = self._format_alignment(idx, alignment)
                    f.write(alignment_text)
                    f.write("\n")

        self.logger.info(f"  Total samples: {len(alignments_data)}")
        return str(output_file)

    def _format_alignment(self, index: int, alignment: Dict) -> str:
        """Format a single alignment as text"""
        lines = []

        # Header
        lines.append("=" * 80)
        lines.append(f"Alignment #{index}: {alignment['query_id']} -> {alignment['subject_id']}")
        lines.append("=" * 80)
        lines.append(f"Identity: {alignment['identity']:.2f}% | "
                    f"Coverage: {alignment['coverage']:.2f}% | "
                    f"E-value: {alignment['evalue']} | "
                    f"Bit Score: {alignment['bitscore']}")
        lines.append(f"Length: {alignment['length']} | "
                    f"Mismatches: {alignment['mismatch']} | "
                    f"Gap Opens: {alignment['gapopen']}")
        lines.append("=" * 80)
        lines.append("")

        # Sequence alignment
        query_seq = alignment.get('query_seq', '')
        subject_seq = alignment.get('subject_seq', '')

        if not query_seq or not subject_seq:
            # No sequence data available
            lines.append("  Sequence alignment not available")
            lines.append("  Note: Re-run BLAST with qseq and sseq in output format")
            lines.append("")
        else:
            # Generate match line
            match_line = self._generate_match_line(query_seq, subject_seq)

            # Format sequence alignment
            width = self.config.alignment_width
            q_start = alignment['qstart']
            s_start = alignment['sstart']

            # Calculate positions (excluding gaps)
            q_pos = q_start
            s_pos = s_start

            for i in range(0, len(query_seq), width):
                q_segment = query_seq[i:i+width]
                s_segment = subject_seq[i:i+width]
                m_segment = match_line[i:i+width]

                # Count non-gap characters in this segment
                q_bases = sum(1 for c in q_segment if c != '-')
                s_bases = sum(1 for c in s_segment if c != '-')

                # Calculate end positions for this segment
                q_end_pos = q_pos + q_bases - 1 if q_bases > 0 else q_pos
                s_end_pos = s_pos + s_bases - 1 if s_bases > 0 else s_pos

                # Format output (using 10-character width for alignment)
                lines.append(f"Query  {q_pos:10d}  {q_segment}  {q_end_pos}")
                lines.append(f"                   {m_segment}")
                lines.append(f"Sbjct  {s_pos:10d}  {s_segment}  {s_end_pos}")
                lines.append("")

                # Update positions (only count non-gap characters)
                q_pos += q_bases
                s_pos += s_bases

            # Statistics
            match_count = match_line.count('|')
            mismatch_count = match_line.count('.')
            gap_count = match_line.count(' ')
            total = len(match_line)

            lines.append("Alignment Statistics:")
            lines.append(f"  Matches: {match_count} / {total} ({match_count/total*100:.1f}%)")
            lines.append(f"  Mismatches: {mismatch_count} ({mismatch_count/total*100:.1f}%)")
            if gap_count > 0:
                lines.append(f"  Gaps: {gap_count} ({gap_count/total*100:.1f}%)")

        lines.append("")
        return "\n".join(lines)

    def _generate_match_line(self, query_seq: str, subject_seq: str) -> str:
        """Generate match line between query and subject sequences"""
        match_line = []
        for q, s in zip(query_seq, subject_seq):
            if q == s:
                match_line.append('|')
            elif q == '-' or s == '-':
                match_line.append(' ')
            else:
                match_line.append('.')
        return ''.join(match_line)

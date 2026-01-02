"""
BLAST Analysis Utility Functions Module
"""

import os
import subprocess
import shutil
from pathlib import Path
from typing import List, Optional


def check_tool_availability(tool_path: str) -> bool:
    """Check if a tool is available in the system"""
    return shutil.which(tool_path) is not None


def get_input_files(input_path: str, suffix: str = "*.fa") -> List[str]:
    """
    Get list of input files

    Args:
        input_path: Input file or directory path
        suffix: File suffix pattern for directory scanning

    Returns:
        List of file paths
    """
    input_path_obj = Path(input_path)

    # If it's a single file, return it directly
    if input_path_obj.is_file():
        return [str(input_path_obj)]

    # If it's a directory, scan for files
    if input_path_obj.is_dir():
        if suffix.startswith('*'):
            pattern = suffix[1:]
            files = list(input_path_obj.glob(f"*{pattern}"))
        else:
            files = list(input_path_obj.glob(suffix))
        return [str(f) for f in files if f.is_file()]

    return []


def extract_sample_name_from_path(file_path: str, pattern: str = r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$') -> str:
    """
    Extract sample name from file path

    Args:
        file_path: File path
        pattern: Regex pattern for extracting sample name

    Returns:
        Sample name
    """
    import re
    filename = os.path.basename(file_path)

    match = re.search(pattern, filename)
    if match:
        return match.group(1)

    return Path(filename).stem


def create_sample_map_file(input_files: List[str], output_file: str, pattern: str) -> int:
    """
    Create sample mapping file

    Args:
        input_files: List of input file paths
        output_file: Output sample map file path
        pattern: Sample name extraction pattern

    Returns:
        Number of samples in the map file
    """
    from datetime import datetime

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("# Auto-generated sample mapping file\n")
        f.write(f"# Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# Format: file_path<TAB>sample_name\n")
        f.write("#" + "="*60 + "\n")

        sample_count = 0
        for input_file in sorted(input_files):
            sample_name = extract_sample_name_from_path(input_file, pattern)
            f.write(f"{input_file}\t{sample_name}\n")
            sample_count += 1

    return sample_count


def load_sample_mapping(map_file: str) -> dict:
    """
    Load sample mapping from file

    Args:
        map_file: Sample mapping file path

    Returns:
        Dictionary mapping filename to sample name
    """
    sample_mapping = {}
    file_paths = []

    if not os.path.exists(map_file):
        raise FileNotFoundError(f"Sample mapping file not found: {map_file}")

    with open(map_file, 'r', encoding='utf-8') as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                file_path = parts[0]
                sample_name = parts[1]

                if not os.path.exists(file_path):
                    continue

                filename = os.path.basename(file_path)
                sample_mapping[filename] = sample_name
                file_paths.append(file_path)

    return {
        'mapping': sample_mapping,
        'file_paths': file_paths
    }


def format_command(cmd: list) -> str:
    """Format command list to string"""
    return ' '.join(str(x) for x in cmd)


def run_command(cmd: list, logger=None, description: str = "") -> bool:
    """
    Run a command

    Args:
        cmd: Command as list of arguments
        logger: Logger object
        description: Command description

    Returns:
        True if successful, False otherwise
    """
    if logger and description:
        logger.info(f"Executing: {description}")

    if logger:
        logger.debug(f"Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        if logger and description:
            logger.info(f"Command completed successfully: {description}")
        return True

    except subprocess.CalledProcessError as e:
        if logger:
            logger.error(f"Command failed: {description}")
            logger.error(f"Return code: {e.returncode}")
            if e.stderr:
                logger.error(f"Error output: {e.stderr}")
        return False


def check_blast_dependencies(config, logger) -> bool:
    """
    Check BLAST dependencies

    Args:
        config: BLAST configuration object
        logger: Logger object

    Returns:
        True if all dependencies available

    Raises:
        RuntimeError: If any dependency is missing
    """
    logger.info("Checking BLAST dependencies")

    # Get the correct BLAST program path
    blast_program = getattr(config, f'{config.blast_type}_path', config.blast_type)

    dependencies = [
        (config.makeblastdb_path, "makeblastdb"),
        (blast_program, config.blast_type.upper())
    ]

    missing_deps = []

    for cmd, name in dependencies:
        if check_tool_availability(cmd):
            logger.info(f"  {name}: available")
        else:
            logger.warning(f"  {name}: NOT found")
            missing_deps.append(name)

    if missing_deps:
        error_msg = f"Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    logger.info("All dependencies available")
    return True

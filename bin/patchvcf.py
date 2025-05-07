#!/usr/bin/env python3
# patchvcf.py - A script to process and modify VCF files using a mapping file

import argparse
import sys
import logging
import os
import gzip
from typing import List, Tuple
import pandas as pd

def parse_args(argv=None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Process and modify VCF files using a mapping file")

    parser.add_argument("-m", "--map", required=True, help="Path to the mapping file")
    parser.add_argument("-v", "--vcf", required=True, help="Path to the input VCF file")
    parser.add_argument("-f", "--fasta", help="Path to the reference FASTA file (optional)")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output files")
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)

def setup_logging(log_level: str, prefix: str) -> logging.Logger:
    """Set up logging configuration."""
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    logger = logging.getLogger('patchvcf')
    logger.setLevel(log_level)

    # Create file handler
    file_handler = logging.FileHandler(f"{prefix}.log")
    file_handler.setLevel(log_level)
    file_handler.setFormatter(logging.Formatter(log_format))

    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(logging.Formatter(log_format))

    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

def read_map_file(file_path: str, logger: logging.Logger) -> pd.DataFrame:
    """Read and parse the mapping file into a pandas DataFrame.

    Args:
        file_path: Path to the mapping file
        logger: Logger object

    Returns:
        A pandas DataFrame containing the mapping information with columns:
        - genome: The genome/reference name
        - old_position: Position in the original sequence
        - nucleotide: The nucleotide at this position
        - new_position: Position in the reference alignment
    """
    logger.info(f"Reading mapping file: {file_path}")

    # Lists to store the parsed data
    genomes = []
    old_positions = []
    nucleotides = []
    new_positions = []

    current_genome = None

    try:
        with open(file_path, 'r') as map_file:
            for line in map_file:
                line = line.strip()

                # Handle genome header lines
                if line.startswith('>'):
                    # Extract the genome ID (everything before the first space)
                    current_genome = line[1:].split()[0]
                    logger.debug(f"Processing mappings for genome: {current_genome}")
                    continue

                # Skip comment lines
                if line.startswith('#'):
                    continue

                # Process mapping entry
                if current_genome and line:
                    try:
                        # Parse "nucleotide, old_position, new_position"
                        parts = line.split(',')
                        if len(parts) == 3:
                            nucleotide = parts[0].strip().lower()
                            old_position = int(parts[1].strip())
                            new_position = int(parts[2].strip())

                            genomes.append(current_genome)
                            old_positions.append(old_position)
                            nucleotides.append(nucleotide)
                            new_positions.append(new_position)
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Skipping invalid line in map file: {line}, error: {e}")

        # Create DataFrame from the collected data
        map_df = pd.DataFrame({
            'genome': genomes,
            'old_position': old_positions,
            'nucleotide': nucleotides,
            'new_position': new_positions
        })

        logger.info(f"Map file loaded successfully with {len(map_df)} entries across {map_df['genome'].nunique()} genomes")
        return map_df
    except Exception as e:
        logger.error(f"Failed to read map file: {e}")
        sys.exit(1)

def read_vcf_file(file_path: str, logger: logging.Logger, fasta_path: str =None ) -> Tuple[List[str], pd.DataFrame]:
    """Read and parse the VCF file."""
    logger.info(f"Reading VCF file: {file_path}")
    header_lines = []
    try:
        # Open the file with the appropriate method based on extension
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = 'rt' if file_path.endswith('.gz') else 'r'

        # Read header lines
        with open_func(file_path, mode) as vcf_file:
            for line in vcf_file:
                line = line.strip()
                if line.startswith('#'):
                    header_lines.append(line)
                else:
                    break  # Stop reading at the first non-header line

        names = header_lines[-1].split("\t") if header_lines else None
        header = None if names else "infer"
        header_lines = header_lines[:-1] if names else header_lines
        vcf_data = pd.read_csv(file_path, sep="\t", compression="infer", comment="#", header=header, names=names)

        logger.debug(f"VCF file loaded with {len(header_lines)} header lines and {len(vcf_data)} data lines")

        # Check for required columns depending on custom VCF format or standard VCF
        if '#CHROM' not in vcf_data.columns:
            logger.warning("No '#CHROM' column found. This might not be a standard VCF file.")
            if fasta_path:
                with open(fasta_path, 'r') as fasta_file:
                    chrom = fasta_file.readline().strip().split()[0][1:]
                    logger.debug(f"Using chromosome name from FASTA file: {chrom}")
                    vcf_data['#CHROM'] = chrom
            else:
                logger.error("No '#CHROM' column found and no FASTA file provided to infer chromosome names.")
                sys.exit(1)

        if 'POS' not in vcf_data.columns:
            if 'Position' in vcf_data.columns:
                vcf_data.rename(columns={'Position': 'POS'}, inplace=True)
                logger.debug("Renamed 'Position' to 'POS'")
            else:
                logger.error("No 'POS' or 'Position' column found in VCF data")
                sys.exit(1)
        if 'REF' not in vcf_data.columns:
            if 'Reference' in vcf_data.columns:
                vcf_data.rename(columns={'Reference': 'REF'}, inplace=True)
                logger.debug("Renamed 'Reference' to 'REF'")
            else:
                logger.error("No 'REF' or 'Reference' column found in VCF data")
                sys.exit(1)

        return header_lines, vcf_data

    except Exception as e:
        logger.error(f"Failed to read VCF file: {e}")
        sys.exit(1)

def process_vcf_with_map(vcf_header: List[str], vcf_data: pd.DataFrame,
                         map_data: pd.DataFrame, logger: logging.Logger) -> Tuple[List[str], pd.DataFrame]:
    """Process VCF data using the mapping information.

    Args:
        vcf_header: List of VCF header lines
        vcf_data: DataFrame containing the VCF data
        map_data: DataFrame containing mapping information
        logger: Logger object

    Returns:
        Tuple containing the VCF header lines and the processed VCF data
    """
    logger.info("Processing VCF data with mapping information")

    # Make a copy to avoid modifying the original data
    processed_data = vcf_data.copy()

    # Extract reference information from VCF
    unique_refs = processed_data['#CHROM'].unique()
    logger.debug(f"Found {len(unique_refs)} unique references in VCF: {', '.join(unique_refs)}")

    # Check if references in VCF exist in the mapping data
    missing_refs = [ref for ref in unique_refs if ref not in map_data['genome'].unique()]
    if missing_refs:
        logger.warning(f"References in VCF not found in mapping data: {', '.join(missing_refs)}")

    # Convert POS to int if it's not already
    if processed_data['POS'].dtype != 'int64':
        processed_data['POS'] = processed_data['POS'].astype(int)

    # Create a merge key by combining genome and position
    processed_data['merge_key'] = list(zip(processed_data['#CHROM'], processed_data['POS']))
    map_data['merge_key'] = list(zip(map_data['genome'], map_data['old_position']))

    # Join the dataframes on the merge key
    merged_data = pd.merge(processed_data,
                           map_data[['merge_key', 'new_position', 'nucleotide']],
                           on='merge_key',
                           how='left')

    # Count mapped and unmapped positions
    mapped = merged_data['new_position'].notna()
    mapped_count = mapped.sum()
    unmapped_count = (~mapped).sum()

    # Check for reference nucleotide mismatches
    ref_mismatches = 0
    if mapped_count > 0:
        # Create lowercase versions for comparison
        merged_data['ref_lower'] = merged_data['REF'].str.lower()
        merged_data['nucleotide_lower'] = merged_data['nucleotide'].str.lower()

        # Identify mismatches
        mismatches = (merged_data['new_position'].notna() &
                      (merged_data['ref_lower'] != merged_data['nucleotide_lower']))
        ref_mismatches = mismatches.sum()

        if ref_mismatches > 0:
            logger.warning(f"Found {ref_mismatches} positions with reference mismatches")
            # For mismatches, we'll leave the position as NA
            merged_data.loc[mismatches, 'new_position'] = None

    # Update positions where there's a valid mapping
    merged_data.loc[merged_data['new_position'].notna(), 'POS'] = merged_data.loc[merged_data['new_position'].notna(), 'new_position']

    # Clean up the dataframe
    result_data = merged_data.drop(columns=['merge_key', 'new_position', 'nucleotide',
                                           'ref_lower', 'nucleotide_lower'], errors='ignore')

    logger.info(f"Processed {len(processed_data)} VCF records: "
               f"{mapped_count - ref_mismatches} positions mapped, {unmapped_count + ref_mismatches} positions unchanged")

    return vcf_header, result_data

def write_output_vcf(vcf_header: List[str], vcf_data: pd.DataFrame, prefix: str, logger: logging.Logger) -> str:
    """Write the processed VCF data to an output file.

    Args:
        vcf_header: List of VCF header lines
        vcf_data: DataFrame containing the processed VCF data
        prefix: Prefix for the output file
        logger: Logger object

    Returns:
        Path to the output file
    """
    output_file = f"{prefix}.vcf"
    logger.info(f"Writing processed data to {output_file}")

    try:
        # Write header lines
        with open(output_file, 'w') as out_file:
            for line in vcf_header:
                out_file.write(f"{line}\n")

        # Append data without writing pandas index
        vcf_data.to_csv(output_file, sep='\t', index=False, mode='a')

        logger.info(f"Successfully wrote VCF with {len(vcf_header)} header lines and {len(vcf_data)} data rows to {output_file}")
        return output_file
    except Exception as e:
        logger.error(f"Failed to write output VCF file: {e}")
        sys.exit(1)

def main(argv=None):
    """Main function to orchestrate the VCF processing workflow."""
    # Parse command line arguments
    args = parse_args(argv)
    # Set up logging
    logger = setup_logging(args.log_level, args.prefix)
    logger.info("Starting VCF patching process")

    # Process files
    try:
        # Read input files
        map_data = read_map_file(args.map, logger)
        vcf_header, vcf_data = read_vcf_file(args.vcf, logger, args.fasta)

        # Process VCF with mapping data
        processed_header, processed_data = process_vcf_with_map(vcf_header, vcf_data, map_data, logger)

        # Write output
        output_file = write_output_vcf(processed_header, processed_data, args.prefix, logger)

        logger.info(f"VCF patching completed successfully. Output written to {output_file}")
        return 0
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())

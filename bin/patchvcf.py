#!/usr/bin/env python3
# patchvcf.py - A script to process and modify VCF files using a mapping file

import argparse
import sys
import logging
import os
import gzip
from typing import List, Tuple, Dict
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
    log_format = '%(asctime)s - %(levelname)s - %(message)s'

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

def reverse_complement(sequence: str) -> str:
    """Generate the reverse complement of a nucleotide sequence.

    Args:
        sequence: Input nucleotide sequence (string of A, T, G, C)

    Returns:
        Reverse complemented sequence
    """
    complement_dict = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                       'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                       'n': 'n', 'N': 'N'}

    # Reverse the sequence and get the complement of each nucleotide
    rev_comp = ''.join(complement_dict.get(base, base) for base in reversed(sequence))
    return rev_comp

def read_map_file(file_path: str, logger: logging.Logger) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Read and parse the mapping file into a pandas DataFrame.

    Args:
        file_path: Path to the mapping file
        logger: Logger object

    Returns:
        A tuple containing:
        - A pandas DataFrame with the mapping information (columns: genome, old_position, nucleotide, new_position)
        - A dictionary mapping reverse complemented genome names to their maximum positions
    """
    logger.info(f"Reading mapping file: {file_path}")

    # Lists to store the parsed data
    data_rows = []

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

                    # Check if this is a reverse complemented genome
                    if current_genome.startswith('_R_'):
                        logger.info(f"Detected reverse complemented genome: {current_genome}")

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
                            data_rows.append({
                                'genome': current_genome,
                                'nucleotide': parts[0].strip().lower(),
                                'old_position': int(parts[1].strip()),
                                'new_position': int(parts[2].strip())
                            })

                    except (ValueError, IndexError) as e:
                        logger.warning(f"Skipping invalid line in map file: {line}, error: {e}")

        # Create DataFrame from the collected data
        map_df = pd.DataFrame(data_rows)
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
                         map_data: pd.DataFrame,logger: logging.Logger) -> Tuple[List[str], pd.DataFrame]:
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
                           map_data[['merge_key', 'new_position', 'nucleotide', 'is_RC']],
                           on='merge_key',
                           how='left')

    merged_data = patch_nucleotide_rc(merged_data,logger)

    # Count mapped and unmapped positions
    mapped = merged_data['new_position'].notna()
    mapped_count = mapped.sum()
    unmapped_count = (~mapped).sum()

    # Check for reference nucleotide mismatches
    if mapped_count > 0:
        # Create lowercase versions for comparison
        merged_data['ref_lower'] = merged_data['REF'].str.lower()
        merged_data['nucleotide_lower'] = merged_data['nucleotide'].str.lower()

        # Identify mismatches
        mismatches = (merged_data['new_position'].notna() &
                      (merged_data['ref_lower'] != merged_data['nucleotide_lower']))

    # Update positions where there's a valid mapping
    merged_data.loc[merged_data['new_position'].notna(), 'POS'] = merged_data.loc[merged_data['new_position'].notna(), 'new_position']

    # Clean up the dataframe
    result_data = merged_data.drop(columns=['merge_key', 'new_position', 'nucleotide',
                                           'ref_lower', 'nucleotide_lower'], errors='ignore')

    logger.info(f"Processed {len(processed_data)} VCF records: "
               f"{mapped_count } positions mapped, {unmapped_count} positions unchanged")

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

def patch_map_rc(map_data: pd.DataFrame) -> pd.DataFrame:
    """Patch genome locations for reverse complemented genomes."""
    # Identify reverse complemented genomes
    df = map_data.copy()
    df['is_RC'] = df['genome'].str.startswith('_R_')

    df.to_csv("reverse_complement_genomes.tsv", index=False, sep="\t", header=True)

    # Remove the '_R_' prefix from genome names
    df['genome'] = df['genome'].str.replace('_R_', '', regex=False)

    # Identify max position
    df['max_position'] = df.groupby('genome')['old_position'].transform('max')

    # Update old positions for reverse complemented genomes
    df['old_position'] = df.apply(
        lambda x: x['max_position'] - x['old_position'] + 1 if x['is_RC'] else x['old_position'], axis=1
    )
    df.drop(columns=['max_position'], inplace=True)

    df.to_csv("updated_map_data.tsv", index=False, sep="\t", header=True)

    return df

def patch_nucleotide_rc(df: pd.DataFrame, logger) -> pd.DataFrame:
    """Patch nucleotide sequences for reverse complemented genomes."""
    # Create a copy to avoid modifying the original dataframe
    df = df.copy()

    logger.debug("%s",df.columns)
    # Create a boolean mask for rows where reverse complement should be applied
    rc_mask = df['is_RC'].fillna(False)

    # If there are no rows to reverse complement, return the dataframe as is
    if not rc_mask.any():
        return df

    logger.debug("Found regions to update due to reverse complement")
    # Apply reverse complement to REF column for RC rows
    if 'REF' in df.columns:
        # Only process rows that need reverse complement
        ref_to_rc = df.loc[rc_mask, 'REF']
        df.loc[rc_mask, 'REF'] = ref_to_rc.apply(reverse_complement)

    # Apply reverse complement to ALT column for RC rows
    if 'ALT' in df.columns:
        # Process each alternative allele
        alt_to_rc = df.loc[rc_mask, 'ALT']
        df.loc[rc_mask, 'ALT'] = alt_to_rc.apply(
            lambda alts: ','.join(reverse_complement(allele) for allele in alts.split(','))
        )

    # Apply reverse complement to Consensus column for RC rows
    if 'Consensus' in df.columns:
        consensus_to_rc = df.loc[rc_mask, 'Consensus']
        df.loc[rc_mask, 'Consensus'] = consensus_to_rc.apply(reverse_complement)

    # Swap complementary nucleotide counts
    nucleotide_pairs = [('A', 'T'), ('C', 'G')]
    for n1, n2 in nucleotide_pairs:
        if n1 in df.columns and n2 in df.columns:
            # Store temporary copies of the values
            temp_n1 = df.loc[rc_mask, n1].copy()
            temp_n2 = df.loc[rc_mask, n2].copy()

            # Swap the values
            df.loc[rc_mask, n1] = temp_n2
            df.loc[rc_mask, n2] = temp_n1

    return df


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
        update_map_data = patch_map_rc(map_data)
        vcf_header, vcf_data = read_vcf_file(args.vcf, logger, args.fasta)

        # Process VCF with mapping data
        processed_header, processed_data = process_vcf_with_map(
            vcf_header, vcf_data, update_map_data, logger
        )

        # Write output
        output_file = write_output_vcf(processed_header, processed_data, args.prefix, logger)

        logger.info(f"VCF patching completed successfully. Output written to {output_file}")
        return 0
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())

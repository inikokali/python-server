# find_amr_genes_module.py

import pandas as pd
import numpy as np
import os
from Bio.Blast.Applications import NcbiblastxCommandline
import io
import subprocess

def get_blast_version():
    try:
        # Run the 'blastx -version' command and capture the output
        result = subprocess.run(["blastx", "-version"], capture_output=True, text=True)
        
        # Print the captured output
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return f"Error: {result.stderr.strip()}"
    except FileNotFoundError:
        return "BLASTX not found. Please ensure BLAST+ is installed and accessible in your PATH."


def run_blastx_to_dataframe(query_file, db_path, num_threads=5):
    """
    Run BLASTX with the specified parameters and return the results as a DataFrame.

    Parameters:
    - query_file: Path to the query FASTA file.
    - db_path: Path to the BLAST database.
    - output_file: Path to the output file.
    - num_threads: Number of threads to use for BLASTX.

    Returns:
    - A pandas DataFrame with BLASTX results.
    """

    #### argis

    blastx_cline = NcbiblastxCommandline(
        query=query_file,
        db=db_path,
        outfmt="6 std qseq sseq slen",
        max_target_seqs=100,
        evalue=1e-10,  # Set the e-value threshold
        num_threads=num_threads
    )
    
    stdout, stderr = blastx_cline()
    
    # Check for errors in stderr
    if stderr:
        raise Exception(f"Error running BLASTX: {stderr}")

    #Use io.StringIO to read the stdout string into a DataFrame
    blast_df = pd.read_csv(io.StringIO(stdout), sep="\t", header=None)

    ####argis

    # Read the BLASTX results into a DataFrame directly from stdout
 #   column_names = [
 #       "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
 #       "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
 #       "qseq", "sseq", "slen"
 #   ]

    
    
    return blast_df


def save_dataframe_to_csv(dataframe, output_csv_path):
    """
    Save the given DataFrame to a CSV file.

    Parameters:
    - dataframe: The pandas DataFrame to save.
    - output_csv_path: Path to the output CSV file.
    """
    dataframe.to_csv(output_csv_path, index=False)


def process_blast_results(protein_annotation_file, assembly_blast_results, file_path, protein_start_filter=20, identity_threshold=70, coverage_threshold=70):
    """
    Process BLAST results, filter based on thresholds, and annotate.

    Parameters:
    - protein_annotation_file: Path to the protein annotation CSV file.
    - file_path: Path to the FASTA file for length calculations.
    - file_output: Path to the output CSV file.
    - protein_start_filter: Threshold for protein start position filter.
    - identity_threshold: Threshold for sequence identity.
    - coverage_threshold: Threshold for coverage percentage.

    Returns:
    - None
    """
    protein_annotation = pd.read_csv(protein_annotation_file, delimiter='\t')
    df = assembly_blast_results

    def read_fasta(file_path):
        sequences = {}
        with open(file_path, 'r') as file:
            entry_name = None
            sequence_length = 0
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if entry_name is not None:
                        sequences[entry_name] = sequence_length
                    entry_name = line[1:]
                    sequence_length = 0
                else:
                    sequence_length += len(line)
            if entry_name is not None:
                sequences[entry_name] = sequence_length
        return sequences

    sequences_dict = read_fasta(file_path)
    df['total_coverage'] = (df[3] / df[14]) * 100

    def get_length(name):
        for key in sequences_dict.keys():
            if name in key:
                return sequences_dict.get(key, 0)
        return 0

    df['partial_coverage'] = np.where(df[6] == 1, ((df[3] / (df[14] - df[3])) * 100), (df[3] / df[14]) * 100)
    df['contig_length'] = df[0].apply(get_length)

    filtered_df = df[(df.iloc[:, 8] <= protein_start_filter) | 
                     (df.iloc[:, 6] == 1) | 
                     (df.iloc[:, 7] == 1) | 
                     (df.iloc[:, 6] == df.loc[:, "contig_length"]) | 
                     (df.iloc[:, 7] == df.loc[:, "contig_length"])]

    filtered_df2 = filtered_df[((filtered_df.iloc[:, 2] >= identity_threshold) & 
                                (filtered_df.loc[:, 'total_coverage'] >= coverage_threshold)) | 
                               ((filtered_df.iloc[:, 6] == 1) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 7] == 1) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 6] == filtered_df.loc[:, "contig_length"]) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) | 
                               ((filtered_df.iloc[:, 7] == filtered_df.loc[:, "contig_length"]) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold))]
    filtered_df2 = filtered_df2.reset_index(drop=True)

    df1 = pd.DataFrame()
    for index, row in filtered_df2.iterrows():
        inner_break = False
        node = row.iloc[0]
        min_value = min(row.iloc[7], row.iloc[6])
        max_value = max(row.iloc[7], row.iloc[6])
        for index_2, row_2 in filtered_df2[filtered_df2[0] == node].iterrows():
            if row.equals(row_2):
                continue
            else:
                min_value_2 = min(row_2.iloc[7], row_2.iloc[6])
                max_value_2 = max(row_2.iloc[7], row_2.iloc[6])
                if (min_value >= max_value_2) or (max_value <= min_value_2):
                    continue
                elif (min_value >= min_value_2) & (min_value <= max_value_2) & (max_value >= max_value_2):
                    ratio = float((max_value_2 - min_value) / (max_value_2 - min_value_2))
                    if ratio >= 0.2:
                        if "_ARO:" in row.iloc[1] and "_ARO:" in row_2.iloc[1]:
                            if row.iloc[10] < row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] <= row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] < row_2.iloc[10] and row.iloc[11] >= row_2.iloc[11]:
                                continue
                            elif (row.iloc[10] == row_2.iloc[10] and 
                                  row.iloc[11] == row_2.iloc[11] and 
                                  row.loc["partial_coverage"] == row_2.loc["partial_coverage"]):
                                continue
                            elif (row.iloc[10] == row_2.iloc[10] and 
                                  row.iloc[11] == row_2.iloc[11] and 
                                  row.loc["partial_coverage"] > row_2.loc["partial_coverage"]):
                                inner_break = True
                                break
                            else:
                                inner_break = True
                                break
                        elif "_ARO:" not in row.iloc[1] and "_ARO:" in row_2.iloc[1]: 
                            continue
                        elif "_ARO:" in row.iloc[1] and "_ARO:" not in row_2.iloc[1]: 
                            inner_break = True
                            break
                        else:
                            if row.iloc[10] < row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] <= row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] < row_2.iloc[10] and row.iloc[11] >= row_2.iloc[11]:
                                continue
                            elif (row.iloc[10] == row_2.iloc[10] and 
                                  row.iloc[11] == row_2.iloc[11] and 
                                  row.loc["partial_coverage"] == row_2.loc["partial_coverage"]):
                                continue
                            elif (row.iloc[10] == row_2.iloc[10] and 
                                  row.iloc[11] == row_2.iloc[11] and 
                                  row.loc["partial_coverage"] > row_2.loc["partial_coverage"]):
                                inner_break = True
                                break
                            else:
                                inner_break = True
                                break
                    else:
                        continue
                elif (min_value <= min_value_2) & (max_value <= max_value_2) & (max_value >= min_value_2):
                    ratio = float((max_value - min_value_2) / (max_value_2 - min_value_2))
                    if ratio >= 0.2:
                        if "_ARO:" in row.iloc[1] and "_ARO:" in row_2.iloc[1]:
                            if row.iloc[10] < row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] <= row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] < row_2.iloc[10] and row.iloc[11] >= row_2.iloc[11]:
                                continue
                            elif (row.iloc[10] == row_2.iloc[10] and 
                                  row.iloc[11] == row_2.iloc[11] and 
                                  row.loc["partial_coverage"] == row_2.loc["partial_coverage"]):
                                continue
                            elif (row.iloc[10] == row_2.iloc[10] and 
                                  row.iloc[11] == row_2.iloc[11] and 
                                  row.loc["partial_coverage"] > row_2.loc["partial_coverage"]):
                                inner_break = True
                                break
                            else:
                                inner_break = True
                                break
                        elif "_ARO:" not in row.iloc[1] and "_ARO:" in row_2.iloc[1]: 
                            continue
                        elif "_ARO:" not in row.iloc[1] and "_ARO:" in row_2.iloc[1]: 
                            inner_break = True
                            break
                        else:
                            if row.iloc[10] < row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] <= row_2.iloc[10] and row.iloc[11] > row_2.iloc[11]:
                                continue
                            elif row.iloc[10] < row_2.iloc[10] and row.iloc[11] >= row_2.iloc[11]:
                                continue
                            elif (row.iloc[10] <= row_2.iloc[10]) and (row.iloc[11] >= row_2.iloc[11]):
                                continue
                            else:
                                inner_break = True
                                break

        if not inner_break:
            row_df = pd.DataFrame([row])
            df1 = pd.concat([df1, row_df], ignore_index=True)

    df1 = df1.drop_duplicates(subset=[11, 12])
    df1['Comments'] = 'Complete gene based on threshold'
    df1.loc[((df1[6] == 1) | (df1[7] == 1)), 'Comments'] = 'Partial gene in the start of Contig'
    df1.loc[((df1[6] == df1['contig_length']) | (df1[7] == df1['contig_length'])), 'Comments'] = 'Partial gene in the end of Contig'
    df1['GeneImportance'] = "Major"
    df1.loc[df1[1].str.contains('_ARO:'), 'GeneImportance'] = 'Minor'

    new_names = {
        0: 'Contig ID',
        1: 'ID',
        2: "Identity",
        3: "Length Alignment",
        4: "Mismatches",
        5: "Gaps",
        6: "Start",
        7: "Stop",
        8: "Start In Reference Protein",
        9: "Stop In Reference Protein",
        10: "E-value",
        11: "BitScore",
        12: "Genome Protein Sequence",
        13: "Protein Reference Sequence",
        14: "Alignment Length", 
        "total_coverage": "% Coverage of Reference Sequence",
        "contig_length": "Length of the Identified Contig"
    }
    df1 = df1.rename(columns=new_names)

    protein_annotation['ID'] = protein_annotation['ID'].astype(str)

    for index, row in df1.iterrows():
        value = row['ID']
        row_extra = protein_annotation[protein_annotation['ID'].apply(lambda x: x in value)]
        
        if not row_extra.empty:
            for column in row_extra.columns:
                new_column_name = f"{column}"
                df1.loc[index, new_column_name] = row_extra[column].values[0]

    df3 = df1.sort_values(by='GeneImportance')
    #df3.to_csv(file_output, index=False)
    #print(f"Final results have been saved to {file_output}.")
    return df3

# find_amr_genes_module.py

import pandas as pd
import numpy as np
import pandas as pd
import os
from Bio.Blast.Applications import NcbiblastxCommandline


def run_blastx_to_dataframe(query_file, db_path, output_file, num_threads=5):
    # Define the blastx command with the provided arguments
    blastx_cline = NcbiblastxCommandline(
        query=query_file,
        db=db_path,
        outfmt="6 std qseq sseq slen",
        max_target_seqs=1000,
        out=output_file,
        num_threads=num_threads
    )
    
    # Execute the command and capture the output
    stdout, stderr = blastx_cline()
    
    # Check if the output file is created
    if not os.path.exists(output_file):
        raise FileNotFoundError(f"BLASTX output file '{output_file}' not found.")
    
    # Read the BLASTX results into a DataFrame
    column_names = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
        "qseq", "sseq", "slen"
    ]
    blast_df = pd.read_csv(output_file, sep="\t", header=None, names=column_names)
    
    return blast_df


def process_blast_results(protein_annotation_file, assembly_blast_results, file_path, file_output, protein_start_filter=20, identity_threshold=70, coverage_threshold=70):
    # Read the CSV file into a DataFrame
    protein_annotation = pd.read_csv(protein_annotation_file, delimiter='\t')

    # Read the results of BLAST for assembly
    df = pd.read_csv(assembly_blast_results, sep="\t", header=None)

    # Function for reading the FASTA file
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

    # Read the FASTA file 
    sequences_dict = read_fasta(file_path)

    # Find the total coverage
    df['total_coverage'] = (df[3] / df[14]) * 100

    # Function to get length from dictionary by matching name to a substring of the keys
    def get_length(name):
        for key in sequences_dict.keys():
            if name in key:  # Check if the name is a substring of the key
                return sequences_dict.get(key, 0)  # Return the corresponding value if found
        return 0  # Default to 0 if no match is found

    # Create the new column using numpy.where()
    df['partial_coverage'] = np.where(df[6] == 1, ((df[3] / (df[14]-df[3])) * 100), (df[3] / df[14]) * 100)

    # Assuming df is your DataFrame
    df['contig_length'] = df[0].apply(get_length)

    # Filter the dataframe df
    filtered_df = df[(df.iloc[:, 8] <= protein_start_filter) | 
                     (df.iloc[:, 6] == 1) | 
                     (df.iloc[:, 7] == 1) | 
                     (df.iloc[:, 6] == df.loc[:, "contig_length"]) | 
                     (df.iloc[:, 7] == df.loc[:, "contig_length"])]

    filtered_df2 = filtered_df[((filtered_df.iloc[:, 2] >= identity_threshold) & 
                                (filtered_df.loc[:,'total_coverage'] >= coverage_threshold)) | 
                               ((filtered_df.iloc[:, 6] == 1) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 7] == 1) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 6] == filtered_df.loc[:, "contig_length"]) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) | 
                               ((filtered_df.iloc[:, 7] == filtered_df.loc[:, "contig_length"]) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold))]
    filtered_df2 = filtered_df2.reset_index(drop=True)

    # Find the best matches
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
                            elif (row.iloc[10


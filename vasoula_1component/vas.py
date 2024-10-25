import find_amr_genes_module

# File paths and parameters
query_file = "S3323_contigs.fasta"  
db_path = "databases/all/all_amr"
num_threads = 4
protein_annotation_file = "databases/genes_annotation_databases.csv"
protein_start_filter = 20  # Threshold for protein start position filter
identity_threshold = 70  # Threshold for sequence identity
coverage_threshold = 70  # Threshold for coverage percentage
output_csv_file = "final_results.csv"  # Path to save the final CSV

# Print the input parameters
print(f"Query file: {query_file}")
print(f"Database path: {db_path}")
print(f"Number of threads: {num_threads}")
print(f"Protein annotation file: {protein_annotation_file}")
print(f"Protein start filter: {protein_start_filter}")
print(f"Identity threshold: {identity_threshold}")
print(f"Coverage threshold: {coverage_threshold}")

# Run BLASTX and store the results in a DataFrame
blast_results_df = find_amr_genes_module.run_blastx_to_dataframe(
    query_file, db_path, num_threads
)

# Process the BLAST results
final_results = find_amr_genes_module.process_blast_results(
    protein_annotation_file, blast_results_df, query_file, 
    protein_start_filter, identity_threshold, coverage_threshold
)

# Save the final results to a CSV file
final_results.to_csv(output_csv_file, index=False)

# print(f"Final results have been saved to {output_csv_file}")

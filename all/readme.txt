take the files from each database 

at protein_fasta_protein_homolog_model.fasta | tr '|' ' ' | sed 's/gb //g' | sed 's/ ARO/_ARO/g' > card_concatenated_prot.fasta


cat resfinder_protein.fasta | awk '/^>/ {if (!match(seq, /\*/)) {print header; print seq} header=$0; seq=""; next} {seq = seq $0} END {if (!match(seq, /\*/)) {print header; print seq}}' 

cat concatenated_amrfinderPlus.fasta card_concatenated_prot.fasta resfinder_protein_filtered.fasta | sed '/^$/d' > all_proteins_amr.fasta 

cd-hit -i all_proteins_amr.fasta -o nr100_all_proteins_amr.fasta -c 1.00 -aL 1.0 -n 5 -M 2000

makeblastdb -in nr100_all_proteins_amr.fasta -parse_seqids -out nr_all_amr -dbtype prot

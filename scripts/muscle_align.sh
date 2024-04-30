for fasta_file in $(ls ../data/manual_checking/*/*/seq4dnds/prot_orths_*.fasta)
do
    # Replace 'prot_orths_orthomcl.fasta' with 'prot_orths_orthomcl.afa' in the output file path
    output_file="${fasta_file%.fasta}.afa"
    
    # Run the muscle command
    muscle -align "$fasta_file" -output "$output_file"
done

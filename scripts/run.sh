#!/bin/bash
set -e
set -x 
mkdir -p "../data"
# Adjust the following lines if you plan to run this script for a new organism
wget -P ../ https://zenodo.org/records/11094116/files/data.tar.gz
tar -xzf ../data.tar.gz -C ..
gff_path="../data/TriTrypDB-65_TbruceiTREU927.gff"
genome_fasta_path="../data/TriTrypDB-65_TbruceiTREU927_Genome.fasta"
proteome_fasta_path="../data/TriTrypDB-65_TbruceiTREU927_AnnotatedProteins.fasta"
genus_proteomes="../data/Trypanosomes_prot.fasta"
genus_cds="../data/Trypanosomes_cds.fasta"
# Note: In the orthologous tables, the first column should correspond to the query organism, and the second column should be the closest organism
orths_tbl_path="../data/rbh_tb927_trypanosomes.tsv"



# Adjust the path to the following tools if needed
diamond_path="diamond"
interproscan_path="interproscan.sh"
muscle_path="muscle"
pal2nal_path="pal2nal.pl"
PhyML_path="PhyML"
hyphy_path="hyphy"

# Adjust the number of available threads:
threads=20


#First, we extract the coding sequence and flanking regions to use as query in the next steps.
python extract_coding_and_flanking.py $gff_path $genome_fasta_path

$interproscan_path -exclappl Coils,AntiFam,MobiDBLite,PANTHER -i $proteome_fasta_path -f tsv -dp -o ../data/ipr_proteome.tsv -cpu $threads

mkdir -p ../data/diamond_db

proteomedb="../data/diamond_db/proteomedb"
queryfasta="../data/cds_and_flanking.fasta"
output_pref="../data/genome_vs_proteome"

$diamond_path makedb --in $proteome_fasta_path -d $proteomedb

# Diamond reports the alignments in different formats. fmt6 will show the alignments in a tabular format. However, if there is a frameshift in the 
# alignments, it can be shown if we use the pairwise format. Therefore, we would get the output of diamond in both pairwise and fmt6 formats and combine them

$diamond_path blastx -d $proteomedb -q $queryfasta -o ${output_pref}_fmt6.txt \
--ultra-sensitive --min-orf 1 -F 15 --masking 0 --evalue 1e-10 \
--outfmt 6 qseqid sseqid qframe pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen 

$diamond_path blastx -d $proteomedb -q $queryfasta -o ${output_pref}_pairwise.txt \
--ultra-sensitive --min-orf 1 -F 15 --masking 0 --outfmt 0 --evalue 1e-10

python convert_dmnd_pairwise2table.py 

python mark4validation_by_intragenome_alignment.py

$interproscan_path -exclappl Coils,AntiFam,MobiDBLite,PANTHER -i ../data/subject_adjusted_proteinseq.fasta -f tsv -dp -o ../data/ipr_sadj_seq.tsv -cpu $threads

python find_genes_with_domcov_inc_aft_correction.py

python collect_orthologs4msa.py $genus_proteomes $genus_cds $orths_tbl_path

for fasta_file in $(ls ../data/manual_checking/*/*/seq4dnds/prot_orths*.fasta)
do
    output_file="${fasta_file%.fasta}.afa"
    $muscle_path -align "$fasta_file" -output "$output_file"
done

python prepare_muscle_out4pal2nal.py

for fasta_file in $(ls ../data/manual_checking/*/*/seq4dnds/prot_orths_4p2n.afa)
do
    tmp_file1="${fasta_file%_4p2n.afa}.nal"
    output_file="${tmp_file1/prot_/}"

    tmp_file2="${fasta_file/prot/cds}"
    cds_file="${tmp_file2/_4p2n.afa/.fasta}"

    $pal2nal_path $fasta_file $cds_file -nogap -output fasta > $output_file
done

python convert_nal2phylip.py

files=$(ls ../data/manual_checking/*/*/seq4dnds/*.PHYLIP)
for file in $files
do
    $PhyML_path -i $file -d nt --sequential --model GTR
done

python backconvert_axidtree.py

ls ../data/manual_checking/*/*/seq4dnds/orths.nal | \
parallel "nwk_file={.}_geneid.nwk; hyphy relax --alignment {} --tree \$nwk_file --test test"

#move relax output to parent directory:
find ../data/manual_checking -type f -name "orths.nal.RELAX.json" -exec sh -c 'mv "$0" "$(dirname "$0")/../"' {} \;

python generate_summary.py $orths_tbl_path

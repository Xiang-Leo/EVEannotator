# This model can download NR database for Diamond blast.

mkdir data database
cd data

wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5
# Taxnomy
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
unzip taxdmp.zip

# Check file integrity 
md5sum -c nr.gz.md5

# construct NR database
diamond makedb --in nr.gz --db ../database/nr --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp
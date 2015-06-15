!#/bin/bash
# This script downloads and formats data from refseq to create training sets for viral classification
cd /global/projectb/scratch/arrivers/rna_virus/results/20150223
mkdir trainingdata
cd trainingdata
#setup
mkdir archaea
mkdir bacteria
mkdir plasmid
mkdir mitochondrion
mkdir plastid
mkdir protozoa
mkdir plant
mkdir invertebrate
mkdir vertebrate_mammalian
mkdir vertebrate_other
mkdir fungi

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER 
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/archaea.*.genomic.fna.gz 
mv archaea.* archaea/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.*.genomic.fna.gz
mv bacteria.* bacteria/ 
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.*.genomic.fna.gz 
mv plasmid.* plasmid/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.*.genomic.fna.gz
mv mitochondrion.* mitochondrion/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.*.genomic.fna.gz
mv plastid.* plastid/ 
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/protozoa.*.genomic.fna.gz
mv protozoa.* protozoa/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/plant.*.genomic.fna.gz
mv plant.* plant/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/invertebrate/invertebrate.*.genomic.fna.gz
mv invertebrate.* invertebrate/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.*.genomic.fna.gz 
mv vertebrate_mammalian* vertebrate_mammalian/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/vertebrate_other.*.genomic.fna.gz 
mv vertebrate_other* vertebrate_other/
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.*.genomic.fna.gz
mv fungal.* fungal/


wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.genomic.fna.gz
gunzip viral.*.genomic.fna.gz
module load jgibio

jgi_tax_filter.pl --rank "no rank" --name "ssRNA viruses" --taxdump /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/ --gi_taxid /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/gi_taxid_nucl.dmp --out ssRNAviruses viral.*.genomic.fna
jgi_tax_filter.pl --rank "no rank" --name "dsRNA viruses" --taxdump /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/ --gi_taxid /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/gi_taxid_nucl.dmp --out dsRNAviruses viral.*.genomic.fna
jgi_tax_filter.pl --rank "no rank" --name "dsDNA viruses, no RNA stage" --taxdump /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/ --gi_taxid /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/gi_taxid_nucl.dmp --out dsDNAviruses viral.*.genomic.fna
jgi_tax_filter.pl --rank "no rank" --name "ssDNA viruses" --taxdump /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/ --gi_taxid /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/gi_taxid_nucl.dmp --out ssDNAviruses viral.*.genomic.fna
jgi_tax_filter.pl --rank "no rank" --name "Retro-transcribing viruses" --taxdump /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/ --gi_taxid /global/dna/projectdirs/MEP/db/NCBI/taxonomy/DEFAULT/gi_taxid_nucl.dmp --out retroviruses viral.*.genomic.fna

cat ssRNAviruses/* dsRNAviruses/* ssDNAviruses/* dsDNAviruses/* retroviruses/* > mainviruses.fna

perl -e ' $count=0; $len=0; while(<>) { s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) { print "\n" } s/ |$/\t/; $count++; $_ .= "\t"; } else { s/ //g; $len += length($_) } print $_; } print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; ' mainviruses.fna > mainviruses.tab
perl -e ' @cols=(0); while(<>) { s/\r?\n//; @F=split /\t/, $_; print join("\t", @F[@cols]), "\n" } warn "\nChose columns ", join(", ", @cols), " for $. lines\n\n" ' mainviruses.tab > mainvirusids.txt 

perl -e ' $count=0; $len=0; while(<>) { s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) { print "\n" } s/ |$/\t/; $count++; $_ .= "\t"; } else { s/ //g; $len += length($_) } print $_; } print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; ' viral.*.genomic.fna > allviruses.tab
perl -e ' @cols=(0); while(<>) { s/\r?\n//; @F=split /\t/, $_; print join("\t", @F[@cols]), "\n" } warn "\nChose columns ", join(", ", @cols), " for $. lines\n\n" ' allviruses.tab > allvirusids.txt 

python invert.py mainvirusids.txt allvirusids.txt >othervirusids.txt

mkdir viruses_other
perl -e ' ($id,$fasta)=@ARGV; open(ID,$id); while (<ID>) { s/\r?\n//; /^>?(\S+)/; $ids{$1}++; } $num_ids = keys %ids; open(F, $fasta); $s_read = $s_wrote = $print_it = 0; while (<F>) { if (/^>(\S+)/) { $s_read++; if($ids{$1}) { $s_wrote++; $print_it = 1; delete $ids{$1} } else { $print_it = 0 } }; if ($print_it) { print $_ } }; END { warn "Searched $s_read FASTA records.\nFound $s_wrote IDs out of $num_ids in the ID list.\n" } ' othervirusids.txt viral.*.genomic.fna > viruses_other/viruses.genomic.fna 
rm mainvirus*
rm othervirus*
rm allvirus*
cd ../..


#setup sqlite db (this takes about 30 min)
module load sqlite3
sqlite3 taxa.db < taxdmpsetup.sql
#create directories for split genomes

mkdir traininggenomes
cd traininggenomes
mkdir archaea
mkdir bacteria
mkdir plasmid
mkdir mitochondrion
mkdir plastid
mkdir protozoa
mkdir plant
mkdir invertebrate
mkdir vertebrate_mammalian
mkdir vertebrate_other
mkdir fungi
mkdir ssDNAviruses
mkdir dsDNAviruses
mkdir ssRNAviruses
mkdir dsRNAviruses
mkdir retroviruses
mkdir viruses_other
cd ..

# create submission scripts for genome classification:
echo "module load biopython; zcat trainingdata/archaea/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/archaea -i -">gqsub1.sh
echo "module load biopython; zcat trainingdata/bacteria/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/bacteria -i -">gqsub2.sh
echo "module load biopython; zcat trainingdata/plasmid/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/plasmid -i -">gqsub3.sh
echo "module load biopython; zcat trainingdata/mitochondrion/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/mitochondrion -i -">gqsub4.sh
echo "module load biopython; zcat trainingdata/plastid/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/plastid -i -">gqsub5.sh
echo "module load biopython; zcat trainingdata/protozoa/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/protozoa -i -">gqsub6.sh
echo "module load biopython; zcat trainingdata/plant/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/plant -i -">gqsub7.sh
echo "module load biopython; zcat trainingdata/invertebrate/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/invertebrate -i -">gqsub8.sh
echo "module load biopython; zcat trainingdata/vertebrate_mammalian/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/vertebrate_mammalian -i -">gqsub9.sh
echo "module load biopython; zcat trainingdata/vertebrate_other/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/vertebrae_other -i -">gqsub10.sh
echo "module load biopython; zcat trainingdata/fungi/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/fungi -i -">gqsub11.sh
echo "module load biopython; cat trainingdata/ssDNAviruses/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/ssDNAviruses -i -">gqsub12.sh
echo "module load biopython; cat trainingdata/dsDNAviruses/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/dsDNAviruses -i -">gqsub13.sh
echo "module load biopython; cat trainingdata/ssRNAviruses/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/ssRNAviruses -i -">gqsub14.sh
echo "module load biopython; cat trainingdata/dsRNAviruses/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/dsRNAviruses -i -">gqsub15.sh
echo "module load biopython; cat trainingdata/retroviruses/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/retroviruses -i -">gqsub16.sh
echo "module load biopython; cat trainingdata/viruses_other/* | ./genome_parserdb.py -t taxa.db -o traininggenomes/viruses_other -i -">gqsub17.sh

# submit scripts
 bash gqsuball.sh

#Cleanup
#mkdir gqsub
#mv gqsub* gqsub

# delete records that are very small
find traininggenomes/archaea/ -name "*.fasta" -size -50k -delete
find traininggenomes/bacteria/ -name "*.fasta" -size -50k -delete
find traininggenomes/plasmid/ -name "*.fasta" -size -5k -delete
find traininggenomes/fungi/ -name "*.fasta" -size -120k -delete
find traininggenomes/mitochondrion/ -name "*.fasta" -size -10k -delete
find traininggenomes/plastid/ -name "*.fasta" -size -50k -delete
find traininggenomes/protozoa/ -name "*.fasta" -size -500k -delete
find traininggenomes/plant/ -name "*.fasta" -size -1000k -delete
find traininggenomes/vertebrate_mammalian/ -name "*.fasta" -size -50k -delete
find traininggenomes/vertebrate_other/ -name "*.fasta" -size -50k -delete
find traininggenomes/invertebrate/ -name "*.fasta" -size -50k -delete

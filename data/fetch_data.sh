echo "Fetching the 1000GP vcf (GRCh37) file"
wget --continue \
ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz \
estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz

echo "Uncompress 1000GP vcf (GRCh37) file"
gunzip estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz


mkdir -p 'reads'

#TODO(gamazeps): make sure to avoid the cram files as well as the exome alignments data.
echo "Fetching the reads (bam) from 1000GP ftp"
set e
wget --recursive --continue \
    --accept='*.mapped.*.bam*' \
    --accept-regex='.*/alignment/' \
    'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/' \
    --directory-prefix='reads

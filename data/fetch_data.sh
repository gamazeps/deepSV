echo "Fetching the 1000GP vcf (GRCh37) file"
wget --continue \
ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz \
estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz

echo "Uncompress 1000GP vcf (GRCh37) file"
gunzip estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz

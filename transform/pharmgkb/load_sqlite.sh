echo 'Loading TSVs into sqlite database...'
python3 transform/pharmgkb/load_sqlite.py

echo 'Loading haplotypes from pharmGKB API...'
# since dvc deletes -o files before running cmd, maintain a backup copy
if [ -f "outputs/pharmgkb/pharmgkb_haplotype_requests_cache.sqlite.bak" ]; then
  cp outputs/pharmgkb/pharmgkb_haplotype_requests_cache.sqlite.bak outputs/pharmgkb/pharmgkb_haplotype_requests_cache.sqlite
fi
python3 transform/pharmgkb/fetch_haplotypes.py
# update it with new file
if [ -f "outputs/pharmgkb/pharmgkb_haplotype_requests_cache.sqlite" ]; then
  cp outputs/pharmgkb/pharmgkb_haplotype_requests_cache.sqlite outputs/pharmgkb/pharmgkb_haplotype_requests_cache.sqlite.bak
fi


echo 'Loading variants from pharmGKB API...'
# since dvc deletes -o files before running cmd, maintain a backup copy
if [ -f "outputs/pharmgkb/pharmgkb_variant_requests_cache.sqlite.bak" ]; then
  cp outputs/pharmgkb/pharmgkb_variant_requests_cache.sqlite.bak outputs/pharmgkb/pharmgkb_variant_requests_cache.sqlite
fi
python3 transform/pharmgkb/fetch_variants.py
# update it with new file
if [ -f "outputs/pharmgkb/pharmgkb_variant_requests_cache.sqlite" ]; then
  cp outputs/pharmgkb/pharmgkb_variant_requests_cache.sqlite outputs/pharmgkb/pharmgkb_variant_requests_cache.sqlite.bak
fi


echo 'Loading refsnp information, this may take a while...'
# since dvc deletes -o files before running cmd, maintain a backup copy
if [ -f "outputs/pharmgkb/refsnp_requests_cache.sqlite.bak" ]; then
  cp outputs/pharmgkb/refsnp_requests_cache.sqlite.bak outputs/pharmgkb/refsnp_requests_cache.sqlite
fi
python3 transform/pharmgkb/load_alleles.py
# update it with new file
if [ ! -f "outputs/pharmgkb/refsnp_requests_cache.sqlite.bak" ]; then
  cp outputs/pharmgkb/refsnp_requests_cache.sqlite outputs/pharmgkb/refsnp_requests_cache.sqlite.bak
fi

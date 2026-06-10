Download ENA metadata tables:

read_run = sequencing runs and FASTQ links
read_experiment = experiment/library metadata
sample = BioSample metadata
assembly = assembly metadata and assembly FASTA links
optionally read_study = project/study metadata
Link tables using accession columns:

run_accession
experiment_accession
sample_accession
study_accession
assembly_accession
Download FASTQs and assemblies using the FTP/HTTPS links in the metadata tables.

Below is a complete practical workflow for Clostridioides difficile, NCBI taxid 1496.

1. Set variables
Use tax_eq(1496) for exact taxon only, or tax_tree(1496) if you want children/subtaxa too.

TAX_QUERY="tax_eq(1496)"
DATE="2026June10"
API="https://www.ebi.ac.uk/ena/portal/api/search"

If you want broader matching:

TAX_QUERY="tax_tree(1496)"

2. Download read runs
This gives run accessions, experiment accessions, sample accessions, study accessions, and FASTQ links.

curl -L -G "$API" \
  --data-urlencode "result=read_run" \
  --data-urlencode "query=${TAX_QUERY}" \
  --data-urlencode "fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,study_accession,secondary_study_accession,scientific_name,tax_id,instrument_platform,instrument_model,library_strategy,library_layout,read_count,base_count,fastq_ftp,fastq_md5,submitted_ftp,sra_ftp,first_public,last_updated" \
  --data-urlencode "format=tsv" \
  --data-urlencode "limit=0" \
  -o ena_cdiff_read_run_${DATE}.tsv

3. Download read experiments
Use read_experiment, not experiment.

curl -L -G "$API" \
  --data-urlencode "result=read_experiment" \
  --data-urlencode "query=${TAX_QUERY}" \
  --data-urlencode "fields=experiment_accession,study_accession,sample_accession,secondary_sample_accession,experiment_title,center_name,instrument_platform,instrument_model,library_strategy,library_source,library_selection,library_layout,first_public,last_updated" \
  --data-urlencode "format=tsv" \
  --data-urlencode "limit=0" \
  -o ena_cdiff_read_experiment_${DATE}.tsv

4. Download samples
curl -L -G "$API" \
  --data-urlencode "result=sample" \
  --data-urlencode "query=${TAX_QUERY}" \
  --data-urlencode "fields=sample_accession,secondary_sample_accession,tax_id,scientific_name,sample_alias,sample_title,center_name,first_public,last_updated,country,collection_date,host,isolation_source" \
  --data-urlencode "format=tsv" \
  --data-urlencode "limit=0" \
  -o ena_cdiff_sample_${DATE}.tsv

If one of those sample fields is not accepted by ENA, check the available fields with:

curl -L "https://www.ebi.ac.uk/ena/portal/api/returnFields?result=sample" \
  -o ena_sample_fields.tsv

less ena_sample_fields.tsv

Then remove unsupported fields from the command.

5. Download assemblies
curl -L -G "$API" \
  --data-urlencode "result=assembly" \
  --data-urlencode "query=${TAX_QUERY}" \
  --data-urlencode "format=tsv" \
  --data-urlencode "limit=0" \
  -o ena_cdiff_assembly_${DATE}.tsv

I recommend first downloading assemblies without specifying fields, because ENA assembly fields vary more often.

Inspect the header:

head -1 ena_cdiff_assembly_${DATE}.tsv | tr '\t' '\n'

Look for fields such as:

assembly_accession
sample_accession
study_accession
scientific_name
tax_id
assembly_name
assembly_level
fasta_ftp
submitted_ftp

If you want to see all available assembly fields:

curl -L "https://www.ebi.ac.uk/ena/portal/api/returnFields?result=assembly" \
  -o ena_assembly_fields.tsv

less ena_assembly_fields.tsv

6. Optional: download studies
curl -L -G "$API" \
  --data-urlencode "result=read_study" \
  --data-urlencode "query=${TAX_QUERY}" \
  --data-urlencode "format=tsv" \
  --data-urlencode "limit=0" \
  -o ena_cdiff_read_study_${DATE}.tsv

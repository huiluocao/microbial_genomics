This includes how to acquire the metadata from enterobase and genome sequnces.
1. Metadata could be downloaded directly from database with "save to local file".
2. Genome/assemble could be downloaded with the script:
before running, need to set: database=clostridium and export ENTEROBASE_API_TOKEN="your_new_enterobase_api_token"(token could be requested from admin of enterobase via email)
e.g. 
python api_clo_dl.py --barcode CLO_AA0007AA_AS --outdir assemblies

Or with a barcode file:
python api_clo_dl.py --barcode-file barcodes.txt --outdir assemblies
Your barcodes.txt should look like:
CLO_AA0007AA_AS
CLO_AA0487AA_AS
CLO_AA1234AA_AS

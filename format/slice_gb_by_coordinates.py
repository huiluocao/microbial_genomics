from Bio import SeqIO

with open('lsaD_contig_coordinates.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:  # Skip empty lines
            continue
            
        i = line.split('\t')
        if len(i) < 4:  # Ensure there are enough columns
            print(f"Skipping malformed line: {line}")
            continue
            
        contig_id = i[0]
        target_id = i[1]
        start = int(i[2])
        end = int(i[3])
        
        print(f"Processing coordinates: {end}")  # More descriptive print
        
        try:
            with open(f'gb/{contig_id}.gb', 'r') as g:
                for rec in SeqIO.parse(g, 'genbank'):
                    if rec.id == target_id:
                        # Create a slice of the record
                        sliced_rec = rec[start:end]
                        sliced_rec.id = f'{contig_id}:{rec.id}:{start}:{end}'
                        sliced_rec.description = ''
                        
                        print(f"Writing sliced record: {sliced_rec.id}")
                        
                        # Use append mode but write each record properly
                        with open('lsaD_contig10k.gb', 'a') as h:
                            SeqIO.write(sliced_rec, h, 'genbank')
        except FileNotFoundError:
            print(f"File not found: gb/{contig_id}.gb")
        except Exception as e:
            print(f"Error processing {contig_id}: {str(e)}")

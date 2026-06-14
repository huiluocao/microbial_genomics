/Applications/anaconda3/bin/poppunk --create-db --r-files 62fna.list --output bf2_poppunk --threads 4 --plot-fit 4 --ignore-length --full-db --overwrite --min-k 15

Huiluos-MacBook-Pro:Bf2 hcao$ /Applications/anaconda3/bin/poppunk --fit-model --distances bf2_poppunk/bf2_poppunk.dists --output bf2_strain_db_62 --full-db --ref-db bf2_poppunk --K 3 --dbscan
PopPUNK (POPulation Partitioning Using Nucleotide Kmers)
Mode: Fitting model to reference database

Failed to find distinct clusters in this dataset
Huiluos-MacBook-Pro:Bf2 hcao$ /Applications/anaconda3/bin/poppunk --fit-model --distances bf2_poppunk/bf2_poppunk.dists --output bf2_strain_db_62 --full-db --ref-db bf2_poppunk --K 3 --cytoscape
PopPUNK (POPulation Partitioning Using Nucleotide Kmers)
Mode: Fitting model to reference database

/Applications/anaconda3/lib/python3.6/site-packages/sklearn/mixture/base.py:268: ConvergenceWarning: Initialization 5 did not converge. Try different init parameters, or increase max_iter, tol or check for degenerate data.
  % (init + 1), ConvergenceWarning)
Fit summary:
	Avg. entropy of assignment	0.0002
	Number of components used	3

Scaled component means:
	[ 0.77373161  0.76647075]
	[ 0.05130682  0.25133931]
	[ 0.5704889   0.65424554]

Network summary:
	Components	32
	Density	0.0471
	Transitivity	1.0000
	Score	0.9529
Writing cytoscape output

Done
Huiluos-MacBook-Pro:Bf2 hcao$ /Applications/anaconda3/bin/poppunk --fit-model --distances bf2_poppunk/bf2_poppunk.dists --output bf2_strain_db_62_k2 --full-db --ref-db bf2_poppunk --K 2 --cytoscape
PopPUNK (POPulation Partitioning Using Nucleotide Kmers)
Mode: Fitting model to reference database

Fit summary:
	Avg. entropy of assignment	0.0032
	Number of components used	2

Scaled component means:
	[ 0.75153973  0.76059142]
	[ 0.31160871  0.40657741]

Network summary:
	Components	12
	Density	0.0925
	Transitivity	0.9426
	Score	0.8554
Writing cytoscape output



for file in fna/*.fna; do mash sketch -S 42 -k 21 -s 1000 -o input $file; mash screen /Users/hcao/Desktop/plasmids/plasmids__2019_03_05__v0/2019_03_05.msh $file -v 0.1 -i 0.99 > $file'_plasmid_mash.tsv'; done &

for file in K2E3.fna Q1F2.fna Q6B8.fna; do mash sketch -S 42 -k 21 -s 1000 -o input $file; mash screen /Users/hcao/Desktop/plasmids/plasmids__2019_03_05__v0/2019_03_05.msh $file -v 0.1 -i 0.99 > $file'_plasmid_mash.tsv'; done &

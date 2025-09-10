IMGDIR = img

mlab = matlab -nosplash -nodesktop -batch


$(IMGDIR)/dimensionality.pdf:	mk_dimensionality.m result_corn_dataset.mat result_pulp_dataset.mat result_gene_dataset.mat
	$(mlab) "mk_dimensionality;"

$(IMGDIR)/mox_corn.pdf:	mk_mox_corn.m moxplot.m dataset_corn.mat
	$(mlab) "mk_mox_corn;"

$(IMGDIR)/mox_pulp.pdf:	mk_mox_pulp.m moxplot.m ../../modes/analysis/modes_kvarnsveden64.mat
	$(mlab) "mk_mox_pulp;"

$(IMGDIR)/mox_gene.pdf:	mk_mox_gene.m data/genes.mat
	$(mlab) "mk_mox_gene;"

result_corn_dataset.mat:	sim_corn_dataset.m
	$(mlab) "sim_corn_dataset;"

result_pulp_dataset.mat:	sim_pulp_dataset.m
	$(mlab) "sim_pulp_dataset;"

result_gene_dataset.mat:	sim_gene_dataset.m
	$(mlab) "sim_gene_dataset;"

dataset_corn.mat:	mk_datasets.m
	$(mlab) "mk_datasets;"


# --- old ---

$(IMGDIR)/benchmark_simulated.eps:	mox_comparisson.m
	$(mlab) "mox_comparisson;"



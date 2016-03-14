import pandas as pd

all_genes = pd.read_csv('../Data/RNASeq/RNAseq_all_genes.csv', sep=',', index_col=0, header=False)
all_genes.index.name = 'Entrez_ID'
effector_genes = pd.read_csv('../Data/RNASeq/rnaseq_effector_path_vals.csv', sep=',', index_col=0, header=False)
mod_activities = pd.read_csv('../Data/RNASeq/rnaseq_metabolic_mod_activities.csv', sep=',', index_col=0, header=False)
mod_genevalues = pd.read_csv('../Data/RNASeq/rnaseq_metabolic_mod_genevalues.csv', sep=',', index_col=0, header=False)
mod_nodevalues = pd.read_csv('../Data/RNASeq/rnaseq_metabolic_mod_nodevalues.csv', sep=',', index_col=0, header=False)
signaling_genes = pd.read_csv('../Data/RNASeq/rnaseq_signaling_genes.csv', sep=',', index_col=0, header=False)
design = pd.read_csv('../Data/RNASeq/rnaseq_design.csv', index_col=1)
design.columns = [u'idx', u'Drug/Chemical', u'Effect', u'Set']
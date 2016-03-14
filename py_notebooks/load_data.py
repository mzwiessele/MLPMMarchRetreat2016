import pandas as pd, sys, os

dtype = sys.argv[1]

if dtype.lower() in 'rnaseq':
    print "loading rnaseq data"
    rnaseq_all_genes = pd.read_csv('../Data/RNASeq/RNAseq_all_genes.csv', sep=',', index_col=0, header=0)
    rnaseq_all_genes.index.name = 'Entrez_ID'
    rnaseq_effector_genes = pd.read_csv('../Data/RNASeq/rnaseq_effector_path_vals.csv', sep=',', index_col=0, header=0)
    rnaseq_mod_activities = pd.read_csv('../Data/RNASeq/rnaseq_metabolic_mod_activities.csv', sep=',', index_col=0, header=0)
    rnaseq_mod_genevalues = pd.read_csv('../Data/RNASeq/rnaseq_metabolic_mod_genevalues.csv', sep=',', index_col=0, header=0)
    rnaseq_mod_nodevalues = pd.read_csv('../Data/RNASeq/rnaseq_metabolic_mod_nodevalues.csv', sep=',', index_col=0, header=0)
    rnaseq_signaling_genes = pd.read_csv('../Data/RNASeq/rnaseq_signaling_genes.csv', sep=',', index_col=0, header=0)
    rnaseq_design = pd.read_csv('../Data/RNASeq/rnaseq_design.csv', index_col=1)
    rnaseq_design.columns = [u'idx', u'Drug/Chemical', u'Effect', u'Set']
    
    print "rnaseq_all_genes:", rnaseq_all_genes.shape
    print "rnaseq_effector_genes:", rnaseq_effector_genes.shape
    print "rnaseq_mod_activities:", rnaseq_mod_activities.shape
    print "rnaseq_mod_genevalues:", rnaseq_mod_genevalues.shape
    print 'rnaseq_mod_nodevalues:', rnaseq_mod_nodevalues.shape
    print 'rnaseq_signaling_genes:', rnaseq_signaling_genes.shape
    print 'rnaseq_design:', rnaseq_design.shape
else:
    print "loading microarray data"
    micro_all_genes = pd.read_csv('../Data/microarray/microarray_all_genes.csv', sep=',', index_col=0, header=0)
    micro_all_genes.index.name = 'Entrez_ID'
    micro_effector_genes = pd.read_csv('../Data/microarray/microarray_effector_path_vals.csv', sep=',', index_col=0, header=0)
    micro_mod_activities = pd.read_csv('../Data/microarray/microarray_metabolic_mod_activities.csv', sep=',', index_col=0, header=0)
    micro_mod_genevalues = pd.read_csv('../Data/microarray/microarray_metabolic_mod_genevalues.csv', sep=',', index_col=0, header=0)
    micro_mod_nodevalues = pd.read_csv('../Data/microarray/microarray_metabolic_mod_nodevalues.csv', sep=',', index_col=0, header=0)
    micro_signaling_genes = pd.read_csv('../Data/microarray/microarray_signaling_genes.csv', sep=',', index_col=0, header=0)
    micro_design = pd.read_csv('../Data/microarray/microarray_design.csv', index_col=1)
    micro_design.columns = [u'idx', u'Drug/Chemical', u'Effect', u'Set']


    print "micro_all_genes:", micro_all_genes.shape
    print "micro_effector_genes:", micro_effector_genes.shape
    print "micro_mod_activities:", micro_mod_activities.shape
    print "micro_mod_genevalues:", micro_mod_genevalues.shape
    print 'micro_mod_nodevalues:', micro_mod_nodevalues.shape
    print 'micro_signaling_genes:', micro_signaling_genes.shape
    print 'micro_design:', micro_design.shape
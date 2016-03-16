import sys

execfile('load_data.py')

all_data = [
    rnaseq_all_genes,#: (11216, 104)
    rnaseq_effector_genes,#: (1044, 104)
    rnaseq_mod_activities,#: (89, 104)
    rnaseq_mod_genevalues,#: (172, 104)
    rnaseq_mod_nodevalues,#: (462, 104)
    rnaseq_signaling_genes,#: (2184, 104)
#    rnaseq_design: (104, 4)
]

names = ['all_genes', 
         'effector_genes', 
         'mod_activities', 
         'mod_genevalues', 
         'mod_nodevalues',
         'signaling_genes'
         ]

for y,n in zip(all_data, names):
    y.name = n

design = rnaseq_design
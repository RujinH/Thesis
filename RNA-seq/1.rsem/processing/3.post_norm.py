#coding=utf-8
"""
post edgeR clean-up and annotation.

"""

import sys, os, glob
from glbase import *



raw_expn = expression(filename="output2/rawtags_gc_normed.tsv", format={"force_tsv": True, "skiplines": -1, "ensg": 0}, expn="column[1:]")
raw_expn.sort_conditions()
print raw_expn.getConditionNames()


user_path = os.path.expanduser("~")
ensg = glload(os.path.join(user_path, "documents/projects/rnaseq/te_corep_andrew/mm10", "mm10_ensembl_v79_enst.glb")) 
ensg = ensg.removeDuplicates("ensg")
mm=ensg.map(key='ensg',genelist=raw_expn)
mm.saveTSV('output3/genes_ntc_expression.tsv')

import sys, os, glob, re, time

from glbase import *

all = None

conds = []

bad_samples = [] 

for f in glob.glob("data/*.genes.results"):
    head = os.path.split(f)[1]	
    base=head
    #emphasis_pattern=re.compile(r"fq\.exp.+?\.genes.results",re.VERBOSE)
    #base=re.sub(emphasis_pattern,'',base)
    base = base.replace('Mm_', '').replace('ESC_','').replace('.genes.results','')
    #print base

    conds.append(base)

    print "...", base
    
    # skip the glbase overhead:
    rsem = []
    oh = open(f, "rU")
    for line in oh:
        ll = line.strip().split("\t")
        if not "gene_id" in ll[0]:
            rsem.append({"ensg": ll[0], base: float(ll[4])})
    oh.close()
    
    
    if not all:
        all = rsem
    else:   
        for idx, g in enumerate(all):
            if g["ensg"] == rsem[idx]["ensg"]:
                g[base] = rsem[idx][base]
            else:
                print f
                print base
                print rsem[idx]["ensg"]
                print "Oh dear, not cool"
                sys.exit()
        #all = all.map(key="ensg", genelist=rsem)

cond_names = conds

expn = expression(loadable_list=all, expn=cond_names)
expn.coerce(int)
expn.saveTSV("output1/norm_input.tsv")
#expn.save("output1/norm_input.glb")



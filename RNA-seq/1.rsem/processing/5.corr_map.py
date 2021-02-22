
import sys, os
from glbase import *

config.draw_mode = "png"

arr = expression(filename="output3/genes_ntc_expression.tsv",format={'force_tsv':True,'name':2},expn='column[5:]')
arr=arr.filter_low_expressed(1,1)
#arr.log(2, 0.1)

l=len(arr.getConditionNames())
print l

# pretty arr names
#arr.setConditionNames([sam_map.sample_description[i] for i in arr.getConditionNames()])

#arr.boxplot(filename="output4/genes_box.png", size=(45,5))
#arr.hist("output4/genes_hist.png", suppress_zeros=False)
arr.correlation_heatmap(filename="output4/without_log2_heatmap.png", 
    size=(12,10),imshow=True,heat_wid=0.6, heat_hei=0.7, col_font_size=5,row_font_size=5)



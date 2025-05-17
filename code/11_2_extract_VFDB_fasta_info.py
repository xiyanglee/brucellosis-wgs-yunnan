import re
import pandas as pd

## https://zhuanlan.zhihu.com/p/651567616
# zcat VFDB_setA_pro.fas.gz | grep '^>' > SetA_anno.txt
# zcat VFDB_setB_pro.fas.gz | grep '^>' > SetB_anno.txt

# python 11_2_extract_VFDB_fasta_info.py

# setA
with open('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly/source/VFDB/SetA_anno.txt',"r") as f:
    data=f.readlines()
anno_list=[]
for line in data:
    info = re.findall('>(\S+)\(gb\|\S+\) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    if len(info)==0:
        info = re.findall('>(\S+) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    tmp = pd.DataFrame(info)
    tmp.columns=['VF_gene_id','VF_gene_symbol','Gene_description','VF_name','VF_id','VF_category_level1','VF_category_id','taxonomy']
    anno_list.append(tmp)
anno = pd.concat(anno_list)
f.close()
anno.to_csv('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly/source/VFDB/SetA_info.tsv',sep='\t',index=0)

#setB
import re
import pandas as pd
with open('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly/source/VFDB/SetB_anno.txt',"r") as f:
    data=f.readlines()
anno_list=[]
for line in data:
    info = re.findall('>(\S+)\(gb\|\S+\) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    if len(info)==0:
        info = re.findall('>(\S+) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    tmp = pd.DataFrame(info)
    tmp.columns=['VF_gene_id','VF_gene_symbol','Gene_description','VF_name','VF_id','VF_category_level1','VF_category_id','taxonomy']
    anno_list.append(tmp)
anno = pd.concat(anno_list)
f.close()
anno.to_csv('/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024/pipelines/denovo_ssembly/source/VFDB/SetB_info.tsv',sep='\t',index=0) 
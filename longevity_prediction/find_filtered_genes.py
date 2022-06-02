from goatools import obo_parser
import os
import wget
import pandas as pd
import re
from Bio.UniProt.GOA import gafiterator
import gzip

go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
folder = os.getcwd()

if(not os.path.isfile('./go-basic.obo')):
    go_obo = wget.download(go_obo_url, 'go-basic.obo')
else:
    go_obo = 'go-basic.obo'

go = obo_parser.GODag(go_obo)

relevant_go_terms = set()

regex = re.compile(r"\b(aging|longevity|lifespan|senescence|stress|ageing|cell death|age)\b", re.IGNORECASE)

for go_term in go.values():
    if (regex.search(go_term.name)):
        relevant_go_terms.add(go_term.id)

all_children = []

for term in relevant_go_terms:
    all_children.append(go[term].get_all_children())

for children in all_children:
    for child in children:
        relevant_go_terms.add(child)

filename = './data/datastore/wb.gaf.gz'

relevant_genes = set()
with gzip.open(filename, 'rt') as fp:
    for annotation in gafiterator(fp):
        if (annotation['GO_ID'] in relevant_go_terms
            and annotation['DB_Object_ID'][0:6] == "WBGene"):
                relevant_genes.add(annotation['DB_Object_ID'])

gene_expression_data = pd.read_csv("../common_datastore/raw_gene_expression_no_outliers_and_singles_for_GO_filtering.csv")

gene_expression_data.columns = ['Sample'] + list(gene_expression_data.columns[1:])
gene_expression_data.set_index("Sample", inplace=True)

relevant_genes = [gene for gene in relevant_genes if gene in gene_expression_data.index]

filtered = gene_expression_data.loc[relevant_genes, :]

filtered.to_csv("../common_datastore/GO_filtered_raw_gene_expression_data.csv")
import pandas as pd
from io import StringIO
import subprocess
import pickle
import os

with open('genage_all_genes.txt', 'r') as f:
    GENAGE_GENES = f.read().split("\n")
    GENAGE_GENES.pop()
    GENAGE_GENES = list(set(GENAGE_GENES))

def create_metadata_df(metadata):
    dfs_concatenated = pd.DataFrame()
    for data in metadata:
        try:
            df = pd.read_csv(StringIO(data.stdout))
            dfs_concatenated = pd.concat([dfs_concatenated, df])
        except pd.errors.EmptyDataError as _:
            pass
        except:
            print('There are other unexpected errors. Exiting.')
            exit()
    return dfs_concatenated


def fetch_metadata(list_of_queries, by_bioproject_id=True):
    if by_bioproject_id:
        def single_metadata_fetcher(bioproject):
            return subprocess.run(f"esearch -db sra -query '{bioproject}[bioproject]' | efetch -format runinfo",
                                  shell=True, capture_output=True, text=True)
    else:
        def single_metadata_fetcher(query):
            return subprocess.run(f"esearch -db gds -query \"{query}\" | efilter -query "
                                  "\"expression profiling by high throughput sequencing "
                                  "[DataSet Type] AND Caenorhabditis elegans [organism] AND 1900/01:2020/08[Publication Date] \" " 
                                  "| elink -target bioproject | elink -target sra | efetch -format runinfo"
                                  , shell=True, capture_output=True, text=True)
    metadata = []
    couldnt_fetch = []
    error_code = '429'

    for i in range(0, len(list_of_queries)):
        metadata.append(single_metadata_fetcher(list_of_queries[i]))

        if error_code in metadata[i].stderr:
            couldnt_fetch.append(list_of_queries[i])

    return metadata, couldnt_fetch

keyword_queries = ['lifespan', 'aging', 'life span', 'longevity', 'senescence', 'ageing',
                   'longlived', 'shortlived', 'long-lived', 'short-lived', 'stress']

if (not os.path.exists('keywords_metadata.p')):
    keywords_metadata, couldnt_fetch_keyword_queries = fetch_metadata(keyword_queries, by_bioproject_id = False)
    pickle.dump( keywords_metadata, open( "keywords_metadata.p", "wb" ) )
    pickle.dump( couldnt_fetch_keyword_queries, open( "couldnt_fetch_keyword_queries.p", "wb" ) )

if (not os.path.exists('genage_genes_metadata.p')):
   genage_genes_metadata, couldnt_fetch_genage = fetch_metadata(GENAGE_GENES, by_bioproject_id=False)
   pickle.dump( couldnt_fetch_genage, open( "couldnt_fetch_genage.p", "wb" ) )
   pickle.dump( genage_genes_metadata, open( "genage_genes_metadata.p", "wb" ) )

# Add any study that you know about that wasn't retrieved
if (not os.path.exists('manually_fetched_metadata.p')):
    manually_fetched_metadata, couldnt_fetch_manual = fetch_metadata(['PRJNA261420'])
    pickle.dump( couldnt_fetch_manual, open( "couldnt_fetch_manual.p", "wb" ) )
    pickle.dump( manually_fetched_metadata, open( "manually_fetched_metadata.p", "wb" ) )

keywords_metadata = pickle.load(open( "keywords_metadata.p", "rb" ))
genage_genes_metadata = pickle.load(open( "genage_genes_metadata.p", "rb" ))
manually_fetched_metadata = pickle.load(open( "manually_fetched_metadata.p", "rb" ))

keywords_metadata_df = create_metadata_df(keywords_metadata)
genage_genes_metadata_df = create_metadata_df(genage_genes_metadata)
manually_fetched_metadata_df = create_metadata_df(manually_fetched_metadata)

combined_dfs_raw = pd.concat([manually_fetched_metadata_df, genage_genes_metadata_df, keywords_metadata_df])
combined_dfs_processed = combined_dfs_raw.copy()
# Remove duplicate runs
combined_dfs_processed = combined_dfs_processed.drop_duplicates('Run')
combined_dfs_processed.set_index('Run', inplace=True)
combined_dfs_processed = combined_dfs_processed.drop(
    (combined_dfs_processed.loc[combined_dfs_processed['LibraryStrategy'] != 'RNA-Seq']).index)

combined_dfs_processed = combined_dfs_processed.loc[(combined_dfs_processed['ScientificName'] == 'Caenorhabditis elegans').values]

if (not os.path.exists("all_bioprojects.txt")):
    with open('all_bioprojects.txt', 'w') as f:
        for index, row in combined_dfs_processed.iterrows():
            f.write(row['BioProject'] + '\n')

# A manual check to retrieve relevant studies (e.g look at abstract, intro, figures, etc...) is inevitable at this point
relevant_bioprojects = ["PRJNA171150",
                        "PRJNA176153",
                        "PRJNA198476",
                        "PRJNA198785",
                        "PRJNA241944",
                        "PRJNA253943",
                        "PRJNA257404",
                        "PRJNA259222",
                        "PRJNA261420",
                        "PRJNA261970",
                        "PRJNA266552",
                        "PRJNA268148",
                        "PRJNA273888",
                        "PRJNA274976",
                        "PRJNA282784",
                        "PRJNA285045",
                        "PRJNA285208",
                        "PRJNA287674",
                        "PRJNA294011",
                        "PRJNA297365",
                        "PRJNA309116",
                        "PRJNA312237",
                        "PRJNA313027",
                        "PRJNA314583",
                        "PRJNA315807",
                        "PRJNA325325",
                        "PRJNA340183",
                        "PRJNA342039",
                        "PRJNA343462",
                        "PRJNA356824",
                        "PRJNA362200",
                        "PRJNA376493",
                        "PRJNA391939",
                        "PRJNA395909",
                        "PRJNA414453",
                        "PRJNA419144",
                        "PRJNA421426",
                        "PRJNA423001",
                        "PRJNA423029",
                        "PRJNA435494",
                        "PRJNA436362",
                        "PRJNA436651",
                        "PRJNA449675",
                        "PRJNA471520",
                        "PRJNA473811",
                        "PRJNA475273",
                        "PRJNA484793",
                        "PRJNA489378",
                        "PRJNA490034",
                        "PRJNA491191",
                        "PRJNA497519",
                        "PRJNA506829",
                        "PRJNA506844",
                        "PRJNA507640",
                        "PRJNA507749",
                        "PRJNA507775",
                        "PRJNA509132",
                        "PRJNA510076",
                        "PRJNA514750",
                        "PRJNA520897",
                        "PRJNA522776",
                        "PRJNA525790",
                        "PRJNA528111",
                        "PRJNA549074",
                        "PRJNA554266",
                        "PRJNA564575",
                        "PRJNA574273",
                        "PRJNA574844",
                        "PRJNA575569",
                        "PRJNA576205",
                        "PRJNA592073",
                        "PRJNA602828",
                        "PRJNA604568",
                        "PRJNA610525"]

# Manual work done to obtain labels.csv and age.csv
# NOTE: labels.csv and age.csv are not needed until we build models. Hence, we put them in
# the directory longevity_prediction/data/datastore/raw/

# The longevity code for labels.csv is as follows:
# long-lived: 0
# normal: 1
# short-lived: 2

# Additionally, for age.csv, we use the following scheme:
# day 1 adults are labeled as 1, day 2 as 2, as so on
# samples for which there is no specific age but are instead labeled as "young adults" are considered to be day 1 adults
# samples for which there is no specific age but are instead labeled as "prefertile young adults" are considered to be day 0.5 adults
# L4 worms are labeled as 0
# L3 worms as -0.5
# L2 worms as -1
# L1 worms as -2
# embryos as -2.5
# samples that have missing age info are considered day 1 adults since it is the mode of the age attribute

# We throw out samples for which we were not able to find longevity information
samples_to_throw_out = ["SRR1578746",
                        "SRR1578747",
                        "SRR6761043",
                        "SRR6761044",
                        "SRR6761045",
                        "SRR8079476",
                        "SRR8079477",
                        "SRR2537188",
                        "SRR2537187",
                        "SRR2537190",
                        "SRR2537189",
                        "SRR2537200",
                        "SRR2537199",
                        "SRR2537202",
                        "SRR2537201",
                        "SRR2537204",
                        "SRR2537203",
                        "SRR2537206",
                        "SRR2537205",
                        "SRR2005821",
                        "SRR11005177",
                        "SRR11005178",
                        "SRR11005179",
                        "SRR11005180",
                        "SRR11005181",
                        "SRR11005182",
                        "SRR11005183",
                        "SRR11005184",
                        "SRR6792637",
                        "SRR6792638",
                        "SRR6792639",
                        "SRR6792640",
                        "SRR6792641",
                        "SRR6792642",
                        "SRR10544936",
                        "SRR10544937",
                        "SRR10544938",
                        "SRR8528174",
                        "SRR8528175",
                        "SRR8528176",
                        "SRR8528177",
                        "SRR2043531",
                        "SRR2043532"]

with open('accession_list.txt', 'w') as f:
    for index, row in combined_dfs_processed.iterrows():
        if (row['BioProject'] in relevant_bioprojects and index not in samples_to_throw_out):
            f.write(index + " " + row['BioProject'] + " " + row['LibraryLayout'].upper() + " " + '\n')

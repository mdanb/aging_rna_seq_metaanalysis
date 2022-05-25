import os
import pandas as pd
from data.utils import ensembl_protein_to_gene_map
import pickle

class GeneInteractionGraph(object):
    def __init__(self):
        self.load_data()

    def load_data(self):
        raise NotImplementedError

class StringDBGraph(GeneInteractionGraph):
    """
    Download link : https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Caenorhabditis+elegans
    """

    def __init__(self):
        self.proteinlinks = "data/datastore/6239.protein.links.detailed.v11.5.txt" # C. elegans StringDB PPI data download from link ^^^
        self.edge_weights = []
        self.graph_edgelist = []

        super(StringDBGraph, self).__init__()

    def load_data(self):
        edgelist_savefile = f"./data/datastore/stringdb_graph_coexpression.p"
        edgeweights_savefile = f"./data/datastore/stringdb_edge_weights_coexpression.p"

        if os.path.isfile(edgelist_savefile) and os.path.isfile(edgeweights_savefile):
            print(f" loading from cache file {edgelist_savefile} and {edgeweights_savefile}")
            self.graph_edgelist = pickle.load(open( f"{edgelist_savefile}", "rb" ))
            self.edge_weights = pickle.load(open( f"{edgeweights_savefile}", "rb" ))

        else:
            print("Building StringDB Graph. It can take a second the first time")
            print("ensp_to_ensg_map")
            ensmap = ensembl_protein_to_gene_map()
            print(" reading self.proteinlinks")
            edges = pd.read_csv(self.proteinlinks, sep=' ')

            selected_edges = edges['coexpression'] != 0

            stringDB_edgelist = edges[selected_edges][["protein1", "protein2", 'coexpression']].values.tolist()
            for edge in stringDB_edgelist:
                for i, node in enumerate(edge[:-1]):
                    if (node.count(".") != 3):
                        edge[i] += ".1"
                if edge[0][5:] in ensmap.keys() and edge[1][5:] in ensmap.keys():
                    self.graph_edgelist.append([ensmap[edge[0][5:]], ensmap[edge[1][5:]]])
                    self.edge_weights.append(edge[2])
            print(" writing gene interaction graph edgelist and edge weights")
            pickle.dump(self.graph_edgelist, open(f"{edgelist_savefile}", "wb"))
            pickle.dump(self.edge_weights, open(f"{edgeweights_savefile}", "wb"))

            print("Graph built!")


# A Meta-Analysis of RNA-Seq Gene Expression Aging Studies
### Mohamad Daniel Bairakdar, Matthias Truttmann, Ambuj Tewari

ADD DOCKER INSTRUCTIONS TO SAVE STUFF SO NOT LOSE STUFF

This repository houses the codebase for the paper **A Meta-Analysis of RNA-Seq Gene Expression Aging Studies**. The aim of the study is to attempt to discover longevity genes/gene modules of C. elegans purely from their RNA-seq gene expression profiles obtained from many different aging studies. 
Here, we present the pipeline to replicate the results of the paper. 

Note that the only piece of software that you'll need to have installed on your machine is Docker. If you're not familiar with Docker and how to use it, fear not. You don't need to know how to use Docker to follow the steps outlined here. All you need to know is that Docker helps you replicate our results with minimum effort. In fact, Docker is the reason why you don't need any packages - other than Docker - installed on your machine: it provides "a software environment" identical to ours. For details on installation and setup, please visit the [official Docker installation instructions](https://docs.docker.com/get-docker/). 

Once you have Docker installed, if you are on Linux, run `sh docker_gpu_setup.sh`, which will set up Docker so that we can access the local machine's GPU. If your machine does not have a GPU, that's fine for everything done in this repo except for training some of the machine learning models. In theory, you could run without a GPU, but in practice, it'll take a very long time to run, especially for the GNN models. So if you don't have access to a GPU, you can skip this step. 

Next, run `docker pull mdanb/aging_rna_seq_metaanalysis` to pull the Docker image associated with this project. Next, start up the software environment (the "container" in Docker speak):
```
docker run -it -e NCBI_API_KEY=<NCBI_KEY> --gpus all aging_rna_seq_metaanalysis
```
where in place of <NCBI_KEY>, you should type in your NCBI API key. Note that using an API key is not strictly necessary, but is highly recommended if you actually want to query the databases and not use the outputs generated by us (see note at the end of this file), as it allows you to post more requests per second to the NCBI databases than if you don't use one. You can read more about NCBI API keys and how to obtain one in the [eutilities guide](https://www.ncbi.nlm.nih.gov/books/NBK179288/). If you choose not to use an API key, leave out `-e NCBI_API_KEY=<NCBI_KEY>`. Also, if you don't have access to a GPU, you can leave out `--gpus all`. 

You are now inside the container, with all necessary dependencies installed (if you know what a virtual machine is, then you can think of the container as a lightweight version of a virtual machine). You'll also notice that the filesystem is separate from your actual filesystem (run `cd ~` and `ls` and you'll notice that you don't see the files you expected to see in your home directory). To exit the container, run `CTRL + D`. 

Inside the container, in the directory `/root`, you'll notice that there are 7 directories. One of them (`common_datastore`) is just a directory to store some data, while 5 of these correspond to the "high-level" tasks that we need to perform:

1. Query (the NCBI databases), download, and process the gene expression data.
2. Exploratory Data Analysis + Data normalization/batch correction
3. (Optional) SAUCIE batch correction
4. (Optional) Build machine learning models to predict the longevity class (long, short, normal lived).
5. (Optional) Perform a traditional differential expression analysis to discover longevity gene/gene modules, bypassing the predictive modeling step. 

The directories `query_download_process_data`, `EDA`, `SAUCIE`, `longevity_prediction`, and `differential_expression_analysis` in this Github repository - which are pretty much equivalent to the directories you see in the container - contain READMEs with instructions for the above mentioned tasks, respectively. The `paper` directory just contains some information used to write the paper.

**Note:** All the output files associated with each step of the workflow are already inside the container in the relevant directories. This way, you do not actually have to run the code (which can be quite time consuming at certain points) to obtain its output, but can instead "pick up" the workflow at any stage you wish.

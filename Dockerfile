FROM pytorch/pytorch:1.8.0-cuda11.1-cudnn8-runtime

WORKDIR /root

SHELL ["/bin/bash", "-c"] 

RUN apt-get update && apt-get install -y \
		git \
		cmake \
		build-essential \
		zlib1g-dev \
		autoconf \
		libhdf5-dev \
		perl \
		curl

RUN mkdir .ncbi
RUN printf $'/LIBS/GUID = "64bcc153-f74b-4505-937d-673f388d0e76"\n\
/config/default = "false"\n\
/libs/cloud/report_instance_identity = "true"\n\
/repository/user/ad/public/apps/file/volumes/flatAd = "."\n\
/repository/user/ad/public/apps/refseq/volumes/refseqAd = "."\n\
/repository/user/ad/public/apps/sra/volumes/sraAd = "."\n\
/repository/user/ad/public/apps/sraPileup/volumes/ad = "."\n\
/repository/user/ad/public/apps/sraRealign/volumes/ad = "."\n\
/repository/user/ad/public/apps/wgs/volumes/wgsAd = "."\n\
/repository/user/ad/public/root = "."\n\
/repository/user/default-path = "/ncbi"\n\
/repository/user/main/public/root = "."' > .ncbi/user-settings.mkfg

RUN mkdir query_download_process_data
WORKDIR /root/query_download_process_data

RUN echo "export PATH=$PATH:$PWD/sratoolkit.2.11.0-ubuntu64/bin" >> ../.bashrc && source ../.bashrc

RUN git clone https://github.com/pachterlab/kallisto.git
RUN cd kallisto && git reset --hard ae81a86 && cd ext/htslib && autoheader && autoconf && cd ../.. \
    && mkdir build && cd build && cmake .. && make && make install

#this is a specific version of the toolkit
COPY query_download_process_data/ ./ 
RUN tar -xvzf sratoolkit.3.0.0-ubuntu64.tar.gz
RUN export PATH=${PATH}:/home/mdanb/sratoolkit.3.0.0-ubuntu64/bin
ENV NCBI_API_KEY=''

RUN apt-get update && apt-get install -y \
   inotify-tools \
   pigz \
   libxml2-dev \
   libcurl4-openssl-dev \
   libssl-dev \
   ncbi-entrez-direct\
   libbz2-dev \
   liblzma-dev \
   gfortran \
   liblapack-dev \
   liblapack3 \
   libopenblas-base \
   libopenblas-dev

RUN apt install -y --no-install-recommends software-properties-common dirmngr 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt install -y -q --no-install-recommends r-base
#=4.1.0-1.2004.0

#RUN pip3 install torch-scatter -f https://pytorch-geometric.com/whl/torch-1.8.0+cu111.html
#RUN pip3 install torch-sparse -f https://pytorch-geometric.com/whl/torch-1.8.0+cu111.html
#RUN pip3 install torch-geometric==2.0.1

RUN pip3 install pandas==1.2.4 tqdm==4.62 scikit-learn==0.24.1 
RUN pip3 install xgboost==1.4.1 
RUN pip3 install pipelinehelper==0.7.8 joblib==1.0.1
RUN pip3 install gtfparse==1.2.1

RUN apt-get install vim -y

RUN R -e "install.packages('BiocManager', version='3.15')"
RUN R -e "BiocManager::install('sva')"
RUN R -e "BiocManager::install('biomaRt')"                                          
RUN R -e "BiocManager::install('edgeR')"                                          
RUN R -e "BiocManager::install('ensembldb')"
RUN R -e "BiocManager::install('tximport')"                                          
RUN R -e "install.packages('tidyverse')"
RUN R -e "BiocManager::install('limma')"
RUN R -e "BiocManager::install('msigdbr')"                                        
RUN R -e "BiocManager::install('gage')" 
RUN R -e "BiocManager::install('fgsea')"
RUN R -e "install.packages('data.table')"                                         

WORKDIR /root
COPY differential_expression_analysis differential_expression_analysis/                                            
COPY EDA EDA/                                        
COPY longevity_prediction longevity_prediction/
COPY SAUCIE SAUCIE/                                        

RUN pip3 install tensorflow==2.9.0
RUN apt-get install -y parallel

FROM continuumio/miniconda3


RUN conda install -c bioconda clustalo=1.2.3 -y 
RUN conda install -c bioconda hmmer=3.3.2 -y
RUN conda install -c bioconda blast=2.5.0 -y
RUN conda install -c bioconda jolytree=1.1b -y
RUN conda install -c bioconda blast-legacy=2.2.26 -y
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/ && \
 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/pro/

RUN wget https://dl.secondarymetabolites.org/releases/latest/antismash-6.0.1.tar.gz
RUN tar -zxvf antismash-6.0.1.tar.gz
RUN pip3 install ./antismash-6.0.1
RUN conda install -c bioconda hmmer2=2.3.2 diamond=0.9.14 fasttree=2.1.10 prodigal=2.6.3 muscle=3.8.1551 glimmerhmm=3.0.4 -y
RUN download-antismash-databases

RUN apt-get update
RUN apt-get install vim -y


RUN pip3 install biopython==1.78
RUN pip3 install pandas==1.3.2
RUN pip3 install openpyxl==3.0.7
RUN pip3 install django-tastypie
RUN pip3 install django
RUN pip3 install ncbi-genome-download

RUN mkdir /opt/DockerImage
COPY  main.py test.json KeywordsSearch.py Database/ /opt/DockerImage/

CMD ["bash"]

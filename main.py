from posixpath import basename
from Bio import Entrez
from xml.dom.minidom import parse
from Bio.Align import MultipleSeqAlignment
from Bio.Entrez.Parser import ValidationError
from Bio.Entrez.Parser import CorruptedXMLError
import requests
import io
import pandas as pd
from pandas.errors import EmptyDataError
import re
import os
from copy import deepcopy
import subprocess
import sys
from urllib.error import HTTPError
from urllib.error import URLError
import time
import ftplib
from io import BytesIO
import ftplib
import gzip
import zlib
import json
import argparse
from numpy import nan
from Bio import SeqIO
import glob
import django
import tastypie
from tastypie import http
from urllib.request import urlopen 
from socket import timeout
import http
from http import client
from http.client import IncompleteRead
import Bio
import matplotlib.pyplot as plt



########### arguments ###########
parser = argparse.ArgumentParser(description="BGCSniffer")
parser.add_argument("-json", "--jsonfile", dest="jsonfile", help= "Please input the complete path of parameter file (JSON format)")
parser.add_argument("-workdir", "--workdir", dest="workdir", help= "Please input your working directory")
args = parser.parse_args()
########### arguments ###########

########### functions ###########
def readPara(Parafile):
    with open (Parafile) as Para:
        parameters = json.load(Para)
        return(parameters)


def extract_query_info(query_string):
    id_list = []
    if ('and' in query_string):
        query_list = re.split(' and ', query_string)
        for parenthesis_ele in query_list:
            parenthesis_ele = parenthesis_ele.replace('(', '')
            parenthesis_ele = parenthesis_ele.replace(')', '')
            if 'AND' in parenthesis_ele:
                for id in re.split(' AND ', parenthesis_ele):
                    id_list.append(id)
            else:
                for id in re.split(' OR ', parenthesis_ele):
                    id_list.append(id)
    elif ('or' in query_string):
        query_list = re.split(' or ', query_string)
        for parenthesis_ele in query_list:
            parenthesis_ele = parenthesis_ele.replace('(', '')
            parenthesis_ele = parenthesis_ele.replace(')', '')
            if 'AND' in parenthesis_ele:
                for id in re.split(' AND ', parenthesis_ele):
                    id_list.append(id)
            else:
                for id in re.split(' OR ', parenthesis_ele):
                    id_list.append(id)
    else:
        query_list = re.split(' or ', query_string)
        for parenthesis_ele in query_list:
            parenthesis_ele = parenthesis_ele.replace('(', '')
            parenthesis_ele = parenthesis_ele.replace(')', '')
            if 'AND' in parenthesis_ele:
                for id in re.split(' AND ', parenthesis_ele):
                    id_list.append(id)
            else:
                for id in re.split(' OR ', parenthesis_ele):
                    id_list.append(id)
    return(id_list)

def searchUniprot(pf_id):
    BASE = 'http://www.uniprot.org'
    KB_ENDPOINT = '/uniprot/'
    QUERY_ID = pf_id
    UNIPROT_URL = BASE + KB_ENDPOINT + '?query=' + QUERY_ID +  '&sort=score&columns=id,entry name,reviewed,protein names,genes,organism,organism-id,length,sequence&format=tab'
    uniprot_result = requests.get(UNIPROT_URL)

    if uniprot_result.ok:
        content = uniprot_result.text
        content = io.StringIO(content)
        uniprot_df = pd.read_csv(content, sep="\t", index_col=0)
        uniprot_df['ProteinFamily'] = pf_id
        return(uniprot_df)
    else:
        return('Something went wrong in Interpro searching ' +  uniprot_result.status_code)

def findMutualProtein(queryString, combined_df):
    if re.search('and', queryString):
        queryString_list = re.split(' and ', queryString)
        mutual_species_id = set(combined_df['Organism ID'].to_list())
        for block in queryString_list:
            block = block.replace('(', '')
            block = block.replace(')', '')
            tmp_block_list = []
            for element in re.split(' OR ', block):
                tmp_speciesID = combined_df[combined_df['ProteinFamily'] == element]['Organism ID'].to_list()
                tmp_block_list.extend(tmp_speciesID)
            tmp_block_list = list(set(tmp_block_list))
            mutual_species_id = mutual_species_id.intersection(tmp_block_list)
        filtered_df = combined_df[combined_df['Organism ID'].isin(list(mutual_species_id)) & ~combined_df['Organism'].isin(['uncultured bacterium','bioreactor metagenome'])]
    else:
        filtered_df = combined_df
    return(filtered_df)

def downloadGenome(dataframe, work_dir):
    species_id_list = []
    organism_list = []
    protein_family_info = {}
    all_assembly_id_list = []

    for row in range(0, dataframe.shape[0]):
        species_name = dataframe['Organism'].iloc[row,]
        try:
            genus = re.split(' ', species_name)[0]
            species = re.split(' ', species_name)[1]
            #query_name = ' '.join([genus,species])
            protein_ID = dataframe.iloc[row,].name
            #protein_sequence = merged_df['Sequence'].iloc[row,]
            organism = dataframe['Organism'].iloc[row,]
            organism_id = dataframe['Organism ID'].iloc[row,]
            ## Get protein family information for each protein
            proteinfamily_type = dataframe['ProteinFamily'].iloc[row,]

            protein_family_info[protein_ID] = proteinfamily_type
            species_id_list.append(organism_id)
            organism_list.append(organism)
        except IndexError:
            next

    # Retrieving genomes from the Entrez databases
    if os.path.exists(work_dir +'/Candidate_genomes/'): 
        subprocess.run(["rm -rf", work_dir + "/Candidate_genomes/"], shell=True)
        subprocess.run(["mkdir", work_dir + "/Candidate_genomes/"], shell=True)
    else:
        pass

    for tax_id in list(set(species_id_list)):
        ## 创建下载目录
        if os.path.isdir(work_dir + '/Candidate_genomes/'+ str(tax_id) + '/'):
            for root, subdir , files in os.walk(work_dir + '/Candidate_genomes/'+ str(tax_id) + '/'):
                for gz_file in files:
                    if gz_file.endswith('gz'):
                        print(gz_file)
                        gzfile_path = os.path.join(root, gz_file)
                        unzip_status = subprocess.run(['gunzip', gzfile_path], stdout=subprocess.PIPE)
                        if unzip_status.returncode:
                            print(gzfile_path + ' unzip failed!!')
                            next
                        else:
                            pass
                    else:
                        next
        else:
            os.system('mkdir -p ' + work_dir + '/Candidate_genomes/'+ str(tax_id) + '/')
            ## 下载 genebankfile
            genebank_flag = 0
            status = subprocess.run(['ncbi-genome-download', '--taxids' , str(tax_id), 'bacteria' ,'--assembly-levels', 'complete,chromosome', '--flat-output','--parallel', '10','-r' , '10', '-o' , work_dir + '/Candidate_genomes/'+ str(tax_id) + '/', '-s', 'genbank'], stdout=subprocess.PIPE)
            if status.returncode:
                print('Something went wrong when downloading the genebank file of tax id ' + str(tax_id) + ' !!')
                next
            else:
                print('tax id ' + str(tax_id) + ' download finished!!')
                genebank_flag = 1
                pass

            ## 下载 fastafile
            fasta_flag = 0
            status = subprocess.run(['ncbi-genome-download', '--taxids' , str(tax_id), 'bacteria' ,'--assembly-levels', 'complete,chromosome', '-F' ,'fasta' ,'--flat-output','--parallel', '10', '-r' , '10','-o' , work_dir + '/Candidate_genomes/'+ str(tax_id) + '/'], stdout=subprocess.PIPE)
            if status.returncode:
                print('Something went wrong when downloading the ref seq file of tax id ' + str(tax_id) + ' !!')
                next
            else:
                print('tax id ' + str(tax_id) + ' download finished!!')
                fasta_flag = 1
                pass
            
            ## 检查 genebankfile和fastafile 是否同时下载完成
            if fasta_flag + genebank_flag > 0:
                pass
            else:
                print('Something went wrong when downloading the files of tax id ' + str(tax_id) + ' !!')
                ## remove downloaded files when download is not completed
                if os.path.isdir(work_dir + '/Candidate_genomes/'+ str(tax_id) + '/'):
                    subprocess.run(['rm', '-rf', work_dir + '/Candidate_genomes/'+ str(tax_id) + '/'])
                else:
                    pass
                next

            ## 解压下载压缩文件
            for root, subdir , files in os.walk(work_dir + '/Candidate_genomes/'+ str(tax_id) + '/'):
                for gz_file in files:
                    if gz_file.endswith('gz'):
                        print(gz_file)
                        gzfile_path = os.path.join(root, gz_file)
                        unzip_status = subprocess.run(['gunzip', gzfile_path], stdout=subprocess.PIPE)
                        if unzip_status.returncode:
                            print(gzfile_path + ' unzip failed!!')
                            next
                        else:
                            pass
                    else:
                        next

        # try:
        #     handle = Entrez.efetch(db="taxonomy", id=str(tax_id))
        #     result = Entrez.read(handle)
        # except (ValidationError, RuntimeError, IncompleteRead, ConnectionResetError, timeout, CorruptedXMLError, HTTPError):
        #         next

        # level = result[0]['Rank']
        # tax_name = result[0]['ScientificName']
        # ## 移除unculture的物种
        # if re.search('[uU]ncluture', tax_name):
        #     next
        # else:
        #     pass
        # ## Start ncbi parsing
        # ## Get maximum retnum
        # try:
        #     tmp_handle = Entrez.esearch(db="assembly", term = tax_name, retmode="xml")
        #     tmp_result = Entrez.read(tmp_handle)
        #     max_record_num = tmp_result['Count']
        #     tmp1_handle = Entrez.esearch(db="assembly", term = tax_name, retmode="xml", retmax=max_record_num)
        #     tmp1_result = Entrez.read(tmp1_handle)
        #     assmbly_id_list = tmp1_result['IdList']
        #     print("\n")
        #     print("物种/菌株" + tax_name + '在NCBI数据库中检索到如下Assembly ID：' + '，'.join(assmbly_id_list))
        #     print("\n")
        #     print("正在根据组装水平对物种/菌株" + tax_name + "对应的" + str(len(assmbly_id_list)) + "个Assembly ID进行过滤：")
        #     print("\n")

        #     filtered_assembly_id_list = []
        #     ## Check if only one species/strain is found for tax_id
        #     if len(assmbly_id_list) == 1:
        #         assmbly_id = ''.join(assmbly_id_list)
        #         assembly_id_handle = Entrez.esummary(db="assembly", id=assmbly_id, report="full")
        #         try:
        #             assembly_id_result = Entrez.read(assembly_id_handle)
        #             ## Filter assembly IDs with AssemblyStatus, discard contig level assemblies
        #             if (assembly_id_result['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus'] == 'Chromosome') | (assembly_id_result['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus'] == 'Complete Genome'):
        #                 filtered_assembly_id_list.append(assmbly_id)
        #             else:
        #                 next
        #         except (ValidationError, RuntimeError, IncompleteRead, ConnectionResetError, timeout, CorruptedXMLError, HTTPError):
        #             next
        #     else:
        #         if len(assmbly_id_list) > 100:
        #             for part in list(range(0,len(assmbly_id_list)+1,100)):
        #                 partial_assmbly_id_list = assmbly_id_list[part:part+100]
        #                 query_assmbly_id = ','.join(partial_assmbly_id_list)
        #                 assembly_id_handle = Entrez.esummary(db="assembly", id=query_assmbly_id, report="full") 
        #                 try:
        #                     assembly_id_result = Entrez.read(assembly_id_handle)
        #                     recoed_list = assembly_id_result['DocumentSummarySet']['DocumentSummary']
        #                     for record in recoed_list:
        #                         assmbly_id = record.attributes['uid']
        #                         if (record['AssemblyStatus'] == 'Chromosome') | (record['AssemblyStatus'] == 'Complete Genome'):
        #                             filtered_assembly_id_list.append(assmbly_id)
        #                         else:
        #                             next
        #                 except (ValidationError, RuntimeError, IncompleteRead, ConnectionResetError, timeout, CorruptedXMLError, HTTPError):
        #                     next
        #         else:
        #             query_assmbly_id = ','.join(assmbly_id_list)
        #             assembly_id_handle = Entrez.esummary(db="assembly", id=query_assmbly_id, report="full") 
        #             try:
        #                 assembly_id_result = Entrez.read(assembly_id_handle)
        #                 recoed_list = assembly_id_result['DocumentSummarySet']['DocumentSummary']
        #                 for record in recoed_list:
        #                     assmbly_id = record.attributes['uid']
        #                     if (record['AssemblyStatus'] == 'Chromosome') | (record['AssemblyStatus'] == 'Complete Genome'):
        #                         filtered_assembly_id_list.append(assmbly_id)
        #                     else:
        #                         next
        #             except (ValidationError, RuntimeError, IncompleteRead, ConnectionResetError, timeout, CorruptedXMLError, HTTPError):
        #                 next
        # except (TimeoutError,RuntimeError,ValidationError,IncompleteRead, ConnectionResetError, timeout, CorruptedXMLError, HTTPError):
        #     next
        # print("物种/菌株" + tax_name + '过滤后剩余' + str(len(filtered_assembly_id_list)) + '个ID：' + '，'.join(filtered_assembly_id_list))
        # print("\n")

        # print("\n")
        # print("正在下载物种/菌株" + tax_name + '过滤后的' + str(len(filtered_assembly_id_list)) + '个complete genomes: ')
        # print("\n")

        # filtered_assembly_id_list = list(set(filtered_assembly_id_list))
        # if len(filtered_assembly_id_list) == 0:
        #     next
        # else:
        #     ## 创建下载目录
        #     os.system('mkdir -p ' + work_dir + '/Candidate_genomes/'+ str(tax_id) + '/')
        #     for genome_id in filtered_assembly_id_list:
        #     #如果基因组已经下载过则next
        #         if genome_id in all_assembly_id_list:
        #             next
        #         else:
        #             pass
                    
                # genome_handle = Entrez.esummary(db="assembly", id=genome_id, report="full")
                # try:
                #     genome_record = Entrez.read(genome_handle)
                #     genome_result = genome_record['DocumentSummarySet']['DocumentSummary'][0]
                #     assembly_accession = genome_result['AssemblyAccession']
                #     assembly_ftp = genome_result['FtpPath_GenBank']
                #     organism = genome_result['Organism']
                #     ## 移除unculture的物种
                #     if re.search('[uU]ncluture', organism):
                #         next
                #     else:
                #         pass
                #     try:
                #         strain_info = genome_result['Biosource']['InfraspeciesList'][0]['Sub_value']
                #     except IndexError:
                #         strain_info = 'unknown'
                #     strain_info = re.sub('[/\|\s:]', '_', strain_info)
                #     assembly_name = re.sub('[/\|\s:]', '_', organism)
                #     assembly_name = re.sub('\(.+\)', '', assembly_name)
                #     assembly_name = assembly_name + '_strain_' + strain_info + '_' + assembly_accession
                #     assembly_name = re.sub(';', '', assembly_name)
                #     assembly_name = re.sub('__', '_', assembly_name)
                #     assembly_name = re.sub(':', '_', assembly_name)

                #     assembly_name = re.sub('\[', '', assembly_name)
                #     assembly_name = re.sub('\]', '', assembly_name)
                #     #如果该物种/strain基因组未下载过则存入所有下载过的基因组id
                #     all_assembly_id_list.extend(genome_id)
                #     ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
                #     ## 当某些基因组没有ftp链接，直接舍弃
                #     try:
                #         ftp_dir = re.split('gov', assembly_ftp)[1]
                #         ftp.login()
                #         ftp.cwd(ftp_dir)
                #         for x in ftp.nlst():
                #             if re.search('[^from]_genomic.fna.gz',x):
                #                 refseq = x
                #             elif x.endswith('_genomic.gbff.gz'):
                #                 gbk = x
                #             else:
                #                 next
                #         ##下载genebank文件
                #         gbk_flag = 0
                #         data = BytesIO()
                #         try:
                #             ftp.retrbinary('RETR ' + gbk, data.write, 1024)
                #             data.seek(0)
                #         except (EOFError,BaseException):
                #             print("Something went wrong in genebank file downloading")
                #             next
                #         try:
                #             uncompressed = gzip.decompress(data.read())
                #             with open(work_dir + '/Candidate_genomes/' + str(tax_id) + '/' + assembly_name + '.gbk', 'wb') as file:
                #                 file.write(uncompressed)
                #             file.close()
                #             gbk_flag = 1
                #         except (zlib.error,OSError,BaseException):
                #             print("Something went wrong in genebank file downloading")
                #             next

                #         ##下载fasta文件
                #         fasta_flag = 0
                #         data = BytesIO()
                #         try:
                #             ftp.retrbinary('RETR ' + refseq, data.write, 1024)
                #             data.seek(0)
                #         except (EOFError,BaseException):
                #             print("Something went wrong in fasta file downloading")
                #             next
                #         try:
                #             uncompressed = gzip.decompress(data.read())
                #             with open(work_dir + '/Candidate_genomes/' + str(tax_id) + '/' + assembly_name + '.fasta', 'wb') as file:
                #                 file.write(uncompressed)
                #             file.close()
                #             fasta_flag = 1
                #         except (zlib.error,OSError,BaseException):
                #             print("Something went wrong in fasta file downloading")
                #             next
                #         # 确保fasta文件和genebank同时下载了
                #         if (fasta_flag + gbk_flag < 2)  and  (fasta_flag + gbk_flag > 0):
                #             try:
                #                 os.system('rm ' + work_dir + '/Candidate_genomes/' + str(tax_id) + '/' + assembly_name + '.*')
                #             except BaseException:
                #                 next
                #         else:
                #             pass
                #     except (IndexError, BaseException):
                #         print("Cant find ftp link for " + str(genome_id))
                #         next

                # except (ValidationError, BaseException): 
                #     print('Something went wrong in Entrez parsing') 
                #     next

            print("\n")
            print("物种/菌株" + str(tax_id) + '下载完成')
            print("\n")
    
    return [species_id_list, organism_list, protein_family_info, all_assembly_id_list]    

def downloadProtein(species_id_list, dataframe, work_dir):
    for tax_id in list(set(species_id_list)):
        if os.path.isdir(work_dir + '/Candidate_genomes/' + str(tax_id) + '/'):
            Organs = dataframe[dataframe['Organism ID']==tax_id]['Organism'].tolist()
            Organs = list(set(Organs))
            Proteins = dataframe[dataframe['Organism ID']==tax_id].index
            Proteins = list(set(Proteins))
            Protein_names = dataframe[dataframe['Organism ID']==tax_id]['Protein names'].tolist()
            Protein_names = list(set(Protein_names))
            content_list = []
            for proteinID in Proteins:
                try:
                    Sequence = dataframe.at[proteinID, 'Sequence'].to_list()
                    Sequence = Sequence[0]
                except AttributeError:
                    Sequence = dataframe.at[proteinID, 'Sequence']
                content = '>' + proteinID + '\n' + Sequence
                content_list.append(content)
            all_protein_seq = '\n'.join(content_list)

            #os.system('mkdir '+ work_dir + '/Candidate_genomes/' + str(tax_id) + '/')

            if os.path.exists(work_dir + '/Candidate_genomes/' + str(tax_id) + '/' + 'proteins.fasta'):
                os.system('rm '+ work_dir + '/Candidate_genomes/' + str(tax_id) + '/' + 'proteins.fasta')
            else:
                pass
            with open (work_dir + '/Candidate_genomes/' + str(tax_id) + '/' + 'proteins.fasta', 'a+') as p_seq:
                p_seq.write(all_protein_seq)
        else:
            next

def sequenceAlignment(sequence_dir):
    dir_list = []
    for root, subdir, files in os.walk(sequence_dir):
        for file in files:
            if file.endswith('gbff'):
                #dir_name = basename(root)
                dir_list.append(root)
    dir_list = list(set(dir_list))

    for dir in dir_list:
        for root, subdir, files in os.walk(dir):
            for file in files:
                if (re.search('[^proteins].fna', file)) and (file.endswith('.fna')):
                    file_path = os.path.join(root,file)
                    query =  os.path.join(root, 'proteins.fasta')
                    os.system("cd "+ dir)
                    os.system("formatdb "+ "-i " + file_path+ ' -p F')
                    os.system("tblastn "+ "-query "  + query + " -outfmt '6 std qlen slen' -db " + file_path + ' -out '+ file_path + '.result.xls')
                else:
                    next


def filterAlignment(df, work_dir, query_string, align_threshold, alignlen_threshold):
    candidate_genomes = []
    protein_family_info = {}
    ## 获得蛋白质所属蛋白家族信息
    for i in df.index.values:
        proteinID = i
        protein_family = df.loc[i, 'ProteinFamily']
        protein_family_info[proteinID] = protein_family
    
    ## 检查过滤后文件是否存在，存在则删除
    os.system("mkdir " + work_dir + "Genome_Alignment/")
    if os.path.isfile(work_dir +  '/Genome_Alignment/Filtered_PF_alignemnt.xls'):
        subprocess.run(["rm", work_dir + "Genome_Alignment/Filtered_PF_alignemnt.xls"], shell=True)
    else:
        pass
    
    ## 开始过滤
    ## 创建存放每个过滤后blalst dataframe的list
    df_list = []
    for root, subdir, files in os.walk(work_dir + 'Candidate_genomes/'):
        for file in files:
            if file.endswith('fna.result.xls'):
                genome_name = re.split('\.fna\.result\.xls', file)[0]
                genome_accession = re.split('_', genome_name)[0:2]
                accessionIDnum = re.split('_', genome_name)[1]
                genome_accession = '_'.join(genome_accession)
                file_path= os.path.join(root,file)
                print(file_path)
                try:
                    blast_result = pd.read_csv(file_path, sep='\t')
                    blast_result.columns = ['proteinID', 'genomeID', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore','qlen', 'slen']
                    query_list = blast_result['proteinID'].tolist()
                    query_list = list(set(query_list))
                    blast_result = blast_result[blast_result['pident'] > align_threshold]
                    blast_result = blast_result[blast_result['length'] > (blast_result['qlen']* alignlen_threshold)]        
                    blast_result = blast_result.dropna(axis=0, how='any')
                    blast_result = blast_result.reset_index(drop=True)
                    blast_result['AccessionID'] = genome_accession
                    blast_result['accessionIDnum'] = accessionIDnum
                    blast_result['ProteinFamily'] = nan
                    print(blast_result)
                    for proteinid in query_list:
                        try:
                            protein_family = protein_family_info[proteinid].to_list()
                            protein_family = protein_family[0]
                        except AttributeError:
                            protein_family = protein_family_info[proteinid]
                        blast_result[blast_result['proteinID'] == proteinid] = blast_result[blast_result['proteinID'] == proteinid].assign(**{'ProteinFamily': protein_family})
                    
                    if blast_result.shape[0] != 0 :
                        pass
                    else:
                        next

                    if ' and ' in query_string:
                        blocks = re.split(' and ', query_string)
                        block_num = len(blocks)
                        block_flag_num = []
                        for block in blocks:
                            block = block.replace('(', '')
                            block = block.replace(')', '')
                            pf_list = re.split(' OR ', block)
                            blast_result_tmp = blast_result[blast_result['ProteinFamily'].isin(pf_list)]
                            if blast_result_tmp.shape[0] > 0:
                                block = 1
                            else:
                                block = 0
                            block_flag_num.append(block)
                        if sum(block_flag_num) < block_num:
                            next
                        else:
                            df_list.append(blast_result)
                            candidate_genomes.append(genome_name)

                    elif ' or ' in query_string:
                        if blast_result.shape[0] != 0 :
                            df_list.append(blast_result)
                            candidate_genomes.append(genome_name)
                        else:
                            next

                    else:
                        if blast_result.shape[0] != 0 :
                            df_list.append(blast_result)
                            candidate_genomes.append(genome_name)
                        else:
                            next
                    
                except EmptyDataError:
                    next
            else:
                next
    try:
        df_all = pd.concat(df_list)
        print(df_all)
        df_all = df_all.drop_duplicates()
        df_all = df_all.reset_index(drop=True)
        ## 增加正负链,蛋白质来源信息库信息以及序列信息
        db_type_list = []
        chain_info_list = []
        for row in range(0,df_all.shape[0],1):
            proteinID = df_all.loc[row,'proteinID']
            proteinID = re.split('\|', proteinID)[0]
            db_type = protein_family_info[proteinID]
            db_type_list.append(db_type)
            sstart = df_all.iloc[row]['sstart']
            send =  df_all.iloc[row]['send']
            if sstart > send:
                tmp = sstart
                df_all.at[row,'sstart'] = send
                df_all.at[row,'send'] = tmp
                chain_info_list.append('-')
            else:
                chain_info_list.append('+')
        df_all['Chain'] = chain_info_list
        df_all['ProteinFamily'] = db_type_list
        df_all.to_csv(work_dir + '/Genome_Alignment/Filtered_PF_alignemnt.xls', sep='\t', index=False)
    except ValueError:
        print('No result left after filtering, try lowering your thresholds !!')

def buildHMMprofile(filtered_blast_result, work_dir, query_string):
    ## 将所有蛋白整合到一起
    complete_fasta = work_dir + 'All_proteins.fasta'

    if os.path.exists(complete_fasta):
        os.system('rm '+ complete_fasta)
    else:
        pass

    for root, subdir, files in os.walk(work_dir + 'Candidate_genomes/'):
        for file in files:
            if file == 'proteins.fasta':
                path = os.path.join(root, file)
                for record in SeqIO.parse(path, "fasta"):
                    id = str(record.id)
                    seq= str(record.seq)
                    content = '>' + id + '\n' + seq + '\n'
                    with open(complete_fasta, 'a+') as out:
                        out.write(content)
                    out.close()
            else:
                next

    if ' and ' in query_string:
        blocks = re.split(' and ', query_string)
        n = 1
        for block in blocks:
            block = block.replace('(', '')
            block = block.replace(')', '')
            pf_list = re.split(' OR ', block)
            tmp_filtered_blast_result = filtered_blast_result[filtered_blast_result['ProteinFamily'].isin(pf_list)]
            block_protein_list = tmp_filtered_blast_result['proteinID'].to_list()

            ## 从下载的蛋白序列中抓取block中包含的蛋白序列

            tmp_block_fasta_file_name = work_dir + 'ProteinFamily' + str(n) + '.fasta'
            n = n +1

            content_list = []
            for record in SeqIO.parse(work_dir + 'All_proteins.fasta', "fasta"):
                id = str(record.id)
                if id in block_protein_list:   
                    content = '>' + id + '\n' + str(record.seq)
                    content_list.append(content)
                    with open(tmp_block_fasta_file_name, 'w+') as out:
                        out.write('\n'.join(content_list))
                    out.close()
                else:
                    next
            

    elif ' or ' in query_string:
        blocks = re.split(' or ', query_string)
        n = 1
        for block in blocks:
            block = block.replace('(', '')
            block = block.replace(')', '')
            pf_list = re.split(' OR ', block)
            tmp_filtered_blast_result = filtered_blast_result[filtered_blast_result['ProteinFamily'].isin(pf_list)]
            block_protein_list = tmp_filtered_blast_result['proteinID'].to_list()

            ## 从下载的蛋白序列中抓取block中包含的蛋白序列

            tmp_block_fasta_file_name = work_dir + 'ProteinFamily' + str(n) + '.fasta'
            n = n +1
            content_list = []
            for record in SeqIO.parse(work_dir + 'All_proteins.fasta', "fasta"):
                id = str(record.id)
                if id in block_protein_list:   
                    content = '>' + id + '\n' + str(record.seq)
                    content_list.append(content)
                    with open(tmp_block_fasta_file_name, 'w+') as out:
                        out.write('\n'.join(content_list))
                    out.close()
                else:
                    next
            
    else:
        n = 1
        query_string = query_string.replace('(', '')
        query_string = query_string.replace(')', '')
        pf_list = re.split(' OR ', query_string)
        tmp_filtered_blast_result = filtered_blast_result[filtered_blast_result['ProteinFamily'].isin(pf_list)]
        block_protein_list = tmp_filtered_blast_result['proteinID'].to_list()

        ## 从下载的蛋白序列中抓取block中包含的蛋白序列

        tmp_block_fasta_file_name = work_dir + 'ProteinFamily' + str(n) + '.fasta'
            
        content_list = []
        for record in SeqIO.parse(work_dir + 'All_proteins.fasta', "fasta"):
            id = str(record.id)
            if id in block_protein_list:   
                content = '>' + id + '\n' + str(record.seq)
                content_list.append(content)
                with open(tmp_block_fasta_file_name, 'w+') as out:
                    out.write('\n'.join(content_list))
                out.close()
            else:
                next

    ## Multiple Seq Alignment and Build HMM profile
    for block_fasta in glob.glob(work_dir + 'ProteinFamily*.fasta'):
        basename = os.path.basename(block_fasta)
        dirname = os.path.dirname(block_fasta)
        basename = re.split('.fasta',basename)[0]
        os.system('clustalo -i ' + block_fasta + ' ' + work_dir + basename + '.st')
        os.system('hmmbuild  --amino ' + os.path.join(dirname, basename + '.hmm ') + ' ' + work_dir + basename + '.st')
        
def runAntismash(query_string, work_dir, Filtered_PF, space_len,neighbour ):
    ## Prepare input files
    block_hmm_list = []
    hmm_file_name_list =[]
    hmm_file_list = []
    for block_hmm in glob.glob(os.path.join(work_dir, '*.hmm')):
        print(block_hmm)
        ## Copy hmm file to antismash folder
        hmm_file_name = os.path.basename(block_hmm)
        hmm_file_list.append(hmm_file_name)
        hmm_file_name = re.split('.hmm',hmm_file_name)[0]
        #os.system('cp ' + block_hmm + ' ' + hmm_folder + '/data/')
        block_hmm_list.append(block_hmm)
        hmm_file_name_list.append(hmm_file_name)
                                ## Define rules in txt files
                                    # RULE loose_test
                                    #     CATEGORY user_define
                                    #     COMMENT test for pilin and sortaseC
                                    #     CUTOFF 20
                                    #     NEIGHBOURHOOD 8
                                    #     CONDITIONS sortaseC and pilin
    if 'and' in query_string:
        # Loose rules file
        content = "RULE" + " " + "loose_rules" + "\n\t"
        content = content + "CATEGORY" +  " " + "user_defined" + "\n\t"
        content = content + "COMMENT" + " " + "user-defined BGC region" + "\n\t"
        content = content + "CUTOFF" + " "  + str(space_len) + "\n\t"
        content = content + "NEIGHBOURHOOD" + " "  + str(neighbour) + "\n\t"
        content = content + "CONDITIONS" +  " "  + " and ".join(hmm_file_name_list)
        with open (work_dir + '/loose.txt', 'w+') as out:
            out.write(content)

        # Relaxed rules file
        content = "RULE" + " " + "relaxed_rules" + "\n\t"
        content = content + "CATEGORY" +  " " + "user_defined" + "\n\t"
        content = content + "COMMENT" + " " + "user-defined BGC region" + "\n\t"
        content = content + "CUTOFF" + " "  + str(space_len) + "\n\t"
        content = content + "NEIGHBOURHOOD" + " "  + str(neighbour) + "\n\t"
        content = content + "CONDITIONS" +  " "  + " and ".join(hmm_file_name_list)
        with open (work_dir + '/relaxed.txt', 'w+') as out:
            out.write(content)

        # Stricted rules file
        content = "RULE" + " " + "strict_rules" + "\n\t"
        content = content + "CATEGORY" +  " " + "user_defined" + "\n\t"
        content = content + "COMMENT" + " " + "user-defined BGC region" + "\n\t"
        content = content + "CUTOFF" + " "  + str(space_len) + "\n\t"
        content = content + "NEIGHBOURHOOD" + " "  + str(neighbour) + "\n\t"
        content = content + "CONDITIONS" +  " "  + " and ".join(hmm_file_name_list)
        with open (work_dir + '/strict.txt', 'w+') as out:
            out.write(content)

    else:
        # Loose rules file
        content = "RULE" + " " + "loose_rules" + "\n\t"
        content = content + "CATEGORY" +  " " + "user_defined" + "\n\t"
        content = content + "COMMENT" + " " + "user-defined BGC region" + "\n\t"
        content = content + "CUTOFF" + " "  + str(space_len) + "\n\t"
        content = content + "NEIGHBOURHOOD" + " "  + str(neighbour) + "\n\t"
        content = content + "CONDITIONS" +  " "  + " or ".join(hmm_file_name_list)
        with open (work_dir + '/loose.txt', 'w+') as out:
            out.write(content)

        # Relaxed rules file
        content = "RULE" + " " + "relaxed_rules" + "\n\t"
        content = content + "CATEGORY" +  " " + "user_defined" + "\n\t"
        content = content + "COMMENT" + " " + "user-defined BGC region" + "\n\t"
        content = content + "CUTOFF" + " "  + str(space_len) + "\n\t"
        content = content + "NEIGHBOURHOOD" + " "  + str(neighbour) + "\n\t"
        content = content + "CONDITIONS" +  " "  + " or ".join(hmm_file_name_list)
        with open (work_dir + '/relaxed.txt', 'w+') as out:
            out.write(content)

        # Stricted rules file
        content = "RULE" + " " + "strict_rules" + "\n\t"
        content = content + "CATEGORY" +  " " + "user_defined" + "\n\t"
        content = content + "COMMENT" + " " + "user-defined BGC region" + "\n\t"
        content = content + "CUTOFF" + " "  + str(space_len) + "\n\t"
        content = content + "NEIGHBOURHOOD" + " "  + str(neighbour) + "\n\t"
        content = content + "CONDITIONS" +  " "  + " or ".join(hmm_file_name_list)
        with open (work_dir + '/strict.txt', 'w+') as out:
            out.write(content)
    

    ## Get path of antismash
    antismash = subprocess.run(['which', 'antismash'], stdout=subprocess.PIPE)
    print(antismash.returncode)
    print('\n')
    if not antismash.returncode:
        antismash= antismash.stdout.decode('utf-8')
        antismash = antismash.replace('\n', '')
        antismash = os.path.dirname(antismash)
        path = subprocess.run(['find', '/opt/', '-name', 'hmm_detection'], stdout=subprocess.PIPE)
        print('hello')
        print('\n')
        print(path)
        print('\n')
        print(path.returncode)
        print('\n')
        if not path.returncode:
            print(path.stdout)
            path = path.stdout.decode('utf-8')
            print(path)
            path = path.replace('\n', '')
            print(path)
        else:
            print("Error found when looking for hmm_detection folder in your docker image, please make sure it's installed correctedly!!")
    else:
        print("Error found when looking for antismash in your docker image, please make sure it's installed correctedly!!")
    if os.path.isdir(os.path.join(os.path.dirname(path), 'hmm_detection_backup/')):
        pass
    else:
        os.system("cp -r {0} {1}".format(path, os.path.join(os.path.dirname(path), 'hmm_detection_backup/')))

    ## Write hmmdetails file
    content_list = []
    for proteinfamily in hmm_file_name_list:
        content = proteinfamily + "\t" + proteinfamily + "\t" + str(space_len) + "\t" + proteinfamily + '.hmm'
        content_list.append(content)
    all_content = '\n'.join(content_list)
    with open(os.path.join(path, 'data','hmmdetails.txt'), 'a+') as HMMDETAIL:
        HMMDETAIL.write(all_content)
    
    os.system('cp ' + os.path.join(path, 'data','hmmdetails.txt') + " " + work_dir)

    ## Write category json file
    category_dict = {
        "user_defined": {
        "description": "user-defined BGC region",
        "version": 1}
    }

    with open(work_dir + "/categories.json","w+") as f:
        json.dump(category_dict,f,indent = 4)

    ## Copy files

    os.system('cp ' + work_dir + "/categories.json " + os.path.join(path, 'data'))
    os.system('cp ' + work_dir +  "strict.txt " + os.path.join(path, 'cluster_rules'))
    os.system('cp ' + work_dir +  "relaxed.txt " + os.path.join(path, 'cluster_rules'))
    os.system('cp ' + work_dir +  "loose.txt " + os.path.join(path, 'cluster_rules'))
    os.system('cp ' + work_dir + "/*hmm " + os.path.join(path, 'data'))
    print('Check param files in antismash')
    os.system('ls ' + os.path.join(path, 'data','*hmm'))
    os.system('ls ' + os.path.join(path, 'cluster_rules','*'))
    print("yes")
    ## run antismash
    candidate_gbk_file_list = list(set(Filtered_PF['accessionIDnum'].to_list()))
    print(candidate_gbk_file_list)
    for accessionIDnum in candidate_gbk_file_list:
        refseq_accession = 'GCF_'+str(accessionIDnum)
        #for gbkfile in glob.glob(os.path.join(work_dir, 'Candidate_genomes', '*','*'+genebak_accession+'*.gbff')):
        for gbkfile in glob.glob(os.path.join(work_dir, 'Candidate_genomes', '*','*'+str(accessionIDnum)+'*.gbff')):
            print(gbkfile)
            subprocess.run(['mkdir','-p',os.path.join(work_dir, 'Anitismash_Result', refseq_accession)])
            os.system('ls -d ' + os.path.join(work_dir, 'Anitismash_Result', refseq_accession))
            antismash_flag = subprocess.run(['antismash', '--hmmdetection-strictness','strict' ,'--taxon', 'bacteria', gbkfile ,'--output-dir', os.path.join(work_dir, 'Anitismash_Result', refseq_accession)])
            if not antismash_flag.returncode:
                pass
            else:
                print('Antismash failed on ' + gbkfile)

def summaryAntismashResult(Result_folder, work_dir):
    all_strain_dict = {}
    strainID2ncID = {}
    assemblyID2ncID = {}
    assemblyIDnumer = {}
    for root, subdir, files in os.walk(Result_folder):
        for file in files:
            if re.search('region', file):
                if file.endswith('gbk'):
                    region_id = os.path.basename(file)
                    strain_id = os.path.basename(root)
                    region_id = re.split('.gbk',region_id )[0]
                    region_dict = {}
                    file_path = os.path.join(root,file)
                    for seq_record in SeqIO.parse(file_path, "genbank"):
                        organism = seq_record.description
                        organism = re.split(',', organism)[0]
                        organism = organism.replace(' plasmid', '')
                        organism = organism.replace(' chromosome', '')
                        ## 获得taxonomy信息
                        try:
                            if len(seq_record.annotations['taxonomy']) < 6:
                                genus = seq_record.annotations['taxonomy'][4] ## Genus level，某些物种缺失某些层级信息
                            else:
                                genus = seq_record.annotations['taxonomy'][5] ## Genus level
                        except IndexError:
                            genus = seq_record.annotations['taxonomy'][-1]
                            print(file_path)
                            print('文件物种层级存在一些问题，请检查')
                            print(seq_record.annotations['taxonomy'])
                        species = ' '.join(re.split(' ', organism)[0:2])
                        ID = seq_record.id

                        for block_hmm in glob.glob(os.path.join(work_dir, '*.hmm')):
                            block_name_num = 0
                            block_name = os.path.basename(block_hmm)
                            block_name = re.split('.hmm', block_name)[0]
                            for feature in seq_record.features:
                                if feature.type == "CDS":
                                    try:
                                        if re.search(block_name, feature.qualifiers['gene_functions'][0]):
                                            block_name_num = block_name_num + 1
                                        else:
                                            next
                                    except KeyError:
                                        next
                                else:
                                    next
                            region_dict[block_name] = block_name_num
                    region_dict['organism'] = organism
                    region_dict['genus'] = genus
                    region_dict['ID'] = ID
                    region_dict['species'] = species
                    #region_dict['StrainID']= strain_id
                    region_dict['regionID'] = region_id
                    all_strain_dict[region_id] = region_dict
                else:
                    next
            else:
                if file.endswith('gbk'):
                    AssemblyID_backup = re.split('_',file)
                    AssemblyID_number = AssemblyID_backup[1]
                    AssemblyID_backup = AssemblyID_backup[0:2]
                    AssemblyID_backup = '_'.join(AssemblyID_backup)
                    #AssemblyID_backup = AssemblyID_backup.replace('.gbk', '')
                    file_path = os.path.join(root,file)
                    for seq_record in SeqIO.parse(file_path, "genbank"):
                        ID = seq_record.id
                        for feature in seq_record.features:
                            try:
                                content = ''.join(feature.qualifiers.get("db_xref"))
                                if re.search('taxon', content):
                                    taxon = content
                                    taxonID = re.split(':',taxon )[-1]
                                else:
                                    next
                            except TypeError:
                                next
                        strainID2ncID[ID] = taxonID
                        Assembly_ID = AssemblyID_backup
                        assemblyID2ncID[ID] = Assembly_ID
                        assemblyIDnumer[ID] = AssemblyID_number


                else:
                    next

    ## 获得物种基因组序列ID对应的taxaID
    Stat_df = pd.DataFrame(all_strain_dict).T
    Stat_df['TaxonID'] = nan
    Stat_df['AssemblyID'] = nan
    Stat_df['AssemblyIDnum'] = nan
    Stat_df['Genus_branches'] = nan

    ## 增加taxonID和AssemblyID信息
    Assembly_id_list = []
    for key in strainID2ncID.keys():
        Stat_df[Stat_df['ID'] == key] = Stat_df[Stat_df['ID'] == key].assign(TaxonID = strainID2ncID[key])
        Stat_df[Stat_df['ID'] == key] = Stat_df[Stat_df['ID'] == key].assign(AssemblyID = assemblyID2ncID[key])
        Stat_df[Stat_df['ID'] == key] = Stat_df[Stat_df['ID'] == key].assign(AssemblyIDnum = assemblyIDnumer[key])
        Assembly_id_list.append(assemblyID2ncID[key])
    Assembly_id_list = list(set(Assembly_id_list))
    # with open ('strain_list.txt', 'w') as strain:
    #     content =  '\n'.join(Assembly_id_list)
    #     strain.write(content)
    #     strain.close

    ## 获取构建进化树的AssemblyID
    Tree_genome_list = []
    Tree_genome_dict = {}
    tree_content = []
    for Genus in list(set(Stat_df['genus'].to_list())):
        tmp_df = Stat_df[Stat_df['genus'] == Genus]
        species_num  = len(set(tmp_df['species'].to_list()))
        species_num = int(species_num)
        Stat_df[Stat_df['genus'] == Genus] = Stat_df[Stat_df['genus'] == Genus].assign(**{"Genus_branches" : species_num})
        ref_assembly_id = tmp_df['AssemblyID'].to_list()[0]
        Tree_genome_list.append(ref_assembly_id)
        Tree_genome_dict[ref_assembly_id]  = Genus

    Stat_df.to_csv(os.path.join(work_dir,'BGC_stat.xls'), sep = '\t')

    for key in Tree_genome_dict.keys():
        line = key + '\t' + Tree_genome_dict[key]
        tree_content.append(line)

    with open (os.path.join(work_dir, 'tree_genomeIDs.txt'), 'w') as genomes:
        content =  '\n'.join(tree_content)
        genomes.write(content)
        genomes.close


def StrainClassification(BGC_stat, db_file):
    df = BGC_stat
    db = pd.read_csv(db_file,sep='\t',encoding = "ISO-8859-1",dtype=str)

    df = df.drop_duplicates()
    df = df.reset_index()
    df['Type'] = nan

    for assemblyid in list(set(df['AssemblyID'].to_list())):
        strain_name = df[df['AssemblyID'] == assemblyid]['organism'].to_list()[0]
        strain_name = re.sub('genome*assembly','',strain_name)
        strain_name = re.sub('complete*genome','',strain_name)
        strain_name = re.sub('DNA','',strain_name)
        strain_name = re.sub('draft*genome','',strain_name)
        strain_name = re.sub('genome','',strain_name)
        strain_name = re.sub('assembly','',strain_name)
        strain_name = strain_name.rstrip()
        strain_name = re.sub(' +', ' ', strain_name)
        strain_name = re.sub('\(','.',strain_name)
        strain_name = re.sub('\)','.',strain_name)
        strain_name = re.sub(r'\s+','.',strain_name)
        taxid = df[df['AssemblyID'] == assemblyid]['TaxonID'].to_list()[0]

        ## 根据assemblyID注释
        if assemblyid in db['assembly_accession'].to_list():
            df[df['AssemblyID'] == assemblyid] = df[df['AssemblyID'] == assemblyid].assign(**{'Type':db[db['assembly_accession'] == assemblyid]['Type'].to_list()[0]})
        elif db[db['genome_name'].str.contains(strain_name, regex=True,case=False) == True].shape[0] > 0:
            ## 根据菌株名字注释
            type = db[db['genome_name'].str.contains(strain_name, regex=True,case=False) == True]['Type'].to_list()[0]
            df[df['AssemblyID'] == assemblyid] = df[df['AssemblyID'] == assemblyid].assign(**{'Type':type})
        else:
            df[df['AssemblyID'] == assemblyid] = df[df['AssemblyID'] == assemblyid].assign(**{'Type':'Others'})

    df['Pathogen'] = ''
    df['Industrial workhorse'] = ''
    df['others'] = ''
    for row in range(0,df.shape[0]):
        if df.loc[row,'Type'] == 'Pathogen':
            df.loc[row,'Pathogen'] = 'yes'
            df.loc[row,'Industrial workhorse'] = 'no'
            df.loc[row,'others'] = 'no'
        elif df.loc[row,'Type'] == 'Industrial workhorse':
            df.loc[row,'Pathogen'] = 'no'
            df.loc[row,'Industrial workhorse'] = 'yes'
            df.loc[row,'others'] = 'no'
        else:
            df.loc[row,'others'] = 'yes'
            df.loc[row,'Pathogen'] = 'no'
            df.loc[row,'Industrial workhorse'] = 'no'

    df.to_csv(os.path.join(work_dir, 'Strain_classification.xls'), sep = '\t')

    Block_names = []

    for col in list(df.columns):
        if 'ProteinFamily' in col:
            Block_names.append(col)
        else:
            next

    for row in range(0,df.shape[0]):
        Patter_num_list = []
        for tmp_name in Block_names:
            tmp_num = df.loc[row, tmp_name]
            Patter_num_list.append(str(tmp_num))
        pattern = ','.join(Patter_num_list)
        df.loc[row,'BGCPattern'] = pattern

    df.to_csv(os.path.join(work_dir, 'Strain_classification.xls'), sep = '\t')


    d = {}
    count = 1
    d['BGCPattern_Name'] = []
    d['Pattern'] = []
    d['BGCPattern number'] = []

    for temp_str in list(set(df['BGCPattern'])):
        d['BGCPattern_Name'].append('BGCPattern'+str(count))
        d['Pattern'].append(temp_str)
        d['BGCPattern number'].append(df[df['BGCPattern'] == temp_str].shape[0])
        count = count + 1
    pattern_df = pd.DataFrame(d)
    pattern_df.to_csv(work_dir + 'Pattern_stat.xls', sep= '\t')

    ## 绘制pielplot
    labels = list(set(df['Type'].to_list()))
    sizes = [(len(df[df['Type'] == x]['AssemblyID'].drop_duplicates().to_list())/(len(df['AssemblyID'].drop_duplicates().to_list()))) * 100 for x in labels]
    #fig, ax = plt.subplots()
    #ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
    #patches, texts  = plt.pie(sizes, startangle=90)
    colours = {
    'Pathogen':'C0',
    'Industrial workhorse':'C1',
    'Others':'C2'
    }
    patches, texts, _  = plt.pie(sizes, autopct='%1.1f%%', colors=[colours[key] for key in labels], startangle=90)
    texts[0].set_fontsize(4)
    plt.legend(patches, labels, loc="best")
    plt.axis('equal')
    plt.savefig(work_dir + 'classification.png', format = 'png', dpi = 250)

    # Construct tree metadata
    for row in range(0,df.shape[0]):
        pattern = df.loc[row,'BGCPattern']
        for p in pattern_df['BGCPattern_Name']:
            tmp_pattern = pattern_df[pattern_df['BGCPattern_Name'] == p]['Pattern'].to_list()[0]
            if pattern == tmp_pattern:
                df.loc[row,'BGCPatternID'] = p
            else:
                next
    print(df)

    for i in pattern_df['BGCPattern_Name']:
        df[i] = nan

    for assemblyid in list(set(df['AssemblyID'].to_list())):
        tmp_df = df[df['AssemblyID'] == assemblyid]
        tmp_pattern_lsit = tmp_df['BGCPatternID'].to_list()
        for i in pattern_df['BGCPattern_Name']:
            count = tmp_pattern_lsit.count(i)
            df[df['AssemblyID'] == assemblyid] = df[df['AssemblyID'] == assemblyid].assign(**{i:count})

    df.to_csv(os.path.join(work_dir, 'Tree_metadata.xls'), sep = '\t')

def buildPhylogenetic_tree(work_dir, input_genomes_folder_name):

    classification_file = os.path.join(work_dir, 'Strain_classification.xls')
    classification_df = pd.read_csv(classification_file, index_col=0, sep = '\t')
    annotated_df = classification_df[classification_df['others'] != 'yes']
    tree_strain_list = list(set(annotated_df['AssemblyID'].to_list()))
    if os.path.isdir(os.path.join(work_dir,'StrainPhylogeneticTree')):
        pass
    else:
        subprocess.run(['mkdir',os.path.join(work_dir,'StrainPhylogeneticTree')])
    
    ## Look for input genomes and copy them into folder
    input_genomes_folder = subprocess.run(['find', '/opt/', '-name', input_genomes_folder_name], stdout=subprocess.PIPE)
    if not input_genomes_folder.returncode:
        input_genomes_folder = input_genomes_folder.stdout.decode('utf-8')
        input_genomes_folder = input_genomes_folder.replace('\n', '')
    else:
        print('Can\'t find input genomes folder in' + input_genomes_folder_name)

    for input_tree_file in glob.glob(os.path.join(input_genomes_folder, '*[fasta,fna]')):
        subprocess.run(['cp', input_tree_file, os.path.join(work_dir,'StrainPhylogeneticTree')])

    for strain_id in tree_strain_list:
        for fasta_file in glob.glob(os.path.join(work_dir, 'Candidate_genomes') + '*/' + strain_id + '*fna'):
            subprocess.run(['cp', fasta_file, os.path.join(work_dir,'StrainPhylogeneticTree')])
    
    ## Create strain tree file
    subprocess.run(['JolyTree.sh', '-i', os.path.join(work_dir,'StrainPhylogeneticTree'), '-b', 'Strain_tree', '-t', '20'])


########### functions ###########

########### Main Run ###########
## Read params
Para_file = args.jsonfile
work_dir = args.workdir
print(Para_file)
params = readPara(Para_file)
align_threshold = params['alignment_threshold']
alignlen_threshold = params['alignment_len_threshold']
Antismash_gap_len = params['Antismash_gap_len']
Antismash_extenson_len = params['Antismash_extenson_len']
Input_genomes = params['Input_genomes']
Email = params['Email']
Entrez.email = Email
query_string = params['Query']
id_list = extract_query_info(query_string)

## Create working dir

if os.path.isfile(work_dir):
    subprocess.run(["mkdir", '-p', work_dir], shell=True)
else:
    pass

## Start Uniprot querying
print("\n")
print('正在根据输入ID检索uniprot数据库')
print("\n")

uniprot_df_list = []
for pf_id in id_list:
   tmp_df = searchUniprot(pf_id)
   if type(tmp_df) == str:
       next
   else:
       uniprot_df_list.append(tmp_df)

combined_df = pd.concat(uniprot_df_list)
combined_df = combined_df.drop_duplicates()
combined_df.to_csv(work_dir + 'Query_result.xls', sep ='\t')
merged_df = findMutualProtein(params['Query'], combined_df)
merged_df = merged_df.drop_duplicates()
merged_df.to_csv(work_dir + 'Filtered_Query_result.xls', sep='\t')

## 下载菌株基因组序列
print("\n")
print('菌株基因组下载开始')
print("\n")

info_list = downloadGenome(merged_df, work_dir)
species_id_list = info_list[0]
organism_list = info_list[1]
protein_family_info = info_list[2]
all_assembly_id_list = info_list[3]

print("\n")
print('菌株基因组下载完成')
print("\n")

## 下载蛋白质序列
print("\n")
print("开始下载物种/菌株过滤后的蛋白质")
print("\n")

downloadProtein(list(set(merged_df['Organism ID'].to_list())), merged_df, work_dir)

print("\n")
print("物种/菌株过滤后的蛋白质下载完成")
print("\n")

## 开始序列比对
print("\n")
print("开始序列比对")
print("\n")

sequenceAlignment(work_dir + '/Candidate_genomes/')

print("\n")
print("序列比对结束")
print("\n")

# 开始比对结果过滤
print("\n")
print("开始序列比对结果过滤")
print("\n")

filterAlignment(merged_df, work_dir, query_string, align_threshold, alignlen_threshold)

print("\n")
print("序列比对过滤结束")
print("\n")

# 开始多序列比对以及构建HMMprofile
print("\n")
print("开始多序列比对以及构建HMMprofile")
print("\n")

Filtered_alignment_file  = work_dir + "Genome_Alignment/Filtered_PF_alignemnt.xls"
Filtered_alignment = pd.read_csv(Filtered_alignment_file, sep = '\t')
buildHMMprofile(Filtered_alignment, work_dir ,query_string )

print("\n")
print("多序列比对以及构建HMMprofile结束")
print("\n")


## 开始运行antismash
print("\n")
print("开始运行antismash")
print("\n")

runAntismash(query_string, work_dir, Filtered_alignment, Antismash_gap_len, Antismash_extenson_len)

print("\n")
print("antismash运行完成")
print("\n")

## 开始总结antismash结果
print("\n")
print("开始总结antismash结果")
print("\n")

Anitismash_Result_path = os.path.join(work_dir, 'Anitismash_Result1')
summaryAntismashResult(Anitismash_Result_path, work_dir)

print("\n")
print("antismash结果总结完成")
print("\n")

## 开始统计菌株分类信息
print("\n")
print("开始统计菌株分类信息")
print("\n")

BGC_stat = pd.read_csv(os.path.join(work_dir, 'BGC_stat.xls'), header= 0, sep = "\t", dtype=str)

db_file = subprocess.run(['find', '/opt/', '-name', 'Database.xls'], stdout=subprocess.PIPE)
if not db_file.returncode:
    db_file = db_file.stdout.decode('utf-8')
    db_file = db_file.replace('\n', '')
else:
    print("Can't find database file, please check")

StrainClassification(BGC_stat, db_file)

print("\n")
print("统计菌株分类信息结束")
print("\n")

## 开始进化树构建
print("\n")
print("开始进化树构建")
print("\n")

buildPhylogenetic_tree(work_dir,Input_genomes)

print("\n")
print("进化树构建结束")
print("\n")

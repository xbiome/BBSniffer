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
import datetime
from subprocess import Popen
from subprocess import TimeoutExpired


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
    BASE = 'http://rest.uniprot.org'
    KB_ENDPOINT = '/uniprotkb/'
    QUERY_ID = pf_id
    UNIPROT_URL = BASE + KB_ENDPOINT + 'search?query=' + QUERY_ID +  '&fields=accession,id,reviewed,protein_name,gene_names,organism_name,organism_id,length,sequence&format=tsv'
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
        mutual_species_id = set(combined_df['Organism (ID)'].to_list())
        for block in queryString_list:
            block = block.replace('(', '')
            block = block.replace(')', '')
            tmp_block_list = []
            for element in re.split(' OR ', block):
                tmp_speciesID = combined_df[combined_df['ProteinFamily'] == element]['Organism (ID)'].to_list()
                tmp_block_list.extend(tmp_speciesID)
            tmp_block_list = list(set(tmp_block_list))
            mutual_species_id = mutual_species_id.intersection(tmp_block_list)
        filtered_df = combined_df[combined_df['Organism (ID)'].isin(list(mutual_species_id)) & ~combined_df['Organism'].isin(['uncultured bacterium','bioreactor metagenome'])]
    else:
        filtered_df = combined_df
    return(filtered_df)

def downloadGenome(dataframe, work_dir):
    species_id_list = []
    organism_list = []
    protein_family_info = {}

    for row in range(0, dataframe.shape[0]):
        species_name = dataframe['Organism'].iloc[row,]
        try:
            genus = re.split(' ', species_name)[0]
            species = re.split(' ', species_name)[1]
            #query_name = ' '.join([genus,species])
            protein_ID = dataframe.iloc[row,].name
            #protein_sequence = merged_df['Sequence'].iloc[row,]
            organism = dataframe['Organism'].iloc[row,]
            organism_id = dataframe['Organism (ID)'].iloc[row,]
            ## Get protein family information for each protein
            proteinfamily_type = dataframe['ProteinFamily'].iloc[row,]

            protein_family_info[protein_ID] = proteinfamily_type
            species_id_list.append(str(organism_id))
            organism_list.append(organism)
        except IndexError:
            next

    # Retrieving genomes from the Entrez databases
    if os.path.exists(work_dir +'/Candidate_genomes/'): 
        os.system("rm -rf " + work_dir + "/Candidate_genomes/*")
        os.system("mkdir " + work_dir + "/Candidate_genomes/")
    else:
        os.system("mkdir" + work_dir + "/Candidate_genomes/")

    content = '\n'.join(list(set(species_id_list)))
    with open(os.path.join(work_dir ,'strains_ID.txt'), 'w+') as out:
        out.write(content)
        

    ## 根据输入文件下载Refseq和genebank
   # subprocess.run(['ncbi-genome-download', '--taxids' , os.path.join(work_dir,'strains_ID.txt'),
   #                 'bacteria' ,'--assembly-levels', 'complete,chromosome',
   #                 '--flat-output','--parallel', '10','-r' , '10', '-o' ,
   #                 work_dir + '/Candidate_genomes/', '-s', 'genbank', '-v',
   #                 '-d'], stdout=subprocess.PIPE)
    subprocess.run(['ncbi-genome-download', '--taxids' , os.path.join(work_dir,'strains_ID.txt'),
                    'bacteria' ,'--assembly-levels', 'complete,chromosome',
                    '-F' ,'fasta,genbank' ,'--flat-output','--parallel', '10', '-r' ,
                    '10','-o' , work_dir + '/Candidate_genomes/', '-v', '-d'], stdout=subprocess.PIPE)

    ## 解压基因组fasta文件
    for gz_file in glob.glob(work_dir + '/Candidate_genomes/*genomic.fna.gz'):
        unzip_status = subprocess.run(['gunzip','-f', gz_file], stdout=subprocess.PIPE)
        if unzip_status.returncode:
            print(gz_file + ' unzip failed!!')
            subprocess.run(["rm", work_dir + "/Candidate_genomes/" + gz_file], shell=True)
            next
        else:
            pass

    ## 检查是否每一个解压后的基因组都有Refseq和对应genebank文件
    for refseq in glob.glob(work_dir + '/Candidate_genomes/*genomic.fna'):
        file_name = os.path.basename(refseq)
        Refseq_ID = re.split('_',file_name)[0:2]
        Refseq_ID = '_'.join(Refseq_ID)
        pattern = work_dir + '/Candidate_genomes/' + Refseq_ID + '*gz'
        file_list = glob.glob(pattern)
        if len(file_list) == 1 :
            pass
        else:
            subprocess.run(['rm', work_dir + '/Candidate_genomes/' + Refseq_ID + '*'])

    return [species_id_list, organism_list, protein_family_info]    

def downloadProtein(dataframe, work_dir):
    dataframe = dataframe.drop_duplicates()
    if os.path.exists(work_dir + '/All_proteins.fasta'):
        os.system('rm'+ work_dir + '/All_proteins.fasta')
    else:
        pass
    id_list = []
    seq_list = []
    content_list = []
    for tmp_id in dataframe.index.values:
        if tmp_id not in id_list:
            id_list.append(tmp_id)
            try:
                tmp_seq = dataframe.loc[tmp_id,'Sequence'].to_list()[0]
                seq_list.append(tmp_seq)
                content ='>' + tmp_id + '\n' + tmp_seq
                content_list.append(content)
            except KeyError:
                tmp_seq = dataframe.loc[tmp_id,'Sequence'].to_list()
                seq_list.append(tmp_seq)
                content ='>' + tmp_id + '\n' + tmp_seq
                content_list.append(content)
            except AttributeError:
                tmp_seq = dataframe.loc[tmp_id,'Sequence']
                seq_list.append(tmp_seq)
                content ='>' + tmp_id + '\n' + tmp_seq
                content_list.append(content)
        else:
            next         
    proteins = '\n'.join(content_list)
    with open (work_dir + '/All_proteins.fasta', 'w+') as p_seq:
        p_seq.write(proteins)

def sortDownloadedData(work_dir):
    ID_dict = {}
    datafolder = os.path.join(work_dir,'Candidate_genomes')
    all_fna_files = glob.glob(os.path.join(datafolder, '*', '*fna'))
    ID_list = []
    for fna_file in all_fna_files:
        refseq_file_name = os.path.basename(fna_file)
        tmp_folder = os.path.dirname(fna_file)
        ID = re.split('_', refseq_file_name)[1]
        ID = re.split('.fna', ID)[0]
        gbff_file_pattern = os.path.join(tmp_folder, '*'+ ID + '*.gbff.gz')
        ## Check if relevant gbff file exists, jump to next id if not.
        for gbff_file in glob.glob(gbff_file_pattern):
            if os.path.isfile(gbff_file):
                flag = 1
            else:
                flag = 0
        if flag :
            ID_list.append(ID)
        else:
            next
            print(fna_file + ' 无对应gbff文件')
    ID_dict['accessionIDnum'] = ID_list
    all_strains_id = pd.DataFrame.from_dict(ID_dict)
    all_strains_id.to_csv(work_dir + '/All_Strains_for_Antismash.xls', sep = '\t')
    return(all_strains_id)



def buildHMMprofile(filtered_query, work_dir, query_string):
    ## 将所有蛋白整合到一起
    complete_fasta = work_dir + '/All_proteins.fasta'

    if ' and ' in query_string:
        blocks = re.split(' and ', query_string)
        n = 1
        for block in blocks:
            block = block.replace('(', '')
            block = block.replace(')', '')
            pf_list = re.split(' OR ', block)
            tmp_filtered_result = filtered_query[filtered_query['ProteinFamily'].isin(pf_list)]
            block_protein_list = tmp_filtered_result['Entry'].to_list()

            ## 从下载的蛋白序列中抓取block中包含的蛋白序列

            tmp_block_fasta_file_name = work_dir + 'ProteinFamily' + str(n) + '.fasta'
            n = n +1

            content_list = []
            for record in SeqIO.parse(complete_fasta, "fasta"):
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
            tmp_filtered_result = filtered_query[filtered_query['ProteinFamily'].isin(pf_list)]
            block_protein_list = tmp_filtered_result['Entry'].to_list()
            print(pf_list)
            print(tmp_filtered_result)
            print(block_protein_list)
            ## 从下载的蛋白序列中抓取block中包含的蛋白序列

            tmp_block_fasta_file_name = work_dir + 'ProteinFamily' + str(n) + '.fasta'
            n = n +1
            content_list = []
            for record in SeqIO.parse(complete_fasta, "fasta"):
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
        tmp_filtered_result = filtered_query[filtered_query['ProteinFamily'].isin(pf_list)]
        block_protein_list = tmp_filtered_result['Entry'].to_list()

        ## 从下载的蛋白序列中抓取block中包含的蛋白序列

        tmp_block_fasta_file_name = work_dir + 'ProteinFamily' + str(n) + '.fasta'
            
        content_list = []
        for record in SeqIO.parse(complete_fasta, "fasta"):
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
        subprocess.run(['clustalo', '-i', block_fasta, '-o', os.path.join(work_dir, basename + '.st'), '--threads','20'])
        os.system('hmmbuild  --amino ' + os.path.join(dirname, basename + '.hmm ') + ' ' + work_dir + basename + '.st')
        
def runAntismash(query_string, work_dir, space_len,neighbour, n_threads):
    ## Prepare input files
    block_hmm_list = []
    hmm_file_name_list =[]
    hmm_file_list = []
    for block_hmm in glob.glob(os.path.join(work_dir, '*.hmm')):
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
    if not antismash.returncode:
        antismash= antismash.stdout.decode('utf-8')
        antismash = antismash.replace('\n', '')
        antismash = os.path.dirname(antismash)
        path = subprocess.run(['find', '/opt/', '-name', 'hmm_detection'], stdout=subprocess.PIPE)
        if not path.returncode:
            path = path.stdout.decode('utf-8')
            path = path.replace('\n', '')
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
    
    ID_list = []
    candidate_gbk_file_list =  glob.glob(os.path.join(work_dir,'Candidate_genomes','*.gbff.gz'))
    for candidate_gbk_file in candidate_gbk_file_list:
        file_name = os.path.basename(candidate_gbk_file)
        ID = re.split('_', file_name)[0:2]
        ID = '_'.join(ID)
        refseq_accession = 'GCF_'+str(ID)
        os.system("mkdir -p " + os.path.join(work_dir, 'Anitismash_Result', refseq_accession))
        ID_list.append(ID)
    
    input = zip(ID_list, candidate_gbk_file_list)
    command_list = []
    for ID, input_file in input:
        refseq_accession = 'GCF_'+str(ID)
        command = ' '.join(['antismash', '--hmmdetection-strictness','strict' , '--genefinding-tool prodigal', '-c', '1', '--taxon', 'bacteria', input_file , '--output-dir', os.path.join(work_dir, 'Anitismash_Result', refseq_accession)])
        command_list.append(command)
    
    for i in range(0, len(command_list), n_threads):
        tmp_command_list = command_list[i:i+n_threads]
        procs = [ Popen(i,  shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE) for i in tmp_command_list]
        for p in procs:
            try:
                outs, errs = p.communicate(timeout = 400)
                if errs:
                    print(errs.decode('ascii'))
            except TimeoutExpired:
                p.kill()
                outs, errs = p.communicate()
                print(errs.decode('ascii'))


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
                                        for function in feature.qualifiers['gene_functions']:
                                            if re.search(block_name, function):
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


def StrainClassification(BGC_stat, db_file, work_dir):
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
        strain_name = re.sub('[-_]','.',strain_name)
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
    classification_df = pd.read_csv(classification_file, index_col=0, sep = '\t', dtype=str)
    #annotated_df = classification_df
    annotated_df = classification_df[classification_df['others'] != 'yes']
    tree_strain_list = list(set(annotated_df['AssemblyIDnum'].to_list()))
    if os.path.isdir(os.path.join(work_dir,'StrainPhylogeneticTree/Input')):
        pass
    else:
        subprocess.run(['mkdir','-p',os.path.join(work_dir,'StrainPhylogeneticTree/Input/')])
    
    ## Look for input genomes and copy them into folder
    input_genomes_folder = subprocess.run(['find', work_dir, '-name', input_genomes_folder_name], stdout=subprocess.PIPE)
    if not input_genomes_folder.returncode:
        input_genomes_folder = input_genomes_folder.stdout.decode('utf-8')
        input_genomes_folder = input_genomes_folder.replace('\n', '')
    else:
        print('Can\'t find input genomes folder in' + input_genomes_folder_name)

    input_tree_list_fna = glob.glob(os.path.join(input_genomes_folder, '*fna'))
    input_tree_list_fasta = glob.glob(os.path.join(input_genomes_folder, '*fasta'))
    input_tree_list = input_tree_list_fna + input_tree_list_fasta
    
    for input_tree_file in input_tree_list:
        print(input_tree_file)
        # 获取Input Genomes的NC ID
        with open(input_tree_file) as handle:
            nc_id_list = []
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in annotated_df['ID'].to_list():
                    next
                else:
                    subprocess.run(['cp', input_tree_file, os.path.join(work_dir,'StrainPhylogeneticTree/Input/')])

    for strain_id in tree_strain_list:
        for fasta_file in glob.glob(os.path.join(work_dir, 'Candidate_genomes') + '/*' + str(strain_id) + '*'):
            if fasta_file.endswith('fasta') or fasta_file.endswith('fna'):
                file_size = os.path.getsize(fasta_file)
                if file_size == 0:
                    next
                else:
                    subprocess.run(['cp', fasta_file, os.path.join(work_dir,'StrainPhylogeneticTree/Input/')])
            else:
                next

    ## Create strain tree file
    subprocess.run(['JolyTree.sh', '-i', os.path.join(work_dir,'StrainPhylogeneticTree/Input/'), '-b', os.path.join(work_dir,'StrainPhylogeneticTree/','Strain_tree'), '-t', '20'])

def Infer_candidate_strains(work_dir):
    annotated_df = pd.read_table(os.path.join(work_dir,'Strain_classification.xls'), sep = '\t',dtype=str)
    # 去掉ID栏多余的前缀
    newID_list = []
    for id in annotated_df['ID'].to_list():
        new_id = re.sub(r'.+_', '',id)
        newID_list.append(new_id)
    annotated_df['ID'] = newID_list

    distance_file = glob.glob(work_dir + '/StrainPhylogeneticTree/' + 'Strain_tree' + '.d')
    distance_file = distance_file[0]
    distance_matrtix = pd.read_csv(distance_file, skiprows = 1, index_col = 0,sep ='\s+',header=None)
    newname = []
    for i in distance_matrtix.index.values:
        newname.append(re.split('_', i)[1])
    name = pd.Series(newname)
    distance_matrtix = distance_matrtix.set_index(name)
    distance_matrtix.columns = newname

    for input_tree_file in glob.glob(os.path.join(work_dir,'Input_genomes/*')):
        if input_tree_file.endswith('fasta') or input_tree_file.endswith('fna'):
        # 获取Input Genomes的NC ID
            with open(input_tree_file) as handle:
                Assemlbly_ID_list = []
                for record in SeqIO.parse(handle, "fasta"):
                    tmp_id = re.split('_',record.id)[1]
                    if tmp_id in annotated_df['ID'].to_list():
                        Assemlbly_ID_num = annotated_df[annotated_df['ID'] == tmp_id]['AssemblyID'].to_list()[0]
                        Assemlbly_ID_num = re.split('_',Assemlbly_ID_num)[1]
                        Assemlbly_ID_list.append(Assemlbly_ID_num)
                    else:
                        next
        else:
            next


    content = '候选菌株信息:\n'

    if len(Assemlbly_ID_list) == 0:
        content = '未能找到候选工程菌'
    else:
        pass

    candidate_dict = {}
    for pathogen in set(Assemlbly_ID_list):
        tmp_pathogen = annotated_df[annotated_df['AssemblyIDnum'] == pathogen]['organism'].to_list()[0]
        #tmp_df = annotated_df
        tmp_df = annotated_df[annotated_df['Type'] == 'Industrial workhorse']    

        content = content + tmp_pathogen + '的候选工程菌为:\n'
        tmp_distance_series = distance_matrtix.loc[pathogen,]
        tmp_distance_series = tmp_distance_series.sort_values(ascending = True)
        tmp_candidate = tmp_distance_series.rename('ID')

        ## 只挑选前10的物种输出到候选菌株信息.txt中
        count = 0
        candidate_list = []
        all_candidate_list = []
        for id in tmp_candidate.index.values:
            if id in tmp_df['AssemblyIDnum'].tolist():
                if id == pathogen:
                    next
                else:
                    if count < 11:
                        candidate = tmp_df[tmp_df['AssemblyIDnum'] == id]['organism'].tolist()[0]
                        candidate_list.append(candidate)
                        count = count + 1
                    else:
                        break
            else:
                next
        try:
            content = content + '\n'.join(candidate_list)
        except KeyError:
            content = '未能找到候选工程菌'
        
        ## 挑选所有物种输出值所有候选菌株信息.xls中
        for id in tmp_candidate.index.values:
            if id in tmp_df['AssemblyIDnum'].tolist():
                candidate = tmp_df[tmp_df['AssemblyIDnum'] == id]['organism'].tolist()[0]
                all_candidate_list.append(candidate)
            else:
                next
        candidate_dict[tmp_pathogen] = all_candidate_list
    ## 将候选菌株写入文件
    candidate_df = pd.DataFrame.from_dict(candidate_dict)
    candidate_df.to_csv(os.path.join(work_dir,'所有候选菌株信息.xls'), sep = '\t')

    with open(os.path.join(work_dir,'候选菌株信息.txt'), 'w') as out:
        out.write(content)

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
n_threads = params['Antismash_threads']
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
print(datetime.datetime.now())
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

print("\n")
print('检索uniprot数据库完成')
print(datetime.datetime.now())
print("\n")

## 下载菌株基因组序列
print("\n")
print('菌株基因组下载开始')
print(datetime.datetime.now())
print("\n")

info_list = downloadGenome(merged_df, work_dir)
species_id_list = info_list[0]
organism_list = info_list[1]
protein_family_info = info_list[2]

print("\n")
print('菌株基因组下载完成')
print(datetime.datetime.now())
print("\n")

## 下载蛋白质序列
print("\n")
print("开始下载物种/菌株过滤后的蛋白质")
print(datetime.datetime.now())
print("\n")

downloadProtein(merged_df, work_dir)

print("\n")
print("物种/菌株过滤后的蛋白质下载完成")
print(datetime.datetime.now())
print("\n")

# 开始整理下载数据并生成统计文件All_Strains_for_Antismash.xls
print("\n")
print("开始整理下载数据并生成统计文件All_Strains_for_Antismash.xls")
print(datetime.datetime.now())
print("\n")

sortDownloadedData(work_dir)

print("\n")
print("整理下载数据并生成统计文件All_Strains_for_Antismash.xls结束")
print(datetime.datetime.now())
print("\n")

# 开始多序列比对以及构建HMMprofile
print("\n")
print("开始多序列比对以及构建HMMprofile")
print(datetime.datetime.now())
print("\n")

# All_Strains_for_Antismash_file  = work_dir + "/All_Strains_for_Antismash.xls"
# All_Strains_for_Antismash = pd.read_csv(All_Strains_for_Antismash_file, sep = '\t', dtype=str)
Filtered_query_file = work_dir + "/Filtered_Query_result.xls"
Filtered_query = pd.read_csv(Filtered_query_file, sep = '\t', dtype=str)
buildHMMprofile(Filtered_query, work_dir ,query_string)

print("\n")
print("多序列比对以及构建HMMprofile结束")
print(datetime.datetime.now())
print("\n")


## 开始运行antismash
print("\n")
print("开始运行antismash")
print(datetime.datetime.now())
print("\n")

runAntismash(query_string, work_dir, Antismash_gap_len, Antismash_extenson_len, n_threads)

print("\n")
print("antismash运行完成")
print(datetime.datetime.now())
print("\n")

## 开始总结antismash结果
print("\n")
print("开始总结antismash结果")
print(datetime.datetime.now())
print("\n")

Anitismash_Result_path = os.path.join(work_dir, 'Anitismash_Result')
summaryAntismashResult(Anitismash_Result_path, work_dir)

print("\n")
print("antismash结果总结完成")
print(datetime.datetime.now())
print("\n")

## 开始统计菌株分类信息
print("\n")
print("开始统计菌株分类信息")
print(datetime.datetime.now())
print("\n")

BGC_stat = pd.read_csv(os.path.join(work_dir, 'BGC_stat.xls'), header= 0, sep = "\t", dtype=str)

db_file = subprocess.run(['find', '/opt/', '-name', 'Database.xls'], stdout=subprocess.PIPE)
if not db_file.returncode:
    db_file = db_file.stdout.decode('utf-8')
    db_file = db_file.replace('\n', '')
else:
    print("Can't find database file, please check")

StrainClassification(BGC_stat, db_file, work_dir)

print("\n")
print("统计菌株分类信息结束")
print(datetime.datetime.now())
print("\n")

## 开始进化树构建
print("\n")
print("开始进化树构建")
print(datetime.datetime.now())
print("\n")

buildPhylogenetic_tree(work_dir,Input_genomes)

print("\n")
print("进化树构建结束")
print(datetime.datetime.now())
print("\n")

## 开始挑选候选菌株
print("\n")
print("开始挑选候选菌株")
print(datetime.datetime.now())
print("\n")

Infer_candidate_strains(work_dir)

print("\n")
print("挑选候选菌株结束")
print(datetime.datetime.now())
print("\n")

import os
import re
import urllib3
from bs4 import BeautifulSoup, NavigableString, Tag
import argparse


########### Parse arguments ###########
parser = argparse.ArgumentParser(description="BGCSniffer")
parser.add_argument("-input_keywords", "--keywords", dest="keywords", help= "Please input the keywords you want to search")
parser.add_argument("-work_dir", "--workdir", dest="workdir", help= "Please specify Working directory")
args = parser.parse_args()
########### Parse arguments ###########


def KeywordsFetching(Input_string):

    # Input handling
    Input_text = re.split('\s+', Input_string)
    Input_text = '+'.join(Input_text)
    
    # Constructt URL string
    http = urllib3.PoolManager()
    url = 'https://pfam.xfam.org/search/keyword?query='
    url = url + Input_text

    # Fetch html file
    response = http.request('GET', url)
    soup = BeautifulSoup(response.data,'html.parser')
    Family_ID = [tag.find("a")["href"] for tag in soup.select("td:has(a)")]
    Family_ID = [re.split('/family/',ID )[1] for ID in Family_ID]

    ## Remove duplicate ID
    Family_ID = list(set(Family_ID))
    Family_ID.sort()
    return(Family_ID)


Family_IDs = KeywordsFetching(args.keywords)
Protein_ID_file = os.path.join(args.workdir, 'ProteinID.list')

content = '\n'.join(Family_IDs)
with open(Protein_ID_file, 'w') as out:
    out.write(content)


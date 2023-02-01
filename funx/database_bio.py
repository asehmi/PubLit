import streamlit as st
import pandas as pd
from Bio import Entrez, Medline
from numpy import empty

# Change this email to your email address
Entrez.email = "A.N.Other@example.com"
# Load the publications based on keywords and value
@st.experimental_memo(suppress_st_warning=True,show_spinner=False)
def get_pub(val,keyword):
    result = Entrez.read(Entrez.esearch(db = "pubmed", retmax = int(val), term = keyword))
    ids = result["IdList"]
    batch_size = 100
    batches = [ids[x: x + 100] for x in range(0, len(ids), batch_size)]   
    record_list = []
    for batch in batches:
        handle = Entrez.efetch(db = "pubmed", id = batch, rettype = "medline", retmode = "text")
        records = Medline.parse(handle)
        record_list.extend(list(records))
    
    # Datas are loaded as dataframe
    publication_data = pd.DataFrame(record_list)
    # Cleaning the data
    publication_data.dropna(subset=['EDAT'], inplace=True)
    publication_data["Year"] = (publication_data["EDAT"].astype(str).str[0:4].astype(int))
                # Keep the ones to display on table
    dic = {'Published on': publication_data["DP"],  'Title' : publication_data['TI'], 'Authors' : publication_data["AU"],'PubchemID': publication_data["PMID"], "DOI":publication_data["LID"] }
    df = pd.DataFrame(dic)
    return publication_data,df

# get idea of max publication with the keyword
@st.cache(suppress_st_warning=True,show_spinner=False)
def avail_pub(keyword):
    search_handle = Entrez.esearch(db = "pubmed", retmax = 10, term = keyword)
    search_result = Entrez.read(search_handle)
    max_pub = search_result["Count"]

    return max_pub
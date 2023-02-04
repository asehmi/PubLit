import streamlit as st
import pandas as pd
from Bio import Entrez, Medline

# Change this email to your email address
Entrez.email = "A.N.Other@example.com"

@st.experimental_memo(suppress_st_warning=True,show_spinner=False)
def get_pub(val, keyword):
    """
    This function retrieves the publication data from Pubmed database by searching the keyword.
    
    Parameters:
    val (int): The maximum number of results to return.
    keyword (str): The keyword to search for in the Pubmed database.
    
    Returns:
    tuple: Tuple of two dataframes - 'publication_data' and 'df'. 
    'publication_data' is the raw data retrieved from the Pubmed database, 
    and 'df' is the cleaned data with columns 'Published on', 'Title', 'Authors', 'PubchemID', 'DOI'.
    
    """
    # Search for the keyword in the Pubmed database and retrieve the maximum number of results specified by 'val'
    result = Entrez.read(Entrez.esearch(db="pubmed", retmax=int(val), term=keyword))
    ids = result["IdList"]
    
    # Batch process the retrieved results to avoid exceeding the API call limit
    batch_size = 100
    batches = [ids[x:x + 100] for x in range(0, len(ids), batch_size)]   
    record_list = []
    for batch in batches:
        handle = Entrez.efetch(db="pubmed", id=batch, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        record_list.extend(list(records))
    
    # Load the data as a dataframe
    publication_data = pd.DataFrame(record_list)
    # Clean the data by removing rows with missing publication dates and adding a "Year" column
    publication_data.dropna(subset=["EDAT"], inplace=True)
    publication_data["Year"] = (publication_data["EDAT"].astype(str).str[0:4].astype(int))
    
    # Create a new dataframe with columns 'Published on', 'Title', 'Authors', 'PubchemID', 'DOI'
    dic = {'Published on': publication_data["DP"], 
        'Title': publication_data['TI'], 'Authors': publication_data["AU"],
        "DOI": publication_data["LID"], 'PubchemID': publication_data["PMID"]}
    df = pd.DataFrame(dic)
    # Clean the 'DOI' column by removing '[doi]' and adding 'https://doi.org/' in front of each entry
    df['DOI'] = df['DOI'].str.replace(' \[doi\]', '')
    df['DOI'] = 'https://doi.org/' + df['DOI']
    
    return publication_data, df


# get idea of max publication with the keyword
@st.cache(suppress_st_warning=True,show_spinner=False)
def avail_pub(keyword):
    """
    This function takes a keyword as input and returns the number of available publications in pubmed database
    that match the keyword.

    Parameters:
    keyword (str): The keyword to be searched in the pubmed database

    Returns:
    int: The number of available publications in pubmed database that match the keyword
    
    Example:
    >>> avail_pub("cancer research")
    59845
    """
    search_handle = Entrez.esearch(db = "pubmed", retmax = 100, term = keyword)
    search_result = Entrez.read(search_handle)
    max_pub = search_result["Count"]

    return max_pub

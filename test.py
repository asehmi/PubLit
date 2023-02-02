import arxiv
import streamlit as st
import pandas as pd
import re
from st_aggrid import AgGrid
def arxiv_search(query,max_results):
    '''
    Search the arXiv repository 
    '''
    search = arxiv.Search(
    query = query ,
    max_results = max_results,
    sort_by = arxiv.SortCriterion.SubmittedDate,
    sort_order = arxiv.SortOrder.Descending
    )

    published_on = []
    title = []
    authors = []
    arxivid = [] # result.entry_id: A url http://arxiv.org/abs/{id
    doi = []
    abstract = []
    pdf_url = []
    journal_details = []
    for result in search.results():
        published_on.append(result.published)
        title.append(result.title)
        authors.append(result.authors)
        arxivid.append(result.entry_id)
        doi.append(result.doi)
        abstract.append(result.summary)
        pdf_url.append(result.pdf_url)
        journal_details.append(result.journal_ref)
    
    dic = {'Published on': published_on,  
            'Title' : title, 
            'Authors' : authors,
            'arxivID': arxivid, 
            "DOI":doi 
            }
    return pd.DataFrame(dic)

df = arxiv_search("phycobilisomes", 10)
df['Authors'] = df['Authors'].astype(str)
df['Authors'] = df['Authors'].apply(lambda x: re.findall(r"\'(.*?)\'", x))
df['Published on'] = pd.to_datetime(df['Published on'])
df['Published on'] = df['Published on'].dt.strftime('%Y %b %d')
# Add 'https://doi.org/' in front of each row of the 'DOI' column
df['DOI'] = 'https://doi.org/' + df['DOI']
AgGrid(df)


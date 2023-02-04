import streamlit as st
from collections import Counter
from itertools import combinations
from stvis import pv_static
from pyvis import network as net

import networkx as nx
import plotly.express as px
import pandas as pd

# Top Authors
def plot_top_authors(pubdf):      
    n_author = 10
    authors_flat = [author for authors in list(pubdf["FAU"].dropna())for author in authors]
    top10authors = pd.DataFrame.from_records(
                Counter(authors_flat).most_common(n_author), columns=["Name", "Count"]
            )
    fig = px.histogram(
                    data_frame = top10authors,
                    x = "Name", y = "Count",
                    title = "Top {} Contributing authors".format(n_author),
                    color = "Name",                                            
                )   
    return fig 

# Top Journals with the key words
def plot_top_journals(pdf):
    top_journals = pd.DataFrame.from_records(
            Counter(pdf["TA"]).most_common(10),
            columns=["Journal", "Count"],
        )   
    fig = px.pie(data_frame = top_journals,
                names = "Journal", 
                values = "Count", 
                title = "Top Journals",
                )
    return fig

# Associated Keywords with the search
def top_keyworkds(pdf):
    if 'OT' in pdf:
        flat_kw = [
                _.lower()
                for kws in list(pdf["OT"].dropna())
                for kw in kws
                for _ in kw.split(" ")
                ] 
        # Top 10 Keywords are choosen - may be user implementation req        
        top_keys = pd.DataFrame.from_records(
                Counter(flat_kw).most_common(10), columns=["Keyword", "Count"]
                )
        fig = px.histogram(
                        data_frame = top_keys,
                        x = "Count", y = "Keyword",
                        title = "Top Keywords related to the articles",
                        color = "Keyword",
                    )
    # Quick bug fix, Not okay                  
    else:
        fig = "Keywords Plot Unavailable"       
    return fig

# Journal Yearly Trends
def year_journal_trend(pdf):
    yearly = pd.DataFrame(pdf["Year"].value_counts().reset_index())
    yearly.columns = ["Year", "Count"]
    fig = px.scatter(yearly, x="Year", y="Count",
	            size="Count", color="Year",
                hover_name = "Year", log_x=True, size_max=60, title= "Yearly trends")
    return fig

# Creating the network map
def plot_connection(pdf):    
    #st.write(pdf)
    authors = pdf["FAU"].dropna()
    author_connections = list(map(lambda x: list(combinations(x[::-1], 2)), authors))
    flat_connections = [item for sublist in author_connections for item in sublist]
    # Create a dataframe with the connections
    df = pd.DataFrame(flat_connections, columns=["From", "To"])
    df_graph = df.groupby(["From", "To"]).size().reset_index()
    df_graph.columns = ["From", "To", "Count"]
    G = nx.from_pandas_edgelist(df_graph, source="From", target="To", edge_attr= "Count")   
    # Limit to TOP 50 authors (needs fix)
    authors_flat = [author for authors in list(pdf["FAU"].dropna())for author in authors]
    top50authors = pd.DataFrame.from_records(Counter(authors_flat).most_common(50), columns=["Name", "Count"])
    top50_nodes = (n for n in list(G.nodes()) if n in list(top50authors["Name"]))
    G_50 = G.subgraph(top50_nodes)
    for n in G_50.nodes():
        G_50.nodes[n]["publications"] = int(top50authors[top50authors["Name"] == n]["Count"])
    #st.write(G_50.nodes(data=True))
    g = net.Network(
        height='600px', width='900px',
        heading='Authors Network',
        font_color='black',
        # bgcolor='#222222',
        notebook = False
        )
    g.from_nx(G_50,default_node_size = 15)
    pv_static(g)   
    return g,top50authors    


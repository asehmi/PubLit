import streamlit as st
from st_aggrid import AgGrid,GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder

from collections import Counter
from itertools import combinations
from stvis import pv_static
from pyvis import network as net

import networkx as nx
import plotly.express as px
import pandas as pd
from numpy import empty

# Ag-Grid Implementation
def grid_table(df):
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_pagination(enabled=True)
    gb.configure_side_bar()
    gb.configure_selection('multiple', use_checkbox=True, groupSelectsChildren=True, groupSelectsFiltered=True)
    gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, editable=True)
    gridOptions = gb.build()
    grid_response = AgGrid(df, 
                gridOptions = gridOptions, 
                enable_enterprise_modules = True,
                fit_columns_on_grid_load = False,
                width='100%',
                theme = "dark",
                update_mode = GridUpdateMode.SELECTION_CHANGED,
                allow_unsafe_jscode=True)
    df = grid_response['data']
    selected_rows = grid_response["selected_rows"]
    return df,selected_rows

# Top Authors
def plot_top_authors(pubdf):      
    n_author = st.sidebar.slider(label = "Contributing Authors", min_value = 0, max_value = 100, value = 5)
    authors_flat = [author for authors in list(pubdf["FAU"].dropna())for author in authors]
    top10authors = pd.DataFrame.from_records(
                Counter(authors_flat).most_common(n_author), columns=["Name", "Count"]
            )
    fig = px.histogram(
                    data_frame = top10authors,
                    x = "Name", y = "Count",
                    title = "Top {} Contributing authors".format(n_author),
                    color = "Name",
                    color_discrete_sequence = px.colors.qualitative.Pastel2,    
                                           
                )    
    return fig 

# Top Journals with the key words
def plot_top_journals(pdf):
    # By default 10 journals are choosen , Need fix
    top_journals = pd.DataFrame.from_records(
            Counter(pdf["TA"]).most_common(10),
            columns=["Journal", "Count"],
        )   
    fig = px.pie(data_frame = top_journals,
                 names = "Journal", 
                 values = "Count", 
                 title = "Top Journals",
                 color_discrete_sequence = px.colors.qualitative.Pastel2)
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
                        color_discrete_sequence = px.colors.qualitative.Pastel2,                                                
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
                color_discrete_sequence = px.colors.qualitative.Pastel2,
                hover_name = "Year", log_x=True, size_max=60, title= "Yearly trends")
    return fig

# Creating the network map
def plot_connection(pdf):    
    authors = pdf["FAU"].dropna()
    author_connections = list(
        map(lambda x: list(combinations(x[::-1], 2)), authors)
    )
    flat_connections = [item for sublist in author_connections for item in sublist]

    # Create a dataframe with the connections
    df = pd.DataFrame(flat_connections, columns=["From", "To"])
    df_graph = df.groupby(["From", "To"]).size().reset_index()
    df_graph.columns = ["From", "To", "Count"]
    G = nx.from_pandas_edgelist(
        df_graph, source="From", target="To", edge_attr="Count"
    )   
    # Limit to TOP 50 authors (needs fix)
    authors_flat = [author for authors in list(pubdf["FAU"].dropna())for author in authors]
    top50authors = pd.DataFrame.from_records(
        Counter(authors_flat).most_common(50), columns=["Name", "Count"]
    )
    top50_nodes = (n for n in list(G.nodes()) if n in list(top50authors["Name"]))
    G_50 = G.subgraph(top50_nodes)
    for n in G_50.nodes():
        G_50.nodes[n]["publications"] = int(
            top50authors[top50authors["Name"] == n]["Count"]
        )
    g = net.Network(height='600px', width='900px',heading='Authors Network',font_color='white',bgcolor='#222222',notebook = True)
    #bgcolor='#222222'
    g.from_nx(G_50,default_node_size= 15,default_edge_weight=3)
    pv_static(g)
    return g,top50authors    


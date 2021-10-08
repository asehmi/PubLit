# Modules
from Bio import Entrez
from Bio import Medline
from numpy import empty
import streamlit as st
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from stqdm import stqdm
from collections import Counter
import pandas as pd
import plotly.express as px
from st_aggrid import AgGrid,GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder
import json
from itertools import combinations
import networkx as nx
from stvis import pv_static
from pyvis import network as net
from streamlit_lottie import st_lottie
import requests
from scholarly import scholarly
###########################

# Load the publications based on keywords and value
@st.experimental_memo(suppress_st_warning=True,show_spinner=False)
def get_pub(val,keyword):
    result = Entrez.read(Entrez.esearch(db = "pubmed", retmax = int(val), term = keyword))
    ids = result["IdList"]
    batch_size = 100
    batches = [ids[x: x + 100] for x in range(0, len(ids), batch_size)]   
    record_list = []
    for batch in stqdm(batches):
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

# Convert Datframe for saving
@st.experimental_memo
def convert_df(df):
    return df.to_csv().encode('utf-8')

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

def load_lottieurl(url: str):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

# Schloarly Implementation (need to work on!)
def author_information(name):
    search_query = scholarly.search_author(name)
    author = scholarly.fill(next(search_query))
    [pub['bib']['title'] for pub in author['publications']]
    pub = scholarly.fill(author['publications'][0])
    return pub




####################
# MAIN APP
#####################

about_info = ''' # Scientific Literatures and Analytics 
Search using specific keyword and obtain relevant informations
and stats related to your search queries. 

Thanks to the Streamlit Community for amazing [`Streamlit components`](https://discuss.streamlit.io/t/streamlit-components-community-tracker/4634)! 

## What else?
- Downlaod batch of articles 
- Check authors network 
- Associated keywords, yearly trends , published journals etc
- Downlaod Abstract for your selected publication
- Filtering Publication Table
- Get specific author's details (under implementation)

## Acknowledgements
### Modules for fetching Publication  
- [`Biopython`](https://biopython.org/docs/1.75/api/index.html)
    - [`Biopython Entrez Packgae`](https://biopython.org/docs/1.75/api/Bio.Entrez.html)
    (Please refer to this page for details) 
    - [`NCBI: Entrez Programming Utilities Help`](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [`Scholarly`](https://scholarly.readthedocs.io/en/latest/index.html)
### Streamlit Components    
- [`Streamlit-AgGrid`](https://github.com/PablocFonseca/streamlit-aggrid)
- [`Streamlit-Pyvis`](https://github.com/napoles-uach/stvis)
- [`Streamlit-Lottie`](https://github.com/andfanilo/streamlit-lottie)

### Blogs to refer 
- [`Network analysis to quickly get insight into an academic field with python`](https://towardsdatascience.com/network-analysis-to-quickly-get-insight-into-an-academic-field-with-python-cd891717d547)
- [`Streamlit-AgGrid usage`](https://towardsdatascience.com/7-reasons-why-you-should-use-the-streamlit-aggrid-component-2d9a2b6e32f0)

*Note : The App is still under devlelopment.*

'''
menu_items = {'About': about_info}
st.set_page_config(layout="wide",menu_items=menu_items,page_title = "SciLit 📑")

url = "https://assets2.lottiefiles.com/temp/lf20_TOE9MF.json"
lottie_url = url
lottie_json = load_lottieurl(lottie_url)
st_cont = st.sidebar.container()
st.sidebar.title("📖 SciLit ")
with st_cont:
    st_lottie(lottie_json,height =150,width =200)
    
# Change this email to your email address
Entrez.email = "A.N.Other@example.com"
st.sidebar.write("Search scientific publications based on keywords.")
# Separate based on Databases
st.sidebar.markdown("### Database ?")
database_choice = st.sidebar.radio(
                label = "💬 Select relevant database for best results",
                options=["PubMed","arXiv","bioRxiv"])
if database_choice == "PubMed":
    st.sidebar.subheader("Search your queries below")
    keyword = st.sidebar.text_input("For Example : Chlorophyll f ")

    if keyword :    
        with st.spinner("Searching query {}....".format(keyword)):
            max_pub = avail_pub(keyword) 
            
            if round(int(max_pub)/2) > 50:
                nSearchQuant = 50
            else:
                nSearchQuant = round(int(max_pub)/2)  

        if int(max_pub) > 0:
            st.sidebar.code("Total hits: "+max_pub)   
            val = st.sidebar.number_input(label = "Number of Publications to Load", min_value = 0, max_value = int(max_pub),value = nSearchQuant )       
            search_val = st.sidebar.checkbox("Start Loading",value = False)
            # Search all the publications
            if search_val:
                with st.spinner("Loading publication..."):
                    fc1, fc2 , fc3 , fc4 = st.columns(4)
                    _cont_tab = st.expander("Publications",expanded = True)
                    with _cont_tab:
                        pubdf , df = get_pub(val,keyword)            
                        res,sel = grid_table(df)
                    with fc2:
                    #url = "https://assets8.lottiefiles.com/packages/lf20_XhrXqc.json"
                        url = "https://assets1.lottiefiles.com/packages/lf20_2sNufo.json"
                        lottie_json = load_lottieurl(url)
                        st_lottie(lottie_json,height =100,width =200)
                    with fc1:
                        st.header(keyword)
                        
                    with fc1:
                        st.code("Search Results ✔️")
                        

                        

                csv = convert_df(df)
                st.download_button(
                    label ="Download Table 📋",
                    data = csv,
                    file_name = keyword+'.csv',
                    mime ='text/csv',
                )

                if sel:
                    exp = st.expander("More information on seleted article(s).")
                    with exp:
                        ids = [] 
                        nID = len(sel)
                        for x in range(nID):
                            ids.append(sel[x]["PubchemID"])       
                        match_df = pubdf[pubdf['PMID'].isin(ids)]               
                        nl = len(match_df)
                        
                        Title = match_df["TI"].to_list()
                        Abstract = match_df["AB"].to_list()
                        PublishingDate =  match_df["DP"].to_list()
                        Journal = match_df["SO"].to_list() 
                        auth = match_df["AU"].to_list()
                        doi = match_df["LID"].to_list()
                        pmid = match_df["PMID"].to_list()
                        
                        for x in range(nl):
                            dic_download = {
                                'Title':Title[x],
                                'Author(s)':auth[x],
                                'Journal Details':Journal[x],
                                'Abstract':Abstract[x]
                            }
                            res = json.dumps(dic_download)
                            st.success(Title[x])
                            st.markdown(auth[x])
                            st.write('Published on: ' +  Journal[x] )
                            st.info(Abstract[x])
                            st.download_button(label = 'Download Content 🗞',data = res,file_name = str(pmid[x])+'.txt', mime='text/csv')                          
                with st.container():
                    _exp = st.expander(label = 'Authors Network')
                    with _exp:
                        g,names = plot_connection(pubdf)
                        choice = st.selectbox("Get Author's Details (Note:may not work in all case!)",list(names["Name"]))
                        
                        disp = st.button("Display")
                        if disp:
                            try:
                                auth = author_information(choice)
                                st.write(auth)
                            except:
                                st.error("Sorry, not found!")
                                #url = "https://assets1.lottiefiles.com/packages/lf20_OT15QW.json"
                                url = "https://assets8.lottiefiles.com/packages/lf20_cwGCWK.json"
                                lottie_json = load_lottieurl(url)
                                st_lottie(lottie_json,height =250,width =250)
                    _exp_plots = st.expander(label = 'Plots  ',expanded = True)
                    with _exp_plots:
                        container_plots = st.container()
                        with container_plots:
                            fig = plot_top_authors(pubdf)
                            fig2 = plot_top_journals(pubdf)
                            fig3 = top_keyworkds(pubdf)
                            fig4 = year_journal_trend(pubdf)
                            st.plotly_chart(fig)
                            st.plotly_chart(fig4)          
                            st.plotly_chart(fig2)
                            try:
                                st.plotly_chart(fig3) 
                            except:
                                st.error('Keywords Stats Unavailable')
                                        
        else:
            st.error("Sorry!Your keyword(s) doesn't contain related publications.")
            url = "https://assets8.lottiefiles.com/temp/lf20_ZGnXlB.json"
            lottie_json = load_lottieurl(url)
            st_lottie(lottie_json,height =150,width =200)
    else:
        w_cont = st.sidebar.container()
        with w_cont:
            url = "https://assets1.lottiefiles.com/packages/lf20_OT15QW.json"
            lottie_json = load_lottieurl(url)
            st_lottie(lottie_json,height =250,width =250)
elif database_choice == "arXiv":
    st.info("Database search not implemented yet!")
else:
    st.info("Database search not implemented yet!")    

# For all Databases, when empty throw waiting cat lottie
#st.sidebar.markdown("### Hi, there!")
end_cont = st.sidebar.container()
with end_cont:  
    '''
    I'm Avra ! *Thanks for visiting my simple app, I'd love feedback on this,*
    *so if you want to reach out or support me, you can find me on - * &nbsp[![Follow](https://img.shields.io/twitter/follow/Avra_b?style=social)](https://www.twitter.com/Avra_b)
    [![Buy me a coffee](https://img.shields.io/badge/Buy%20me%20a%20coffee--yellow.svg?logo=buy-me-a-coffee&logoColor=orange&style=social)](https://www.buymeacoffee.com/AvraCodes) 

    '''
    st.markdown("<br>", unsafe_allow_html=True)


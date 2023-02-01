# Modules
import streamlit as st
import json
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
# Custom functions implementation
from funx import custom_funx, database_bio, plot_utils
import openai

def generate_response(token, prompt , output, temparature):
    openai.api_key = token
    response = openai.Completion.create(
                engine = "text-davinci-003",
                prompt = prompt,
                max_tokens = output,
                temperature = temparature,
            )
            # Print the generated summary
    res = response["choices"][0]["text"]
    return res

# Main App 
st.set_page_config(layout="wide",
                page_title = "PubLit 📑 - ✨ Scientific searches & insights all-in-one place.")
st_cont = st.sidebar.container()
st.sidebar.title("📖 PubLit ")
st.sidebar.markdown("✨ Scientific searches & insights all-in-one place.")


st.info("Search scientific publications based on keywords.", icon="🔎")
keyword = st.sidebar.text_input(label = ":blue[Search your queries below] ", 
                                placeholder= "Chlorophyll f")
if keyword :    
    with st.spinner("Searching query {}....".format(keyword)):
    # Obtain the maximum publications
        max_pub = database_bio.avail_pub(keyword) 
    # Default Serach Quantity Choosen from app
    if round(int(max_pub)/2) > 50:
        nSearchQuant = 50
    else:
        nSearchQuant = round(int(max_pub)/2)  

    if int(max_pub) > 0:
        st.sidebar.code("Total hits: "+max_pub)   
        val = st.sidebar.number_input(label = ":blue[Number of Publications to Load]", min_value = 0, max_value = int(max_pub),value = nSearchQuant )       
        search_val = st.sidebar.checkbox("⏳ Start Loading",
                                        value = False,
                                        help= "The App uses PubMed database to retrieve publications.")
        # Search all the publications
        if search_val:
            with st.spinner("Loading publication..."):
                pubdf , df = database_bio.get_pub(val,keyword)  
                with st.expander(label=f"{keyword}: Search Results", expanded=True):
                    res,sel = plot_utils.grid_table(df)   
                    csv = custom_funx.convert_df(df)
                    st.download_button(
                        label ="📋 Download Table",
                        data = csv,
                        file_name = keyword+'.csv',
                        mime ='text/csv',
                    ) 
            if st.sidebar.checkbox("📊 Visualize",value = False):
            # Figures            
                with st.expander("Visualize", expanded=True):
                    fig = plot_utils.plot_top_authors(pubdf)
                    fig2 = plot_utils.plot_top_journals(pubdf)
                    fig3 = plot_utils.top_keyworkds(pubdf)
                    fig4 = plot_utils.year_journal_trend(pubdf)
                    tab1, tab2, tab3, tab4, tab5 = st.tabs(["Top Author(s)",
                                                    "Journal(s) Trend", 
                                                    "Top Keyword(s)", 
                                                    "Top Journal(s)",
                                                    # "Journal(s) Trend",
                                                    "Author(s) Network"])
                    tab1.plotly_chart(fig, theme='streamlit')
                    tab2.plotly_chart(fig4, theme='streamlit')          
                    tab3.plotly_chart(fig2, theme='streamlit')
                    with tab5:
                        g,names = plot_utils.plot_connection(pubdf)
                    try:
                        tab4.plotly_chart(fig3, theme='streamlit') 
                    except:
                        tab4.error('Keywords Stats Unavailable')
                
            if sel:
                ids = [] 
                nID = len(sel)
                for x in range(nID):
                    ids.append(sel[x]["PubchemID"])       
                match_df = pubdf[pubdf['PMID'].isin(ids)]   # Get the matching rows           
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
                    st.sidebar.download_button(label = '🗞 Download Selected Content',  
                                    data = res,file_name = str(pmid[x])+'.txt',
                                    mime='text/csv',
                                    help = Journal[x] + " : [click to download].")   
                with st.form("AI support"):
                    st.subheader(Title[x])
                    st.markdown(' , '.join(auth[x]))
                    st.markdown('Published on: ' +  Journal[x] )
                    st.info(Abstract[x])                   
                    opt = st.radio(":blue[🆕 OpenAI Support ✨]",
                                ["Article Background", 'Summarize Abstract'],
                                horizontal=True, 
                                help = "How to use? Refer to this article - "+"https://medium.com/@avra42/summarizing-scientific-articles-with-openai-and-streamlit-fdee12aa1a2b")
                    
                    st.markdown(
                        '''
                        ``` 
                        Article Background : The app uses the openAI model to look for
                        more relevant information based on the article's title. The output
                        is likely to be obtained in point-wise format with relevant citations provided.
                        
                        Summarize Abstract : The app uses the openAI model to provide a concise 
                        summary of the abstract in point-wise format.
                        
                        Open AI's API Secret : You need to use your own API key from [here](https://platform.openai.com/account/api-keys).
                        The app **doesnot** store any sorts of information. But still, you are highly recommended to delete your used API keys after usage.
                        ```
                        
                        '''
                        
                        
                    )
                    
                    secret = st.text_input(":blue[Open AI's API Secret]",
                                    type="password",
                                    help = "Your API key is NOT stored in any form. However, highly recommended to delete once used."
                                    )
                    temp = st.slider("Relevance",value = 0.7,min_value=0.0,max_value=1.0,step=0.1 )
                    output_size = st.slider("Output Length",value=1024, min_value=64, max_value=4000)
                    if opt == "Article Background":    
                        fprompt = f" Give relevant scientific background on this topic: ' {Title[x]} ' bullet point-wise. Add scientific citations wherever necessary."  
                    else :
                        fprompt = f" Give scientific summarization on this topic: ' {Abstract[x]} ' bullet point-wise."  

                    submit_btn = st.form_submit_button("Submit",
                                    type = "primary")
                    if submit_btn :
                        res = generate_response(token = secret, 
                                                prompt = fprompt , 
                                                output = output_size,
                                                temparature = temp)
                        st.markdown('''------------''')
                        st.markdown(res)
                        st.warning('''
                                   Disclaimer: Results are AI generated by the openAI model, 
                                   cannot guarantee scientific relevance in all instance.
                                   Be cautious in further usage.
                                   ''', icon="🤖"
                            
                        )
            else:
                st.info("Select any article to know more.",icon="ℹ️")
                                                      
    else:
        st.error("Sorry!Your keyword(s) doesn't contain related publications.")


st.sidebar.markdown(
    '''
    --------------------
    Incase of feedbacks, you can DM me over [![Twitter](https://img.shields.io/badge/Twitter-%231DA1F2.svg?style=for-the-badge&logo=Twitter&logoColor=white)](https://www.twitter.com/Avra_b) or write to me at : avrab.yt@gmail.com
    '''
)
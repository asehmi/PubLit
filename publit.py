# Modules
import streamlit as st
import json
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
# Importing custom modules
from utils import custom_funx, database_bio, plot_utils
from utils.feedback import feedback, embedTweet
import time
from streamlit.components.v1 import html

def load_started():
    st.session_state.load_pub = True
def change_load_status():
    st.session_state.load_pub = False
# Main App 
st.set_page_config(layout="wide",page_title = "📑 PubLit AI - ✨ Scientific searches & insights all-in-one place.")

# Initiate SessionState
if 'load_pub' not in st.session_state:
    st.session_state['load_pub'] = False
if 'load_button_type' not in st.session_state:
    st.session_state.load_button_type = 'secondary'

if st.session_state['load_pub']:
    st.session_state.load_button_type = 'secondary'
else:
    st.session_state.load_button_type = 'primary'

st.sidebar.title("📖 PubLit AI ")
st.info('''
        ✨ Scientific searches & insights all-in-one place  🆕 AI based summarization 
        🛢arXiv database support coming soon
        ''')

with st.sidebar:
    st.markdown(''' 🔎 Search scientific publications based on keywords. ''')
    # Database
    # selected_datbase = st.radio(label = ":blue[Select database]", options =["⚛ PubMed", "ㄨArXiv"], horizontal=True)
    selected_datbase = "⚛ PubMed"
    # Seach / keywords Query 
    keyword = st.text_input(label = ":blue[Search your queries below] ", 
                                placeholder= "Example : Chlorophyll f", 
                                on_change = change_load_status
                                )

# When Keyword / Search query is present
if keyword :    
    with st.spinner("Searching query {}....".format(keyword)):
        if selected_datbase == "⚛ PubMed" :
            # Obtain the maximum publications available
            max_pub = database_bio.avail_pub(keyword)           
            # Total hits can be deprecated 
            if int(max_pub) > 0:               
                st.sidebar.code("Total hits: "+max_pub)  
        else: # Database ArXiv
            # Impement search for keywords less than 10
            max_pub = 50 
    # Default Serach Quantity Choosen from app
    if round(int(max_pub)/2) > 50:
        nSearchQuant = 50
    else:
        nSearchQuant = round(int(max_pub)/2)  

    if int(max_pub) > 0:  
        val = st.sidebar.number_input(label = ":blue[Number of Publications to Load]", min_value = 0, max_value = int(max_pub),value = nSearchQuant )       
        # search_val = st.sidebar.checkbox("⏳ Start Loading", value = st.session_state.load_pub,
                                        # help = "The App uses PubMed database to retrieve publications.")
        search_val = st.sidebar.button("⏳ Load Publications", on_click=load_started,type = st.session_state.load_button_type,
                                        help = "The App uses PubMed database to retrieve publications.")
        # Search all the publications
        if search_val or st.session_state.load_pub:
            with st.spinner("Loading publication..."):
                pubdf , df = database_bio.get_pub(val,keyword)  
                with st.expander(label=f"Search Results : {keyword}", expanded=True):
                    t1 , t2 = st.tabs(["Publications", "Visualize"])
                    # Tab with Data Frame
                    with t1:
                        res,sel = custom_funx.grid_table(df)  
                        # Button to dowload table
                        # Need some pause, because AgGrid takes time to dump data
                        time.sleep(0.5)
                        csv = custom_funx.convert_df(df)
                        st.download_button(
                            label ="📋 Download Table",
                            data = csv,
                            file_name = keyword+'.csv',
                            mime ='text/csv',
                        )
                        info_box_sel = st.empty()
                    
                    # t2 -> Tab with Figures            
                    with t2:
                        fig = plot_utils.plot_top_authors(pubdf)
                        fig2 = plot_utils.plot_top_journals(pubdf)
                        fig3 = plot_utils.top_keyworkds(pubdf)
                        fig4 = plot_utils.year_journal_trend(pubdf)
                        tab1, tab2, tab3, tab4, tab5 = st.tabs(["Top Author(s)",
                                                        "Journal(s) Trend", 
                                                        "Top Journal(s)",
                                                        "Top Keyword(s)", 
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
            
            # Any selection on the AgGrid data table   
            if sel:
                # This can be rewritten for better implementation
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
                                    help = 'Content of:' + Journal[x] + " : [click to download].")   
                with st.form("ai_support"):
                    st.subheader(Title[x])
                    st.markdown(' , '.join(auth[x]))
                    # l = 'https://doi.org/' + match_df['LID'].str.replace(' \[doi\]', '')
                    # st.markdown(l.to_list())
                    # match_df["SO"]
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
                    
                    secret = st.text_input(":blue[OpenAI API Key]",
                                    type="password",
                                    # value = st.session_state.user_secret,
                                    placeholder="Paste your OpenAI API key here (sk-...)",
                                    help = '''You can get your API key from https://platform.openai.com/account/api-keys.
                                            Your API key is NOT stored in any form. 
                                            However, highly recommended to delete once used.
                                            '''
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
                        empty_info = st.empty()
                        
                        if not secret or not secret.startswith('sk-'):
                            st.warning("Invalid : OpenAI API key not found.")
                            st.markdown(
                                    "#### Usage:\n"
                                    "1. Get your [OpenAI API key](https://platform.openai.com/account/api-keys) from here.\n"
                                    "2. Paste your OpenAI API key above in the `OpenAI API Key` text box.\n"
                                    "3. Change the parameters based on your needs.\n"
                                    "4. Learn more about the parameters in this [blog](https://medium.com/@avra42/summarizing-scientific-articles-with-openai-and-streamlit-fdee12aa1a2b)." 
    )
                            
                        else:
                            empty_info.text("🤖 AI magic initiated 🪄 ...") 
                            with st.spinner("Performing `{}` method.".format(opt)):
                                res = custom_funx.generate_response(token = secret, 
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
                            empty_info.empty()
            else:
                info_box_sel.info("Select any publication to know more.",icon="ℹ️")                                                  
    else:
        st.error("Sorry!Your keyword(s) doesn't contain related publications.")

else: # When no keyword is present in the search box widget
    st.text("⬅️ Search your query using the search box present in the sidebar."
            )

    st.sidebar.markdown(
        '''
        --------------------  
        Incase of feedbacks, you can DM me over [![Twitter](https://img.shields.io/badge/Twitter-%231DA1F2.svg?style=for-the-badge&logo=Twitter&logoColor=white)](https://www.twitter.com/Avra_b) or write to me at : avrab.yt@gmail.com
        '''
    )
    tweet_url = "https://twitter.com/Avra_b/status/1620606960709300224?s=20&t=34ICyDyZIx6RE0HxB1Lhlw"
    st.markdown("#### For a quick-start refer to my tweet [here]({}).".format(tweet_url))
    res = embedTweet(tweet_url=tweet_url)
    html(res,height=700)
with st.sidebar:
    st.markdown('')
    emailID = st.text_input(label=":blue[Notify me with future release.]",
                        placeholder="Please enter your email ID here.")
    if st.button("📩 Subscribe"):
        if len(emailID) > 0:
            feedback(email = emailID, sheetname = "PubLitAI")
# button(username="AvraCodes", floating=True, bg_color="#BCD1E6", width=220)


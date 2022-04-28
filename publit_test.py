# Modules
import streamlit as st
from streamlit_lottie import st_lottie
import json
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# Custom functions implementation
from publit import about_the_app,custom_funx,plot_utils,database_bio,database_scholarly

# Main App 

about_info = about_the_app.text_about()
menu_items = {'About': about_info}
st.set_page_config(layout="wide",menu_items=menu_items,page_title = "SciLit 📑")
lottie_url = "https://assets2.lottiefiles.com/temp/lf20_TOE9MF.json"
lottie_json = custom_funx.load_lottieurl(lottie_url)
st_cont = st.sidebar.container()
st.sidebar.title("📖 SciLit ")
with st_cont:
    st_lottie(lottie_json,height =150,width =200)
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
            # Obtain the maximum publications
            max_pub = plot_utils.avail_pub(keyword) 
            # Default Serach Quantity Choosen from app
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
                        pubdf , df = database_bio.get_pub(val,keyword)            
                        res,sel = plot_utils.grid_table(df)
                    with fc2:
                    #url = "https://assets8.lottiefiles.com/packages/lf20_XhrXqc.json"
                        url = "https://assets1.lottiefiles.com/packages/lf20_2sNufo.json"
                        lottie_json = custom_funx.load_lottieurl(url)
                        st_lottie(lottie_json,height =100,width =200)
                    with fc1:
                        st.header(keyword)   
                    with fc1:
                        st.text("Search Results ✔️")
                
                csv = custom_funx.convert_df(df)
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
                        g,names = database_bio.plot_connection(pubdf)
                        choice = st.selectbox("Get Author's Details (Note:may not work in all case!)",list(names["Name"]))
                        
                        disp = st.button("Display")
                        if disp:
                            try:
                                auth = database_scholarly.author_information(choice)
                                st.write(auth)
                            except:
                                st.error("Sorry, not found!")
                                #url = "https://assets1.lottiefiles.com/packages/lf20_OT15QW.json"
                                url = "https://assets8.lottiefiles.com/packages/lf20_cwGCWK.json"
                                lottie_json = custom_funx.load_lottieurl(url)
                                st_lottie(lottie_json,height =250,width =250)
                    _exp_plots = st.expander(label = 'Plots  ',expanded = True)
                    with _exp_plots:
                        container_plots = st.container()
                        with container_plots:
                            fig = plot_utils.plot_top_authors(pubdf)
                            fig2 = plot_utils.plot_top_journals(pubdf)
                            fig3 = plot_utils.top_keyworkds(pubdf)
                            fig4 = plot_utils.year_journal_trend(pubdf)
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
            lottie_json = custom_funx.load_lottieurl(url)
            st_lottie(lottie_json,height =150,width =200)
    else:
        w_cont = st.sidebar.container()
        with w_cont:
            url = "https://assets1.lottiefiles.com/packages/lf20_OT15QW.json"
            lottie_json = custom_funx.load_lottieurl(url)
            st_lottie(lottie_json,height =250,width =250)
elif database_choice == "arXiv":

    st.info("arXiv database search not implemented yet!")
    st.sidebar.subheader("Search your queries below")
    keyword = st.sidebar.text_input("For Example : Quantum")
else:
    st.info("bioRxiv database search not implemented yet!")
    st.sidebar.subheader("Search your queries below")
    keyword = st.sidebar.text_input("For Example : CRISPAR CAS9")  

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


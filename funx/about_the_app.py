import streamlit as st
def text_about(): 
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
    return about_info 

def foot_note():
    
    note = '''
    *I'm Avra ! Thanks for visiting my simple app, I'd love feedback on this,
    so if you want to reach out or support me, you can find me here*[![Follow](https://img.shields.io/twitter/follow/Avra_b?style=social)](https://www.twitter.com/Avra_b)
    [![Buy me a coffee](https://img.shields.io/badge/Buy%20me%20a%20coffee--yellow.svg?logo=buy-me-a-coffee&logoColor=orange&style=social)](https://www.buymeacoffee.com/AvraCodes) 
    '''
    return note 
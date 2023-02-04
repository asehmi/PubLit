import arxiv
import streamlit as st
import pandas as pd
import re
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import GridOptionsBuilder
import ssl
ssl._create_default_https_context = ssl._create_unverified_context


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
            "DOI":doi,
            }
    return pd.DataFrame(dic)

s = st.text_input(label='Search', value="phycobilisomes")
df = arxiv_search(s, 10)
df['Authors'] = df['Authors'].astype(str)
df['Authors'] = df['Authors'].apply(lambda x: re.findall(r"\'(.*?)\'", x))
df['Published on'] = pd.to_datetime(df['Published on'])
df['Published on'] = df['Published on'].dt.strftime('%Y %b %d')
# Add 'https://doi.org/' in front of each row of the 'DOI' column
# df['DOI'] = 'https://doi.org/' + df['DOI']
gb = GridOptionsBuilder.from_dataframe(df)
gb.configure_default_column(editable=True)
gb.configure_selection('single', use_checkbox=True)
gridOptions = gb.build()

grid_response = AgGrid(df, gridOptions = gridOptions)
selected_rows = grid_response["selected_rows"]
if selected_rows:   
    url = selected_rows[0]['arxivID']
    selected_rows[0]['Title']
    id_string = url.split("/")[-1]
    id_string
    paper = next(arxiv.Search(id_list=[id_string]).results())
    
    if st.button("Download"):
        dFilename = id_string + ".pdf"
        info = st.empty()
        info.info(f"Downloading {paper} ...")
        paper.download_pdf(filename= id_string + ".pdf")

# paper = next(arxiv.Search(id_list=["1605.08386v1"]).results())
# paper.download_source(filename="download.pdf")

        import base64
        info.info(f"Trying to encod {paper} ...")
        with open(dFilename, "rb") as file:
            pdf_data = file.read()
            encoded_pdf = base64.b64encode(pdf_data).decode("utf-8")

        # pdf_data = base64.b64decode(encoded_pdf.encode("utf-8"))
        # print(pdf_data)
        # pdf_data
        # encoded_pdf
        # In your Streamlit front-end, you can use the encoded data to display the PDF.
        # html_template = """
        # <iframe
        #     src="data:application/pdf;base64,{encoded_pdf}"
        #     width="100%"
        #     height="500"
        # >
        # </iframe>
        # """

        # st.markdown(html_template.format(encoded_pdf=encoded_pdf), unsafe_allow_html=True)   
        pdf_display = F'<iframe src="data:application/pdf;base64,{encoded_pdf}" width="700" height="1000" type="application/pdf"></iframe>'
        info.info("Embedding PDF {dFilename} now ...")
        st.markdown(pdf_display, unsafe_allow_html=True)
        info.empty()
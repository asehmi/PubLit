# PubLit
Publit provides you to search for scientific publictions with specific keywords and obtain useful anlaytics related to the :
- authors
- journals, 
- author's network  
- yearly trends and many more

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/avratanubiswas/publit/main/publit.py)

## How it works ?
PubLit uses, [PubMed](https://pubmed.ncbi.nlm.nih.gov/) database to search for any user-defined keywords and grasp the related informations. Finally the information extracted is thrown to the front-end. Subsequent analysis are made with common data analysis library such as pandas and visulaized using plotly. Data anlaysis motivation behind this basic UI is based on this [article](https://towardsdatascience.com/network-analysis-to-quickly-get-insight-into-an-academic-field-with-python-cd891717d547). You will find more useful references related to the Streamlit UI in the **About** section of the application.

## Future plans?
Well, I have recieved a couple of insightful comments and feedbacks which I would love to work on in my spare time. For example, polishing the UI, charts, using semantischolar API or making the code a bit more universal for different databases.
>(Thanks to Nikolay Karelin for his valuable feedbacks !!)

Therefore, the plan is to integrate and take advantage of more such public databases, for example, [arXiv](https://arxiv.org/help/api/index) or [bioRxiv](https://api.biorxiv.org/), as a result making the application available for science-enthusiast from various backgrounds.

At this moment PubLit is in a work in progress and a fun project to work on. I'll be happy to recive feedbacks and collaborate with others. 



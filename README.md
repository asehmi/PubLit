# PubLit
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/avratanubiswas/publit/main/publit.py)
With PubLit, search for scientific publictions with specific keywords and obtaining useful anlaytics gets easier. What analytics does PubLit offer ?
- authors details
- journals related informations
- author's network  
- yearly trends and many more

*In future release, you can opt for specific datbase which suits your keywords.*

## How it works ?
PubLit currently uses, [PubMed](https://pubmed.ncbi.nlm.nih.gov/) database to search for any user-defined keywords and grasp the related informations. Finally the information extracted is thrown to the Streamlit front-end. Subsequent analysis are made with common data analysis library such as pandas and visulaized using plotly. Data anlaysis motivation behind this basic UI is based on this [article](https://towardsdatascience.com/network-analysis-to-quickly-get-insight-into-an-academic-field-with-python-cd891717d547). You will find more useful references related to the Streamlit UI in the *About* section of the application.

## Future plans?
Well, I have recieved a couple of insightful comments and feedbacks which I would love to work on during my spare time. For example, polishing the UI, plots, using semantischolar API or making the code a bit more universal for different databases.
>(Special mention to Nikolay Karelin for his valuable feedbacks. Thank you!!)

Therefore, the future immediate plan is to integrate and take advantage of more such public databases, for example, [arXiv](https://arxiv.org/help/api/index) or [bioRxiv](https://api.biorxiv.org/), as a result making the application available for science-enthusiast from various backgrounds.

At this moment PubLit is in a work in progress and a fun project to work on. I'll be happy to receive feedbacks and collaborate with others. Feel free to ping me here or in [twitter](https://twitter.com/Avra_b).



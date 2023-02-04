from scihub import SciHub
import streamlit as st
sh = SciHub()

# retrieve 5 articles on Google Scholars related to 'bittorrent'
result = sh.fetch('http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=1648853')
st.write(result.doi)
# download the papers; will use sci-hub.io if it must
# for paper in results['papers']:
# 	sh.download(paper['url'])
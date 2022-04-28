
import requests
import streamlit as st

# Convert Datframe for saving
@st.experimental_memo
def convert_df(df):
    return df.to_csv().encode('utf-8')

def load_lottieurl(url: str):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

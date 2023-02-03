import streamlit as st
import openai
# Convert Datframe for saving
@st.experimental_memo
def convert_df(df):
    return df.to_csv().encode('utf-8')

def generate_response(token, prompt , output, temparature):
    """
    Generates a response using OpenAI's API.
    
    Parameters
    ----------
    token : str
        The API key for OpenAI.
    prompt : str
        The text prompt to generate the response from.
    output : int
        The maximum number of tokens in the generated response.
    temperature : float
        A value controlling the randomness of the generated response.
        
    Returns
    -------
    str
        The generated response.
    """
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

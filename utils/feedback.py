from google.oauth2 import service_account
import gspread
import streamlit as st
import re
import requests


def feedback(email, sheetname):
    """
    This function takes an email address and a sheet name as input and writes the email address to a Google Sheet.

    Parameters:
    email (str): The email address to be written to the Google Sheet.
    sheetname (str): The name of the worksheet in the Google Sheet to which the email address will be written.

    Returns:
    None
    """
    # Create a Google Authentication connection object
    scope = ['https://spreadsheets.google.com/feeds',
            'https://www.googleapis.com/auth/drive']
    credentials = service_account.Credentials.from_service_account_info(
                    st.secrets["gcp_service_account"], scopes = scope)
    
    # Regular expression pattern for validating email addresses
    regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'         
    
    # Check if the email address is valid
    if re.fullmatch(regex, email):
        try:
            # Authorize the Google connection
            gc = gspread.authorize(credentials)
            
            # Open the Google Sheet
            sh = gc.open("Python101-course")
            
            # Get the specified worksheet
            worksheet = sh.worksheet(sheetname)
            
            # Get the values in the first column of the worksheet
            values_list = worksheet.col_values(1)
            
            # Get the number of rows in the worksheet
            length_row = len(values_list)
            
            # Write the email address to the next row in the first column of the worksheet
            worksheet.update_cell(length_row + 1, 1, email)
            
            # Display a success message
            st.info("Subscribed", icon="✅")
        except:
            # Display an error message if the operation failed for any unknown reason
            st.warning("Submission failed for unknown reason(s). Kindly DM over Twitter or write to here - avrab.yt@gmail.com", icon="🚨")
    else:
        # Display an error message if the email address is not valid
        st.error("Invalid email ID - `{}` error.".format(email), icon="⚠️")
        return

def embedTweet(tweet_url):
    api = "https://publish.twitter.com/oembed?url={}".format(tweet_url)
    response = requests.get(api)
    res = response.json()["html"]
    return res
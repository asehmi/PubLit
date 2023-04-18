import streamlit as st
import openai
from st_aggrid import AgGrid,GridUpdateMode, JsCode
from st_aggrid.grid_options_builder import GridOptionsBuilder

# Convert Datframe for saving
@st.cache_data
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

def grid_table(df):
    """
    This function takes a Pandas DataFrame as input and returns a modified DataFrame and a list of selected rows.

    Parameters:
    df (Pandas DataFrame): The input DataFrame.

    Returns:
    df (Pandas DataFrame): The modified DataFrame.
    selected_rows (list): A list of selected rows in the DataFrame.
    """
    # Create a GridOptionsBuilder object from the input DataFrame
    gb = GridOptionsBuilder.from_dataframe(df)

    # Configure pagination
    gb.configure_pagination(enabled=True, 
                            paginationAutoPageSize=True, 
                            paginationPageSize=5)

    # Configure the side bar
    gb.configure_side_bar()

    # Configure selection
    gb.configure_selection('single', use_checkbox=True)

    # Configure default column options
    gb.configure_default_column(groupable=True, 
                                value = True, 
                                enableRowGroup=True, 
                                editable= True,
                                tooltipField = 'Title'
                                )
    cell_renderer = JsCode(
        """
            function(params) {
                return '<a href=' + params.value + ' target="_blank">'+ params.value +'</a>'
                }
            """
    )
    gb.configure_column(
        "DOI",
        headerName="DOI",
        cellRenderer=cell_renderer,
    )
    # Build the grid options
    gridOptions = gb.build()

    # Create a grid response
    grid_response = AgGrid(df, 
                gridOptions = gridOptions, 
                enable_enterprise_modules = True,
                width = 500, height= 400,
                theme = "balham",
                update_mode = GridUpdateMode.SELECTION_CHANGED,
                allow_unsafe_jscode=True
                )

    # Get the modified DataFrame and the selected rows
    df = grid_response['data']
    selected_rows = grid_response["selected_rows"]

    # Return the modified DataFrame and the list of selected rows
    return df, selected_rows

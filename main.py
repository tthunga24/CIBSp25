import streamlit as st

st.set_page_config(page_title="My App", layout="centered")

st.title("📊 Options Strategy Tool")

# Textbox input
user_input = st.text_input("Enter stock ticker:")

# Dropdown / selectbox input
strategy = st.selectbox(
    "Choose a strategy:",
    ["Calendar Spread", "Straddle", "Iron Condor", "Butterfly", "Custom"]
)

# Display user input
st.write(f"Ticker entered: `{user_input}`")
st.write(f"Selected strategy: **{strategy}**")
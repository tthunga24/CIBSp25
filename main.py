import streamlit as st
import yfinance as yf
import datetime

st.set_page_config(page_title="Options Chain Viewer", layout="centered")

st.title("📈 Options Chain Lookup")

# Input for stock ticker
ticker_symbol = st.text_input("Enter stock ticker:", "AAPL")

# Dropdown for strategy (placeholder for now)
strategy = st.selectbox("Choose a strategy:", ["Calendar Spread", "Straddle", "Long Call", "Custom"])

# Fetch data when a ticker is entered
if ticker_symbol:
    ticker = yf.Ticker(ticker_symbol.upper())
    
    try:
        # Get all expirations
        expirations = ticker.options
        if not expirations:
            st.warning("No options data available for this ticker.")
        else:
            # Find expiration date ~3 weeks out
            target_date = datetime.date.today() + datetime.timedelta(weeks=3)
            exp_date = min(expirations, key=lambda d: abs(datetime.datetime.strptime(d, "%Y-%m-%d").date() - target_date))

            st.write(f"📅 Closest expiration to 3 weeks out: `{exp_date}`")

            # Get calls for that expiration
            opt_chain = ticker.option_chain(exp_date)
            calls = opt_chain.calls

            st.subheader("Call Options Chain")
            st.dataframe(calls)

    except Exception as e:
        st.error(f"Error fetching options data: {e}")

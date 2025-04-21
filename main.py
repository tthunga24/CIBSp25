import streamlit as st
import yfinance as yf
import datetime
import numpy as np
from scipy.stats import norm

# ----- Black-Scholes Greeks -----
def black_scholes_greeks(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2)*T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    delta = norm.cdf(d1)
    gamma = norm.pdf(d1) / (S * sigma * np.sqrt(T))
    theta = (-S * norm.pdf(d1) * sigma / (2 * np.sqrt(T))) - r * K * np.exp(-r * T) * norm.cdf(d2)
    vega = S * norm.pdf(d1) * np.sqrt(T)
    
    return delta, gamma, theta, vega

# ----- Streamlit UI -----
st.set_page_config(page_title="Options Chain Viewer", layout="centered")
st.title("📈 Options Chain Lookup with Greeks")

ticker_symbol = st.text_input("Enter stock ticker:", "AAPL")
strategy = st.selectbox("Choose a strategy:", ["Calendar Spread", "Straddle", "Long Call", "Custom"])

if ticker_symbol:
    ticker = yf.Ticker(ticker_symbol.upper())
    
    try:
        current_price = ticker.history(period="1d")["Close"].iloc[-1]
        st.write(f"📌 Current Stock Price: **${current_price:.2f}**")
        
        expirations = ticker.options
        if not expirations:
            st.warning("No options data available.")
        else:
            target_date = datetime.date.today() + datetime.timedelta(weeks=3)
            exp_date = min(expirations, key=lambda d: abs(datetime.datetime.strptime(d, "%Y-%m-%d").date() - target_date))
            st.write(f"📅 Closest expiration to 3 weeks out: `{exp_date}`")

            opt_chain = ticker.option_chain(exp_date)
            calls = opt_chain.calls

            available_strikes = calls['strike'].sort_values().tolist()
            selected_strike = st.selectbox("🎯 Select Short Leg (Call):", available_strikes)

            selected_option = calls[calls['strike'] == selected_strike]
            if not selected_option.empty:
                opt = selected_option.iloc[0]
                st.markdown("### Greeks:")

                # Inputs for BS model
                S = current_price
                K = opt['strike']
                T = (datetime.datetime.strptime(exp_date, "%Y-%m-%d").date() - datetime.date.today()).days / 365
                r = 0.05  # assume 5% risk-free rate
                sigma = opt['impliedVolatility'] if not np.isnan(opt['impliedVolatility']) else 0.4

                delta, gamma, theta, vega = black_scholes_greeks(S, K, T, r, sigma)

                st.write(f"**Strike:** {K}")
                st.write(f"**IV:** {sigma:.2%}")
                st.write(f"**Delta:** {delta:.4f}")
                st.write(f"**Gamma:** {gamma:.6f}")
                st.write(f"**Theta:** {theta / 365:.4f}")
                st.write(f"**Vega:** {vega / 100:.4f}")
            else:
                st.warning("Could not find selected option.")

            with st.expander("📂 View Full Call Chain Table"):
                st.dataframe(calls)

    except Exception as e:
        st.error(f"Error: {e}")

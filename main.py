import streamlit as st
import yfinance as yf
import datetime
import numpy as np
from scipy.stats import norm
from math import isclose

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
st.title("📈 Gamma Neutral Calendar Spread Builder")

ticker_symbol = st.text_input("Enter stock ticker:", "SPY")

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
            # Get current stock price
            current_price = ticker.history(period="1d")["Close"].iloc[-1]

            # Sort strikes and find the ATM strike (closest to spot)
            strikes_sorted = sorted(calls['strike'].tolist())
            atm_strike = min(strikes_sorted, key=lambda x: abs(x - current_price))

            # Find index of ATM strike
            atm_index = strikes_sorted.index(atm_strike)

            # Select ATM and 5 OTM strikes
            selected_strikes = strikes_sorted[atm_index:atm_index + 6]  # atm + 5 higher

            # Filter calls to only include these strikes
            filtered_calls = calls[calls['strike'].isin(selected_strikes)]

            # Dropdown to select from these
            selected_strike = st.selectbox("🎯 Select ATM or OTM Strike:", filtered_calls['strike'].sort_values())

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

                delta *= -1
                gamma *= -1
                theta *= -1
                vega *= -1

                st.write(f"**Strike:** {K}")
                st.write(f"**Price: {opt['lastPrice']}**")
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

if not selected_option.empty:
    if st.button("Generate Gamma-Neutral Calendar Spread"):
        st.session_state.generate_long_leg = True

if st.session_state.get("generate_long_leg", False):
    long_target = datetime.date.today() + datetime.timedelta(weeks=9)
    long_exp = min(expirations, key=lambda d: abs(datetime.datetime.strptime(d, "%Y-%m-%d").date() - long_target))
    st.write(f"📅 Long leg expiration selected: `{long_exp}`")

    long_chain = ticker.option_chain(long_exp).calls

    # Grab closest strike
    long_strikes = long_chain['strike'].tolist()
    closest_strike = min(long_strikes, key=lambda x: abs(x - selected_strike))
    long_option = long_chain[long_chain['strike'] == closest_strike]

    if not long_option.empty:
        long_opt = long_option.iloc[0]
        T_long = (datetime.datetime.strptime(long_exp, "%Y-%m-%d").date() - datetime.date.today()).days / 365
        sigma_long = long_opt['impliedVolatility'] if not np.isnan(long_opt['impliedVolatility']) else 0.4

        delta_L, gamma_L, theta_L, vega_L = black_scholes_greeks(
            S=current_price,
            K=closest_strike,
            T=T_long,
            r=0.05,
            sigma=sigma_long
        )

        

        st.markdown("### 🔍 Long Leg Option (≈45 Days)")
        st.write(f"**Closest Available Strike:** {closest_strike}")
        st.write(f"**Price:** {long_opt['lastPrice']}")
        st.write(f"**IV:** {sigma_long:.2%}")
        st.write(f"**Delta:** {delta_L:.4f}")
        st.write(f"**Gamma:** {gamma_L:.6f}")
        st.write(f"**Theta (per day):** {theta_L / 365:.4f}")
        st.write(f"**Vega:** {vega_L / 100:.4f}")
    else:
        st.warning("⚠️ Could not find a suitable strike for long leg.")

    gamma_short = gamma

    max_short_qty = 50
    tolerance = 0.005
    solution_found = False

    for q_s in range(1, max_short_qty + 1):
        q_l_exact = -q_s * (gamma_short / gamma_L)
        q_l_rounded = round(q_l_exact)

        if isclose(q_l_exact, q_l_rounded, abs_tol=tolerance) and q_l_rounded > 0:
            net_gamma = gamma_short * q_s + gamma_L * q_l_rounded
            net_delta = (delta * q_s) + (delta_L * q_l_rounded)
            net_vega = (vega * q_s) + (vega_L * q_l_rounded)
            net_theta = (theta * q_s) + (theta_L * q_l_rounded)

            st.markdown("### Final Gamma Neutral Position")
            st.write(f"**Short Calls:** {q_s}")
            st.write(f"**Long Calls:** {q_l_rounded}")
            st.write(f"**Net Delta:** {net_delta:.6f}")
            st.write(f"**Net Gamma:** {net_gamma:.6f}")
            st.write(f"**Net Theta:** {net_theta / 365:.6f}")
            st.write(f"**Net Vega:** {net_vega / 100:.6f}")
            st.write(f"**To neutralize delta, equity position:** {(net_delta * -100):.0f}")

            solution_found = True
            break

    r = 0.05
    S = selected_strike
    K_short = selected_strike
    K_long = closest_strike
    T_long = (datetime.datetime.strptime(long_exp, "%Y-%m-%d").date() - datetime.date.today()).days / 365
    qty_short = q_s
    qty_long = q_l_rounded
    dte_min = 1
    dte_max = int(T_long * 365)
    entry_cost = opt['lastPrice'] * qty_short + long_opt['lastPrice'] * qty_long

    # Define the Black-Scholes call formula in Mathematica syntax
    bs_call_mma = lambda S, K, T, r, sigma: (
        f"{S} * CDF[NormalDistribution[0, 1], "
        f"(Log[{S}/{K}] + ({r} + 0.5 * {sigma}^2) * {T}) / ({sigma} * Sqrt[{T}])] "
        f"- {K} * Exp[-{r} * {T}] * CDF[NormalDistribution[0, 1], "
        f"(Log[{S}/{K}] + ({r} - 0.5 * {sigma}^2) * {T}) / ({sigma} * Sqrt[{T}])]"
    )

    # Build full expression for spread PnL in terms of iv and t
    short_expr = bs_call_mma(S, K_short, "(dte/365)", r, "iv")
    long_expr = bs_call_mma(S, K_long, T_long, r, "iv")
    spread_expr = f"({qty_long} * ({long_expr}) - {qty_short} * ({short_expr})) - {entry_cost}"

    # Final Mathematica command string
    plot_cmd = (
    f'Plot3D[{spread_expr}, '
    f'{{iv, 0.1, 0.8}}, '
    f'{{dte, {dte_min}, {dte_max}}}, '
    f'PlotRange -> Full, '
    f'AxesLabel -> {{"IV", "Days to Expiry", "PnL"}}, '
    f'ColorFunction -> "TemperatureMap", Mesh -> None]'
    )

    st.code(plot_cmd)

    if st.button("Final position and visualize"):
        st.session_state.show_vis = True

    if not solution_found:
        st.warning("⚠️ Couldn't find integer quantities to gamma hedge within 50 short calls.")
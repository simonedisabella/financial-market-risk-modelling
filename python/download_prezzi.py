import yfinance as yf
import pandas as pd

tickers = ["AAPL", "MSFT", "GOOGL", "NVDA", "JPM", "GS", "XOM", "KO", "GLD", "TLT"]

prezzi = yf.download(tickers, end="2026-01-01", start="2020-01-01")

prezzi_chiusura = prezzi["Close"]
print(prezzi_chiusura)

prezzi_chiusura.to_csv("prezzi_chiusura.csv")

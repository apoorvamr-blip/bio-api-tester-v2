import streamlit as st
import requests

st.set_page_config(
    page_title="Bio API Tester",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Bio API Tester")
st.write("Enter a Locus ID to retrieve biological sequence information.")

locus_id = st.text_input(
    "Locus ID",
    placeholder="Example: At1g01010"
)

if locus_id:
    st.success(f"Locus ID entered: {locus_id}")

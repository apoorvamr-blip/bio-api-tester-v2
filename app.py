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

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        "db": "gene",
        "term": f"{locus_id}[Gene Name]",
        "retmode": "json"
    }

    response = requests.get(url, params=params)

    if response.status_code == 200:

        data = response.json()

        ids = data["esearchresult"]["idlist"]

       if len(ids) == 1:
    gene_id = ids[0]
    st.success(f"NCBI Gene ID: {gene_id}")

    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    summary_params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }

    summary_response = requests.get(
        summary_url,
        params=summary_params
    )

    if summary_response.status_code == 200:

        summary_data = summary_response.json()
        gene_data = summary_data["result"][gene_id]

        st.subheader("Gene Information")

        st.write("Gene name:", gene_data.get("name"))
        st.write("Description:", gene_data.get("description"))
        st.write("Organism:", gene_data.get("organism", {}).get("scientificname"))
        st.write("Chromosome:", gene_data.get("chromosome"))
        st.write("Taxonomy ID:", gene_data.get("taxid"))

    else:
        st.error("Could not retrieve NCBI gene information."))

        elif len(ids) == 0:
            st.error("No NCBI Gene record found.")

        else:
            st.warning("Multiple NCBI Gene records found.")

    else:
        st.error("NCBI request failed.")

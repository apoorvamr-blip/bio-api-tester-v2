import streamlit as st
import requests

# -----------------------------
# Page configuration
# -----------------------------

st.set_page_config(
    page_title="Bio API Tester",
    page_icon="🧬",
    layout="wide"
)

# -----------------------------
# App title
# -----------------------------

st.title("🧬 Bio API Tester")
st.write(
    "Enter a Locus ID to retrieve biological sequence information."
)

# -----------------------------
# User input
# -----------------------------

locus_id = st.text_input(
    "Locus ID",
    placeholder="Example: At1g01010"
)

# -----------------------------
# NCBI lookup
# -----------------------------

if locus_id:

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        "db": "gene",
        "term": f"{locus_id}[Gene Name]",
        "retmode": "json"
    }

    response = requests.get(
        url,
        params=params
    )

    if response.status_code == 200:

        data = response.json()

        ids = data["esearchresult"]["idlist"]

        # -----------------------------
        # Exactly one NCBI result
        # -----------------------------

        if len(ids) == 1:

            gene_id = ids[0]

            st.success(
                f"NCBI Gene ID: {gene_id}"
            )

            # -----------------------------
            # Get NCBI gene summary
            # -----------------------------

            summary_url = (
                "https://eutils.ncbi.nlm.nih.gov/"
                "entrez/eutils/esummary.fcgi"
            )

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

                # -----------------------------
                # DEBUG: Show NCBI response
                # -----------------------------

                st.write("DEBUG - NCBI gene data:")
                st.json(gene_data)

                # -----------------------------
                # Gene information
                # -----------------------------

                st.subheader("Gene Information")

                gene_name = gene_data.get(
                    "name",
                    "Not available"
                )

                description = gene_data.get(
                    "description",
                    "Not available"
                )

                chromosome = gene_data.get(
                    "chromosome",
                    "Not available"
                )

                organism = gene_data.get(
                    "organism",
                    {}
                )

                organism_name = organism.get(
                    "scientificname",
                    "Not available"
                )

                taxonomy_id = (
                    gene_data.get("taxid")
                    or gene_data.get("tax_id")
                    or "Not available"
                )

                st.write(
                    "Gene name:",
                    gene_name
                )

                st.write(
                    "Description:",
                    description
                )

                st.write(
                    "Organism:",
                    organism_name
                )

                st.write(
                    "Chromosome:",
                    chromosome
                )

                st.write(
                    "Taxonomy ID:",
                    taxonomy_id
                )

            else:

                st.error(
                    "Could not retrieve NCBI gene information."
                )

        # -----------------------------
        # No results
        # -----------------------------

        elif len(ids) == 0:

            st.error(
                "No NCBI Gene record found for this Locus ID."
            )

        # -----------------------------
        # Multiple results
        # -----------------------------

        else:

            st.warning(
                "Multiple NCBI Gene records were found. "
                "The Locus ID is ambiguous."
            )

    else:

        st.error(
            "NCBI request failed."
        )

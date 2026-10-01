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

                taxonomy_id = organism.get(
                    "taxid",
                    "Not available"
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

                # -----------------------------
                # Link Gene to nucleotide records
                # -----------------------------

                st.subheader("Nucleotide Records")

                link_url = (
                    "https://eutils.ncbi.nlm.nih.gov/"
                    "entrez/eutils/elink.fcgi"
                )

                link_params = {
                    "dbfrom": "gene",
                    "db": "nuccore",
                    "id": gene_id,
                    "retmode": "json"
                }

                link_response = requests.get(
                    link_url,
                    params=link_params
                )

                if link_response.status_code == 200:

                    link_data = link_response.json()

                    try:

                        nucleotide_ids = (
                            link_data["linksets"][0]
                            ["linksetdbs"][0]
                            ["links"]
                        )

                        st.write(
                            f"Found {len(nucleotide_ids)} "
                            "linked nucleotide records."
                        )

                        # -----------------------------
                        # Get nucleotide summaries
                        # -----------------------------

                        nucleotide_summary_url = (
                            "https://eutils.ncbi.nlm.nih.gov/"
                            "entrez/eutils/esummary.fcgi"
                        )

                        nucleotide_summary_params = {
                            "db": "nuccore",
                            "id": ",".join(nucleotide_ids),
                            "retmode": "json"
                        }

                        nucleotide_response = requests.get(
                            nucleotide_summary_url,
                            params=nucleotide_summary_params
                        )

                        if nucleotide_response.status_code == 200:

                            nucleotide_data = (
                                nucleotide_response.json()
                            )

                            nucleotide_results = (
                                nucleotide_data.get(
                                    "result",
                                    {}
                                )
                            )

                            # -----------------------------
                            # Find RefSeq mRNA records
                            # -----------------------------

                            refseq_records = []

                            for nucleotide_id in nucleotide_ids:

                                record = nucleotide_results.get(
                                    nucleotide_id,
                                    {}
                                )

                                accession = record.get(
                                    "accessionversion",
                                    ""
                                )

                                title = record.get(
                                    "title",
                                    ""
                                )

                                if accession.startswith(
                                    "NM_"
                                ):

                                    refseq_records.append(
                                        {
                                            "id": nucleotide_id,
                                            "accession": accession,
                                            "title": title
                                        }
                                    )

                            # -----------------------------
                            # Display RefSeq records
                            # -----------------------------

                            if refseq_records:

                                st.subheader(
                                    "RefSeq mRNA Candidates"
                                )

                                for record in refseq_records:

                                    st.write(
                                        f"**{record['accession']}**"
                                    )

                                    st.write(
                                        record["title"]
                                    )

                                    st.write(
                                        f"NCBI nucleotide ID: "
                                        f"{record['id']}"
                                    )

                            else:

                                st.warning(
                                    "No RefSeq mRNA records "
                                    "were identified."
                                )

                        else:

                            st.error(
                                "Could not retrieve nucleotide "
                                "record summaries."
                            )

                    except (
                        KeyError,
                        IndexError
                    ):

                        st.warning(
                            "No linked nucleotide records "
                            "were found for this gene."
                        )

                else:

                    st.error(
                        "Could not retrieve linked "
                        "nucleotide records."
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

import streamlit as st
import requests
import re

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
                        # Inspect GenBank records
                        # -----------------------------

                        refseq_records = []

                        efetch_url = (
                            "https://eutils.ncbi.nlm.nih.gov/"
                            "entrez/eutils/efetch.fcgi"
                        )

                        for nucleotide_id in nucleotide_ids:

                            efetch_params = {
                                "db": "nuccore",
                                "id": nucleotide_id,
                                "rettype": "gb",
                                "retmode": "text"
                            }

                            nucleotide_response = requests.get(
                                efetch_url,
                                params=efetch_params
                            )

                            if nucleotide_response.status_code != 200:
                                continue

                            genbank_text = nucleotide_response.text

                            # -----------------------------
                            # Find accession
                            # -----------------------------

                            accession_match = re.search(
                                r"ACCESSION\s+(\S+)",
                                genbank_text
                            )

                            if not accession_match:
                                continue

                            accession = (
                                accession_match.group(1)
                            )

                            # -----------------------------
                            # Find version
                            # -----------------------------

                            version_match = re.search(
                                r"VERSION\s+(\S+)",
                                genbank_text
                            )

                            if version_match:
                                accession_version = (
                                    version_match.group(1)
                                )
                            else:
                                accession_version = accession

                            # -----------------------------
                            # Find definition
                            # -----------------------------

                            definition_match = re.search(
                                r"DEFINITION\s+(.+)",
                                genbank_text
                            )

                            if definition_match:
                                title = (
                                    definition_match.group(1)
                                    .strip()
                                )
                            else:
                                title = "Not available"

                            # -----------------------------
                            # Keep RefSeq mRNA records
                            # -----------------------------

                            if accession.startswith("NM_"):

                                refseq_records.append(
                                    {
                                        "id": nucleotide_id,
                                        "accession": accession,
                                        "version": accession_version,
                                        "title": title,
                                        "genbank": genbank_text
                                    }
                                )

                        # -----------------------------
                        # Display RefSeq candidates
                        # -----------------------------

                        if refseq_records:

                            st.subheader(
                                "RefSeq mRNA Candidates"
                            )

                            for record in refseq_records:

                                st.write(
                                    f"**{record['version']}**"
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

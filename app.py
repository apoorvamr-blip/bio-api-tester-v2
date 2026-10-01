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

    # --------------------------------------------------
    # STEP 1: NCBI Gene Search
    # --------------------------------------------------

    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    search_params = {
        "db": "gene",
        "term": f"{locus_id}[Gene Name]",
        "retmode": "json"
    }

    search_response = requests.get(
        search_url,
        params=search_params,
        timeout=20
    )

    if search_response.status_code != 200:
        st.error("Could not search NCBI Gene database.")
        st.stop()

    search_data = search_response.json()

    gene_ids = search_data["esearchresult"]["idlist"]

    if not gene_ids:
        st.error(f"No NCBI Gene record found for {locus_id}.")
        st.stop()

    if len(gene_ids) > 1:
        st.warning("Multiple NCBI Gene records were found.")
        st.write(gene_ids)
        st.stop()

    gene_id = gene_ids[0]

    # --------------------------------------------------
    # STEP 2: NCBI Gene Summary
    # --------------------------------------------------

    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    summary_params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }

    summary_response = requests.get(
        summary_url,
        params=summary_params,
        timeout=20
    )

    if summary_response.status_code != 200:
        st.error("Could not retrieve NCBI Gene information.")
        st.stop()

    summary_data = summary_response.json()

    gene_data = summary_data["result"][gene_id]

    organism = gene_data.get("organism", {})

    gene_name = gene_data.get("name", "Not available")
    description = gene_data.get("description", "Not available")
    organism_name = organism.get("scientificname", "Not available")
    taxonomy_id = organism.get("taxid", "Not available")
    chromosome = gene_data.get("chromosome", "Not available")

    # --------------------------------------------------
    # STEP 3: Display Gene Information
    # --------------------------------------------------

    st.subheader("Gene Information")

    col1, col2 = st.columns(2)

    with col1:
        st.write(f"**Gene name:** {gene_name}")
        st.write(f"**Description:** {description}")
        st.write(f"**Organism:** {organism_name}")

    with col2:
        st.write(f"**Chromosome:** {chromosome}")
        st.write(f"**Taxonomy ID:** `{taxonomy_id}`")
        st.write(f"**NCBI Gene ID:** `{gene_id}`")

    # --------------------------------------------------
    # STEP 4: Find RefSeq mRNA
    # --------------------------------------------------

    st.subheader("RefSeq mRNA Candidates")

    nuccore_search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    )

    nuccore_params = {
        "db": "nuccore",
        "term": (
            f'"{locus_id}"[All Fields] '
            f'AND {taxonomy_id}[TaxID] '
            f'AND srcdb_refseq[PROP] '
            f'AND biomol_mrna[PROP]'
        ),
        "retmode": "json",
        "retmax": 20
    }

    nuccore_response = requests.get(
        nuccore_search_url,
        params=nuccore_params,
        timeout=20
    )

    if nuccore_response.status_code != 200:
        st.error("Could not search NCBI nucleotide database.")
        st.stop()

    nuccore_data = nuccore_response.json()

    nuccore_ids = nuccore_data["esearchresult"]["idlist"]

    if not nuccore_ids:
        st.error("No RefSeq mRNA was found.")
        st.stop()

    st.write(f"Found {len(nuccore_ids)} RefSeq mRNA candidate(s).")

    # --------------------------------------------------
    # STEP 5: Retrieve GenBank record
    # --------------------------------------------------

    selected_id = nuccore_ids[0]

    efetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    )

    efetch_params = {
        "db": "nuccore",
        "id": selected_id,
        "rettype": "gb",
        "retmode": "text"
    }

    with st.spinner("Retrieving GenBank record..."):

        efetch_response = requests.get(
            efetch_url,
            params=efetch_params,
            timeout=30
        )

    if efetch_response.status_code != 200:
        st.error("Could not retrieve GenBank record.")
        st.stop()

    genbank_text = efetch_response.text

    # --------------------------------------------------
    # STEP 6: Display GenBank record
    # --------------------------------------------------

    st.subheader("GenBank Record")

    with st.expander("View GenBank record"):
        st.code(genbank_text)

    # --------------------------------------------------
    # STEP 7: Find CDS feature
    # --------------------------------------------------

    cds_start = None
    cds_end = None

    for line in genbank_text.splitlines():

        stripped_line = line.strip()

        if stripped_line.startswith("CDS"):
            location = stripped_line.split()[1]

            if ".." in location:

                start, end = location.split("..")

                start = start.replace("<", "")
                end = end.replace(">", "")

                if start.isdigit() and end.isdigit():
                    cds_start = int(start)
                    cds_end = int(end)

                    break

    if cds_start is None or cds_end is None:
        st.error("Could not find a CDS feature in the GenBank record.")
        st.stop()

    st.success(
        f"CDS detected: {cds_start}..{cds_end}"
    )

    st.write(f"**CDS start:** {cds_start}")
    st.write(f"**CDS end:** {cds_end}")
    st.write(f"**CDS length:** {cds_end - cds_start + 1} nt")

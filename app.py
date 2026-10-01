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

    # ==================================================
    # STEP 1: NCBI GENE SEARCH
    # ==================================================

    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    )

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

    try:
        search_data = search_response.json()
    except ValueError:
        st.error("NCBI Gene search did not return valid JSON.")
        st.code(search_response.text[:1000])
        st.stop()

    gene_ids = search_data["esearchresult"]["idlist"]

    if not gene_ids:
        st.error(f"No NCBI Gene record found for {locus_id}.")
        st.stop()

    if len(gene_ids) > 1:
        st.warning("Multiple NCBI Gene records were found.")
        st.write(gene_ids)
        st.stop()

    gene_id = gene_ids[0]

    # ==================================================
    # STEP 2: NCBI GENE SUMMARY
    # ==================================================

    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    )

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

    try:
        summary_data = summary_response.json()
    except ValueError:
        st.error("NCBI Gene summary did not return valid JSON.")
        st.code(summary_response.text[:1000])
        st.stop()

    gene_data = summary_data["result"][gene_id]

    organism = gene_data.get("organism", {})

    gene_name = gene_data.get(
        "name",
        "Not available"
    )

    description = gene_data.get(
        "description",
        "Not available"
    )

    organism_name = organism.get(
        "scientificname",
        "Not available"
    )

    taxonomy_id = organism.get(
        "taxid",
        "Not available"
    )

    chromosome = gene_data.get(
        "chromosome",
        "Not available"
    )

    # ==================================================
    # STEP 3: GENE INFORMATION
    # ==================================================

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

    # ==================================================
    # STEP 4: SEARCH FOR REFSEQ mRNA
    # ==================================================

    st.subheader("RefSeq mRNA Candidates")

    nuccore_search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    )

    # --------------------------------------------------
    # First search: locus ID
    # --------------------------------------------------

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

    with st.spinner("Searching NCBI for RefSeq mRNA..."):

        nuccore_response = requests.get(
            nuccore_search_url,
            params=nuccore_params,
            timeout=30
        )

    # --------------------------------------------------
    # Check response
    # --------------------------------------------------

    if nuccore_response.status_code != 200:

        st.error(
            "NCBI nucleotide search returned an HTTP error."
        )

        st.write(
            f"HTTP status: {nuccore_response.status_code}"
        )

        with st.expander("NCBI response"):
            st.code(nuccore_response.text[:2000])

        st.stop()

    # --------------------------------------------------
    # Parse JSON safely
    # --------------------------------------------------

    try:

        nuccore_data = nuccore_response.json()

    except ValueError:

        st.error(
            "NCBI nucleotide search did not return valid JSON."
        )

        with st.expander("NCBI response"):
            st.code(nuccore_response.text[:2000])

        st.stop()

    # --------------------------------------------------
    # Extract nucleotide IDs
    # --------------------------------------------------

    nuccore_result = nuccore_data.get(
        "esearchresult",
        {}
    )

    nuccore_ids = nuccore_result.get(
        "idlist",
        []
    )

    # --------------------------------------------------
    # If locus search finds nothing, try gene name
    # --------------------------------------------------

    if not nuccore_ids:

        st.info(
            "No RefSeq mRNA found using the locus ID. "
            "Trying the NCBI gene name..."
        )

        fallback_params = {
            "db": "nuccore",
            "term": (
                f'"{gene_name}"[All Fields] '
                f'AND {taxonomy_id}[TaxID] '
                f'AND srcdb_refseq[PROP] '
                f'AND biomol_mrna[PROP]'
            ),
            "retmode": "json",
            "retmax": 20
        }

        fallback_response = requests.get(
            nuccore_search_url,
            params=fallback_params,
            timeout=30
        )

        if fallback_response.status_code != 200:

            st.error(
                "The fallback NCBI nucleotide search failed."
            )
            st.stop()

        try:

            fallback_data = fallback_response.json()

        except ValueError:

            st.error(
                "The fallback NCBI search did not return JSON."
            )

            with st.expander("NCBI response"):
                st.code(fallback_response.text[:2000])

            st.stop()

        nuccore_ids = (
            fallback_data
            .get("esearchresult", {})
            .get("idlist", [])
        )

    # ==================================================
    # STEP 5: CHECK CANDIDATES
    # ==================================================

    if not nuccore_ids:

        st.warning(
            "No RefSeq mRNA candidates were found."
        )

        st.stop()

    st.success(
        f"Found {len(nuccore_ids)} RefSeq mRNA candidate(s)."
    )

    # ==================================================
    # STEP 6: RETRIEVE GENBANK RECORD
    # ==================================================

    # For now, use the first candidate.
    # We will improve candidate selection later.

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

        st.error(
            "Could not retrieve the GenBank record."
        )

        st.stop()

    genbank_text = efetch_response.text

    # ==================================================
    # STEP 7: DISPLAY GENBANK RECORD
    # ==================================================

    st.subheader("GenBank Record")

    with st.expander("View GenBank record"):
        st.code(genbank_text)

    # ==================================================
    # STEP 8: FIND CDS
    # ==================================================

    cds_start = None
    cds_end = None

    for line in genbank_text.splitlines():

        stripped_line = line.strip()

        if stripped_line.startswith("CDS"):

            parts = stripped_line.split()

            if len(parts) >= 2:

                location = parts[1]

                if ".." in location:

                    start, end = location.split("..")

                    start = start.replace("<", "")
                    end = end.replace(">", "")

                    if start.isdigit() and end.isdigit():

                        cds_start = int(start)
                        cds_end = int(end)

                        break

    # ==================================================
    # STEP 9: DISPLAY CDS INFORMATION
    # ==================================================

    if cds_start is None or cds_end is None:

        st.error(
            "Could not find a CDS feature in the GenBank record."
        )

        st.stop()

    cds_length = cds_end - cds_start + 1

    st.subheader("Coding Sequence (CDS)")

    st.success(
        f"CDS detected: {cds_start}..{cds_end}"
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("CDS start", cds_start)

    with col2:
        st.metric("CDS end", cds_end)

    with col3:
        st.metric("CDS length", f"{cds_length} nt")

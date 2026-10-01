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
    # STEP 4: LINK GENE TO NUCLEOTIDE RECORDS
    # ==================================================

    st.subheader("Nucleotide Records")

    elink_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    )

    elink_params = {
        "dbfrom": "gene",
        "db": "nuccore",
        "id": gene_id,
        "retmode": "json"
    }

    with st.spinner("Finding linked nucleotide records..."):

        elink_response = requests.get(
            elink_url,
            params=elink_params,
            timeout=30
        )

    if elink_response.status_code != 200:
        st.error("Could not retrieve linked nucleotide records.")
        st.stop()

    try:
        elink_data = elink_response.json()
    except ValueError:
        st.error("NCBI ELink did not return valid JSON.")
        st.stop()

    # Extract nucleotide IDs
    linksets = elink_data.get("linksets", [])

    if not linksets:
        st.error("No nucleotide records were linked to this gene.")
        st.stop()

    linksetdbs = linksets[0].get("linksetdbs", [])

    nuccore_ids = []

    for linkset in linksetdbs:

        if linkset.get("dbto") == "nuccore":

            nuccore_ids.extend(
                linkset.get("links", [])
            )

    if not nuccore_ids:
        st.error("No nucleotide records were found.")
        st.stop()

    st.write(
        f"Found {len(nuccore_ids)} linked nucleotide records."
    )

    # ==================================================
    # STEP 5: RETRIEVE LINKED RECORDS
    # ==================================================

    efetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    )

    efetch_params = {
        "db": "nuccore",
        "id": ",".join(nuccore_ids),
        "rettype": "gb",
        "retmode": "text"
    }

    with st.spinner("Searching linked records for RefSeq mRNA..."):

        efetch_response = requests.get(
            efetch_url,
            params=efetch_params,
            timeout=60
        )

    if efetch_response.status_code != 200:
        st.error("Could not retrieve nucleotide records.")
        st.stop()

    genbank_text = efetch_response.text

    # ==================================================
    # STEP 6: SPLIT GENBANK RECORDS
    # ==================================================

    records = genbank_text.split("\n//")

    refseq_candidates = []

    for record in records:

        if "LOCUS" not in record:
            continue

        # Look for RefSeq mRNA records.
        # RefSeq mRNA accessions normally begin with NM_.

        accession = None
        definition = None

        for line in record.splitlines():

            if line.startswith("ACCESSION"):
                parts = line.split()

                if len(parts) >= 2:
                    accession = parts[1]

            elif line.startswith("DEFINITION"):
                definition = line.replace(
                    "DEFINITION", ""
                ).strip()

        if accession and accession.startswith("NM_"):

            refseq_candidates.append({
                "accession": accession,
                "definition": definition,
                "record": record
            })

    # ==================================================
    # STEP 7: CHECK REFSEQ CANDIDATES
    # ==================================================

    st.subheader("RefSeq mRNA Candidates")

    if not refseq_candidates:

        st.warning(
            "No RefSeq mRNA records were found among the linked "
            "nucleotide records."
        )

        st.stop()

    st.success(
        f"Found {len(refseq_candidates)} RefSeq mRNA candidate(s)."
    )

    # ==================================================
    # STEP 8: DISPLAY CANDIDATES
    # ==================================================

    for candidate in refseq_candidates:

        st.write(
            f"**{candidate['accession']}**"
        )

        if candidate["definition"]:
            st.write(candidate["definition"])

    # ==================================================
    # STEP 9: SELECT RECORD
    # ==================================================

    if len(refseq_candidates) > 1:

        candidate_names = [
            candidate["accession"]
            for candidate in refseq_candidates
        ]

        selected_accession = st.selectbox(
            "Select a RefSeq mRNA:",
            candidate_names
        )

        selected_candidate = next(
            candidate
            for candidate in refseq_candidates
            if candidate["accession"] == selected_accession
        )

    else:

        selected_candidate = refseq_candidates[0]

    selected_accession = selected_candidate["accession"]
    selected_genbank = selected_candidate["record"]

    # ==================================================
    # STEP 10: GENBANK RECORD
    # ==================================================

    st.subheader("GenBank Record")

    st.write(
        f"**Selected RefSeq:** `{selected_accession}`"
    )

    with st.expander("View GenBank record"):
        st.code(selected_genbank)

    # ==================================================
    # STEP 11: FIND CDS
    # ==================================================

    cds_start = None
    cds_end = None

    for line in selected_genbank.splitlines():

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

    if cds_start is None or cds_end is None:

        st.error(
            "Could not find a CDS feature in the selected "
            "GenBank record."
        )

        st.stop()

    # ==================================================
    # STEP 12: CDS INFORMATION
    # ==================================================

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

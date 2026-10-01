import streamlit as st
import requests

# ==================================================
# PAGE CONFIGURATION
# ==================================================

st.set_page_config(
    page_title="Bio API Tester",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Bio API Tester")
st.write("Enter a Locus ID to retrieve biological sequence information.")


# ==================================================
# USER INPUT
# ==================================================

locus_id = st.text_input(
    "Locus ID",
    placeholder="Example: At1g01010"
)


# ==================================================
# MAIN WORKFLOW
# ==================================================

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

    with st.spinner("Searching NCBI Gene database..."):

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

        with st.expander("NCBI response"):
            st.code(search_response.text[:2000])

        st.stop()

    gene_ids = (
        search_data
        .get("esearchresult", {})
        .get("idlist", [])
    )

    # --------------------------------------------------
    # Check gene results
    # --------------------------------------------------

    if not gene_ids:

        st.error(
            f"No NCBI Gene record found for `{locus_id}`."
        )

        st.stop()

    if len(gene_ids) > 1:

        st.warning(
            "Multiple NCBI Gene records were found."
        )

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

    with st.spinner("Retrieving gene information..."):

        summary_response = requests.get(
            summary_url,
            params=summary_params,
            timeout=20
        )

    if summary_response.status_code != 200:

        st.error(
            "Could not retrieve NCBI Gene information."
        )

        st.stop()

    try:

        summary_data = summary_response.json()

    except ValueError:

        st.error(
            "NCBI Gene summary did not return valid JSON."
        )

        with st.expander("NCBI response"):
            st.code(summary_response.text[:2000])

        st.stop()

    gene_data = (
        summary_data
        .get("result", {})
        .get(gene_id)
    )

    if not gene_data:

        st.error("NCBI Gene summary was empty.")

        st.stop()


    # ==================================================
    # EXTRACT GENE INFORMATION
    # ==================================================

    organism = gene_data.get(
        "organism",
        {}
    )

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
    # STEP 3: DISPLAY GENE INFORMATION
    # ==================================================

    st.subheader("Gene Information")

    col1, col2 = st.columns(2)

    with col1:

        st.write(
            f"**Gene name:** {gene_name}"
        )

        st.write(
            f"**Description:** {description}"
        )

        st.write(
            f"**Organism:** {organism_name}"
        )

    with col2:

        st.write(
            f"**Chromosome:** {chromosome}"
        )

        st.write(
            f"**Taxonomy ID:** `{taxonomy_id}`"
        )

        st.write(
            f"**NCBI Gene ID:** `{gene_id}`"
        )


    # ==================================================
    # STEP 4: LINK GENE → NUCLEOTIDE
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

    with st.spinner(
        "Finding linked nucleotide records..."
    ):

        elink_response = requests.get(
            elink_url,
            params=elink_params,
            timeout=30
        )

    if elink_response.status_code != 200:

        st.error(
            "Could not retrieve linked nucleotide records."
        )

        st.stop()

    try:

        elink_data = elink_response.json()

    except ValueError:

        st.error(
            "NCBI ELink did not return valid JSON."
        )

        with st.expander("NCBI response"):
            st.code(elink_response.text[:2000])

        st.stop()


    # ==================================================
    # EXTRACT NUCLEOTIDE IDs
    # ==================================================

    linksets = elink_data.get(
        "linksets",
        []
    )

    if not linksets:

        st.error(
            "No nucleotide records were linked to this gene."
        )

        st.stop()

    linksetdbs = linksets[0].get(
        "linksetdbs",
        []
    )

    nuccore_ids = []

    for linkset in linksetdbs:

        if linkset.get("dbto") == "nuccore":

            nuccore_ids.extend(
                linkset.get("links", [])
            )


    # Remove duplicates while preserving order

    nuccore_ids = list(
        dict.fromkeys(nuccore_ids)
    )


    if not nuccore_ids:

        st.error(
            "No nucleotide records were found."
        )

        st.stop()


    st.write(
        f"Found {len(nuccore_ids)} linked nucleotide records."
    )


    # ==================================================
    # STEP 5: RETRIEVE LINKED GENBANK RECORDS
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

    with st.spinner(
        "Searching linked records for RefSeq mRNA..."
    ):

        efetch_response = requests.get(
            efetch_url,
            params=efetch_params,
            timeout=60
        )

    if efetch_response.status_code != 200:

        st.error(
            "Could not retrieve nucleotide records."
        )

        st.write(
            f"HTTP status: {efetch_response.status_code}"
        )

        with st.expander("NCBI response"):
            st.code(efetch_response.text[:2000])

        st.stop()


    genbank_text = efetch_response.text


    if not genbank_text.strip():

        st.error(
            "NCBI returned an empty GenBank response."
        )

        st.stop()


    # ==================================================
    # STEP 6: SPLIT GENBANK RECORDS
    # ==================================================

    records = genbank_text.split("\n//")

    refseq_candidates = []


    for record in records:

        if "LOCUS" not in record:
            continue


        accession = None
        version = None
        definition = None


        # ----------------------------------------------
        # Find ACCESSION
        # ----------------------------------------------

        for line in record.splitlines():

            if line.startswith("ACCESSION"):

                parts = line.split()

                if len(parts) >= 2:

                    accession = parts[1]


            # ------------------------------------------
            # Find VERSION
            # ------------------------------------------

            elif line.startswith("VERSION"):

                parts = line.split()

                if len(parts) >= 2:

                    version = parts[1]


            # ------------------------------------------
            # Find DEFINITION
            # ------------------------------------------

            elif line.startswith("DEFINITION"):

                definition = line.replace(
                    "DEFINITION",
                    "",
                    1
                ).strip()


        # ----------------------------------------------
        # Keep RefSeq mRNA records
        # ----------------------------------------------

        if accession and accession.startswith("NM_"):

            refseq_candidates.append(
                {
                    "accession": accession,
                    "version": version,
                    "definition": definition,
                    "record": record
                }
            )


    # ==================================================
    # STEP 7: REFSEQ CANDIDATES
    # ==================================================

    st.subheader("RefSeq mRNA Candidates")


    if not refseq_candidates:

        st.warning(
            "No RefSeq mRNA records were found among "
            "the linked nucleotide records."
        )

        st.stop()


    st.success(
        f"Found {len(refseq_candidates)} "
        f"RefSeq mRNA candidate(s)."
    )


    # ==================================================
    # STEP 8: DISPLAY CANDIDATES
    # ==================================================

    for candidate in refseq_candidates:

        accession_display = candidate["accession"]

        if candidate["version"]:

            accession_display = candidate["version"]

        st.write(
            f"**{accession_display}**"
        )

        if candidate["definition"]:

            st.write(
                candidate["definition"]
            )


    # ==================================================
    # STEP 9: SELECT CANDIDATE
    # ==================================================

    if len(refseq_candidates) > 1:

        candidate_names = []

        for candidate in refseq_candidates:

            if candidate["version"]:

                candidate_names.append(
                    candidate["version"]
                )

            else:

                candidate_names.append(
                    candidate["accession"]
                )


        selected_name = st.selectbox(
            "Select a RefSeq mRNA:",
            candidate_names
        )


        selected_candidate = None

        for candidate in refseq_candidates:

            candidate_name = (
                candidate["version"]
                if candidate["version"]
                else candidate["accession"]
            )

            if candidate_name == selected_name:

                selected_candidate = candidate

                break

    else:

        selected_candidate = refseq_candidates[0]


    if selected_candidate is None:

        st.error(
            "Could not determine the selected RefSeq record."
        )

        st.stop()


    selected_accession = selected_candidate[
        "accession"
    ]

    selected_version = selected_candidate[
        "version"
    ]

    selected_genbank = selected_candidate[
        "record"
    ]


    # ==================================================
    # STEP 10: DISPLAY SELECTED RECORD
    # ==================================================

    st.subheader("Selected RefSeq mRNA")

    if selected_version:

        st.write(
            f"**Accession:** `{selected_version}`"
        )

    else:

        st.write(
            f"**Accession:** `{selected_accession}`"
        )


    # Find NCBI nucleotide ID

    selected_nuccore_id = None

    for record_id in nuccore_ids:

        # We will determine the matching record
        # by checking the fetched GenBank record.
        #
        # The selected record itself does not always
        # expose the numeric ID directly in the text,
        # so this remains informational for now.

        pass


    with st.expander("View GenBank record"):

        st.code(
            selected_genbank
        )


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


                # Handle simple locations:
                #
                # 130..1419
                #
                # More complex CDS locations such as
                # join(...) will be handled later.

                if ".." in location:

                    start, end = location.split(
                        "..",
                        1
                    )


                    start = start.replace(
                        "<",
                        ""
                    )

                    end = end.replace(
                        ">",
                        ""
                    )


                    if (
                        start.isdigit()
                        and end.isdigit()
                    ):

                        cds_start = int(start)
                        cds_end = int(end)

                        break


    # ==================================================
    # STEP 12: CHECK CDS
    # ==================================================

    if cds_start is None or cds_end is None:

        st.error(
            "Could not find a CDS feature in the "
            "selected GenBank record."
        )

        st.stop()


    cds_length = (
        cds_end - cds_start + 1
    )


    # ==================================================
    # STEP 13: DISPLAY CDS INFORMATION
    # ==================================================

    st.subheader(
        "Coding Sequence (CDS)"
    )


    st.success(
        f"CDS detected: "
        f"{cds_start}..{cds_end}"
    )


    col1, col2, col3 = st.columns(3)


    with col1:

        st.metric(
            "CDS start",
            cds_start
        )


    with col2:

        st.metric(
            "CDS end",
            cds_end
        )


    with col3:

        st.metric(
            "CDS length",
            f"{cds_length} nt"
        )

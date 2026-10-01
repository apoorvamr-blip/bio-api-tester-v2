import time

import requests
import streamlit as st


# ==================================================
# PAGE CONFIGURATION
# ==================================================

st.set_page_config(
    page_title="Bio API Tester",
    page_icon="🧬",
    layout="wide"
)


# ==================================================
# NCBI CONFIGURATION
# ==================================================

NCBI_BASE_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
)

# Identify your application to NCBI.
# Replace the email with your own email before deployment.
NCBI_EMAIL = "your_email@example.com"

NCBI_TOOL = "BioAPITester"


# ==================================================
# HELPER FUNCTION
# ==================================================

def ncbi_get(endpoint, params, timeout=30):
    """
    Send a request to NCBI E-utilities.

    Returns the response object or raises a
    user-friendly Streamlit error.
    """

    params = params.copy()

    params["tool"] = NCBI_TOOL
    params["email"] = NCBI_EMAIL

    url = f"{NCBI_BASE_URL}/{endpoint}"

    response = requests.get(
        url,
        params=params,
        timeout=timeout
    )

    # --------------------------------------------------
    # Rate limit
    # --------------------------------------------------

    if response.status_code == 429:

        st.error(
            "NCBI rate limit reached. "
            "Please wait a few seconds and try again."
        )

        st.stop()

    # --------------------------------------------------
    # Other HTTP errors
    # --------------------------------------------------

    if response.status_code != 200:

        st.error(
            f"NCBI request failed "
            f"(HTTP {response.status_code})."
        )

        st.stop()

    return response


# ==================================================
# STEP 1: NCBI GENE SEARCH
# ==================================================

@st.cache_data(ttl=3600)
def search_gene(locus_id):

    params = {
        "db": "gene",
        "term": f"{locus_id}[Gene Name]",
        "retmode": "json"
    }

    response = ncbi_get(
        "esearch.fcgi",
        params,
        timeout=20
    )

    try:
        data = response.json()
    except ValueError:

        st.error(
            "NCBI Gene search returned an unexpected response."
        )

        st.stop()

    return data


# ==================================================
# STEP 2: GENE SUMMARY
# ==================================================

@st.cache_data(ttl=3600)
def get_gene_summary(gene_id):

    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }

    response = ncbi_get(
        "esummary.fcgi",
        params,
        timeout=20
    )

    try:
        data = response.json()
    except ValueError:

        st.error(
            "NCBI Gene summary returned an unexpected response."
        )

        st.stop()

    return data


# ==================================================
# STEP 3: GENE → NUCLEOTIDE LINKS
# ==================================================

@st.cache_data(ttl=3600)
def get_nucleotide_links(gene_id):

    params = {
        "dbfrom": "gene",
        "db": "nuccore",
        "id": gene_id,
        "retmode": "json"
    }

    response = ncbi_get(
        "elink.fcgi",
        params,
        timeout=30
    )

    try:
        data = response.json()
    except ValueError:

        st.error(
            "NCBI ELink returned an unexpected response."
        )

        st.stop()

    return data


# ==================================================
# STEP 4: FETCH GENBANK RECORDS
# ==================================================

@st.cache_data(ttl=3600)
def fetch_genbank_records(nuccore_ids):

    params = {
        "db": "nuccore",
        "id": ",".join(nuccore_ids),
        "rettype": "gb",
        "retmode": "text"
    }

    response = ncbi_get(
        "efetch.fcgi",
        params,
        timeout=60
    )

    return response.text


# ==================================================
# MAIN APPLICATION
# ==================================================

st.title("🧬 Bio API Tester")

st.write(
    "Enter a Locus ID to retrieve biological "
    "sequence information."
)


locus_id = st.text_input(
    "Locus ID",
    placeholder="Example: At1g01010"
)


if locus_id:

    # Normalize input

    locus_id = locus_id.strip()


    # ==================================================
    # GENE SEARCH
    # ==================================================

    with st.spinner(
        "Searching NCBI Gene database..."
    ):

        search_data = search_gene(
            locus_id
        )


    gene_ids = (
        search_data
        .get("esearchresult", {})
        .get("idlist", [])
    )


    if not gene_ids:

        st.error(
            f"No NCBI Gene record found for `{locus_id}`."
        )

        st.stop()


    # --------------------------------------------------
    # Ambiguous gene
    # --------------------------------------------------

    if len(gene_ids) > 1:

        st.warning(
            "Multiple NCBI Gene records were found. "
            "We will handle candidate selection in a "
            "later step."
        )

        st.write(gene_ids)

        st.stop()


    gene_id = gene_ids[0]


    # ==================================================
    # GENE SUMMARY
    # ==================================================

    with st.spinner(
        "Retrieving gene information..."
    ):

        summary_data = get_gene_summary(
            gene_id
        )


    gene_data = (
        summary_data
        .get("result", {})
        .get(gene_id)
    )


    if not gene_data:

        st.error(
            "NCBI Gene summary was empty."
        )

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
    # GENE INFORMATION
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
    # GENE → NUCLEOTIDE
    # ==================================================

    st.subheader("Nucleotide Records")


    with st.spinner(
        "Finding linked nucleotide records..."
    ):

        elink_data = get_nucleotide_links(
            gene_id
        )


    linksets = elink_data.get(
        "linksets",
        []
    )


    if not linksets:

        st.error(
            "No nucleotide records were linked "
            "to this gene."
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


    # Remove duplicates

    nuccore_ids = list(
        dict.fromkeys(nuccore_ids)
    )


    if not nuccore_ids:

        st.error(
            "No nucleotide records were found."
        )

        st.stop()


    st.write(
        f"Found {len(nuccore_ids)} "
        f"linked nucleotide records."
    )


    # ==================================================
    # SMALL DELAY BEFORE EFETCH
    # ==================================================

    # Helps avoid immediately hitting NCBI's
    # request-rate limit.

    time.sleep(0.5)


    # ==================================================
    # FETCH GENBANK
    # ==================================================

    with st.spinner(
        "Retrieving linked GenBank records..."
    ):

        genbank_text = fetch_genbank_records(
            nuccore_ids
        )


    if not genbank_text.strip():

        st.error(
            "NCBI returned an empty GenBank response."
        )

        st.stop()


    # ==================================================
    # SPLIT GENBANK RECORDS
    # ==================================================

    records = genbank_text.split(
        "\n//"
    )


    refseq_candidates = []


    for record in records:

        if "LOCUS" not in record:
            continue


        accession = None
        version = None
        definition = None


        for line in record.splitlines():

            # ------------------------------------------
            # ACCESSION
            # ------------------------------------------

            if line.startswith(
                "ACCESSION"
            ):

                parts = line.split()

                if len(parts) >= 2:

                    accession = parts[1]


            # ------------------------------------------
            # VERSION
            # ------------------------------------------

            elif line.startswith(
                "VERSION"
            ):

                parts = line.split()

                if len(parts) >= 2:

                    version = parts[1]


            # ------------------------------------------
            # DEFINITION
            # ------------------------------------------

            elif line.startswith(
                "DEFINITION"
            ):

                definition = line.replace(
                    "DEFINITION",
                    "",
                    1
                ).strip()


        # ----------------------------------------------
        # RefSeq mRNA
        # ----------------------------------------------

        if (
            accession
            and accession.startswith("NM_")
        ):

            refseq_candidates.append(
                {
                    "accession": accession,
                    "version": version,
                    "definition": definition,
                    "record": record
                }
            )


    # ==================================================
    # REFSEQ CANDIDATES
    # ==================================================

    st.subheader(
        "RefSeq mRNA Candidates"
    )


    if not refseq_candidates:

        st.warning(
            "No RefSeq mRNA records were found "
            "among the linked nucleotide records."
        )

        st.stop()


    st.success(
        f"Found {len(refseq_candidates)} "
        f"RefSeq mRNA candidate(s)."
    )


    # ==================================================
    # CANDIDATE SELECTION
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


        selected_candidate = next(
            candidate
            for candidate in refseq_candidates
            if (
                candidate["version"]
                or candidate["accession"]
            ) == selected_name
        )

    else:

        selected_candidate = (
            refseq_candidates[0]
        )


    # ==================================================
    # DISPLAY SELECTED REFSEQ
    # ==================================================

    selected_accession = (
        selected_candidate["accession"]
    )


    selected_version = (
        selected_candidate["version"]
    )


    selected_genbank = (
        selected_candidate["record"]
    )


    st.subheader(
        "Selected RefSeq mRNA"
    )


    if selected_version:

        st.write(
            f"**Accession:** `{selected_version}`"
        )

    else:

        st.write(
            f"**Accession:** `{selected_accession}`"
        )


    if selected_candidate["definition"]:

        st.write(
            selected_candidate["definition"]
        )


    # ==================================================
    # GENBANK RECORD
    # ==================================================

    with st.expander(
        "View GenBank record"
    ):

        st.code(
            selected_genbank
        )


    # ==================================================
    # FIND CDS
    # ==================================================

    cds_start = None
    cds_end = None


    for line in selected_genbank.splitlines():

        stripped_line = line.strip()


        if stripped_line.startswith(
            "CDS"
        ):

            parts = stripped_line.split()


            if len(parts) >= 2:

                location = parts[1]


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
    # CDS RESULT
    # ==================================================

    if (
        cds_start is None
        or cds_end is None
    ):

        st.error(
            "Could not find a CDS feature "
            "in the selected GenBank record."
        )

        st.stop()


    cds_length = (
        cds_end - cds_start + 1
    )


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

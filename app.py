import time
import requests
import streamlit as st


# -----------------------------
# Page configuration
# -----------------------------

st.set_page_config(
    page_title="Bio API Tester",
    page_icon="🧬",
    layout="wide"
)


# -----------------------------
# NCBI configuration
# -----------------------------

NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Replace this with your real email before deployment.
NCBI_EMAIL = "your_email@example.com"

NCBI_TOOL = "BioAPITester"


# -----------------------------
# NCBI request helper
# -----------------------------

def ncbi_get(endpoint, params, timeout=30):

    params = params.copy()

    params["tool"] = NCBI_TOOL
    params["email"] = NCBI_EMAIL

    # Stay below NCBI's public request-rate limit.
    time.sleep(0.4)

    url = f"{NCBI_BASE_URL}/{endpoint}"

    response = requests.get(
        url,
        params=params,
        timeout=timeout
    )

    if response.status_code == 429:

        st.error(
            "NCBI rate limit reached. "
            "Please wait a few seconds and try again."
        )

        st.stop()

    if response.status_code != 200:

        st.error(
            f"NCBI request failed. "
            f"HTTP status: {response.status_code}"
        )

        st.stop()

    return response


# -----------------------------
# NCBI Gene search
# -----------------------------

@st.cache_data(ttl=3600)
def search_gene(locus_id):

    params = {
        "db": "gene",
        "term": f"{locus_id}[Gene Name]",
        "retmode": "json"
    }

    response = ncbi_get(
        "esearch.fcgi",
        params
    )

    data = response.json()

    return data["esearchresult"]["idlist"]


# -----------------------------
# NCBI Gene summary
# -----------------------------

@st.cache_data(ttl=3600)
def get_gene_summary(gene_id):

    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }

    response = ncbi_get(
        "esummary.fcgi",
        params
    )

    data = response.json()

    return data["result"][gene_id]


# -----------------------------
# Gene → nucleotide links
# -----------------------------

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
        params
    )

    data = response.json()

    nucleotide_ids = []

    try:

        linksets = data["linksets"]

        for linkset in linksets:

            for link in linkset.get(
                "linksetdbs",
                []
            ):

                if link.get("dbto") == "nuccore":

                    nucleotide_ids.extend(
                        link.get("links", [])
                    )

    except (KeyError, TypeError):

        pass

    return nucleotide_ids


# -----------------------------
# Lightweight nucleotide summaries
# -----------------------------

@st.cache_data(ttl=3600)
def get_nucleotide_summaries(nuccore_ids):

    params = {
        "db": "nuccore",
        "id": ",".join(nuccore_ids),
        "retmode": "json"
    }

    response = ncbi_get(
        "esummary.fcgi",
        params
    )

    data = response.json()

    result = data.get(
        "result",
        {}
    )

    summaries = []

    for uid in result.get(
        "uids",
        []
    ):

        item = result.get(
            uid,
            {}
        )

        accession = (
            item.get("accessionversion")
            or item.get("caption")
            or ""
        )

        title = item.get(
            "title",
            ""
        )

        length = item.get(
            "slen"
        )

        summaries.append(
            {
                "uid": uid,
                "accession": accession,
                "title": title,
                "length": length
            }
        )

    return summaries


# -----------------------------
# Fetch ONE GenBank record
# -----------------------------

@st.cache_data(ttl=3600)
def fetch_one_genbank(nuccore_id):

    params = {
        "db": "nuccore",
        "id": nuccore_id,
        "rettype": "gb",
        "retmode": "text"
    }

    response = ncbi_get(
        "efetch.fcgi",
        params
    )

    return response.text


# -----------------------------
# Extract CDS from GenBank
# -----------------------------

def extract_cds(genbank_record):

    lines = genbank_record.splitlines()

    cds_start = None
    cds_end = None

    # Find CDS coordinates
    for line in lines:

        if line.startswith("     CDS"):

            location = line[21:].strip()

            if ".." in location:

                start, end = location.split("..")

                cds_start = int(start)
                cds_end = int(end)

            break

    if cds_start is None or cds_end is None:

        return None, None, None

    # Find ORIGIN
    origin_index = None

    for i, line in enumerate(lines):

        if line.startswith("ORIGIN"):

            origin_index = i

            break

    if origin_index is None:

        return None, None, None

    # Collect nucleotide sequence
    sequence_parts = []

    for line in lines[origin_index + 1:]:

        if line.startswith("//"):

            break

        parts = line.split()

        # Remove the position number
        if parts:

            sequence_parts.extend(
                parts[1:]
            )

    nucleotide_sequence = "".join(
        sequence_parts
    ).upper()

    # Extract CDS
    cds_sequence = nucleotide_sequence[
        cds_start - 1:cds_end
    ]

    return (
        cds_sequence,
        cds_start,
        cds_end
    )


# -----------------------------
# Header
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
# Main workflow
# -----------------------------

if locus_id:

    locus_id = locus_id.strip()

    if not locus_id:

        st.warning(
            "Please enter a Locus ID."
        )

        st.stop()

    # -----------------------------
    # Step 1: Find gene
    # -----------------------------

    st.subheader("1️⃣ Gene Search")

    with st.spinner(
        "Searching NCBI Gene..."
    ):

        gene_ids = search_gene(
            locus_id
        )

    if not gene_ids:

        st.error(
            f"No NCBI Gene record found "
            f"for {locus_id}."
        )

        st.stop()

    if len(gene_ids) > 1:

        st.warning(
            "Multiple NCBI Gene records were found. "
            "Please select one."
        )

        selected_gene_id = st.selectbox(
            "Select Gene ID",
            gene_ids
        )

    else:

        selected_gene_id = gene_ids[0]

    # -----------------------------
    # Step 2: Gene information
    # -----------------------------

    with st.spinner(
        "Retrieving gene information..."
    ):

        gene_data = get_gene_summary(
            selected_gene_id
        )

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

    st.subheader(
        "Gene Information"
    )

    col1, col2, col3 = st.columns(3)

    with col1:

        st.metric(
            "Gene",
            gene_name
        )

    with col2:

        st.metric(
            "Organism",
            organism_name
        )

    with col3:

        st.metric(
            "Taxonomy ID",
            taxonomy_id
        )

    st.write(
        f"**Description:** {description}"
    )

    st.write(
        f"**Chromosome:** {chromosome}"
    )

    st.write(
        f"**NCBI Gene ID:** {selected_gene_id}"
    )

    # -----------------------------
    # Step 3: Linked nucleotide records
    # -----------------------------

    st.subheader(
        "2️⃣ Nucleotide Records"
    )

    with st.spinner(
        "Finding linked nucleotide records..."
    ):

        nuccore_ids = get_nucleotide_links(
            selected_gene_id
        )

    if not nuccore_ids:

        st.warning(
            "No nucleotide records were linked "
            "to this Gene record."
        )

        st.stop()

    st.write(
        f"Found {len(nuccore_ids)} "
        f"linked nucleotide records."
    )

    # -----------------------------
    # Step 4: Lightweight summaries
    # -----------------------------

    with st.spinner(
        "Identifying the RefSeq transcript..."
    ):

        summaries = get_nucleotide_summaries(
            tuple(nuccore_ids)
        )

    # Find NM_ RefSeq mRNA records
    refseq_candidates = [
        record
        for record in summaries
        if record["accession"].startswith("NM_")
    ]

    if not refseq_candidates:

        st.warning(
            "No RefSeq NM_ mRNA record "
            "was identified."
        )

        st.write(
            "Available nucleotide records:"
        )

        for record in summaries:

            st.write(
                f"- {record['accession']} — "
                f"{record['title']}"
            )

        st.stop()

    # -----------------------------
    # Step 5: Select transcript
    # -----------------------------

    if len(refseq_candidates) == 1:

        selected_record = refseq_candidates[0]

    else:

        st.warning(
            "Multiple RefSeq mRNA records "
            "were found."
        )

        options = [
            record["accession"]
            for record in refseq_candidates
        ]

        selected_accession = st.selectbox(
            "Select RefSeq transcript",
            options
        )

        selected_record = next(
            record
            for record in refseq_candidates
            if record["accession"]
            == selected_accession
        )

    st.write(
        f"**Selected transcript:** "
        f"{selected_record['accession']}"
    )

    st.write(
        f"**Transcript length:** "
        f"{selected_record['length']} nt"
    )

    st.write(
        f"**NCBI nucleotide ID:** "
        f"{selected_record['uid']}"
    )

    # -----------------------------
    # Step 6: Fetch selected GenBank record
    # -----------------------------

    with st.spinner(
        f"Retrieving "
        f"{selected_record['accession']}..."
    ):

        genbank_record = fetch_one_genbank(
            selected_record["uid"]
        )

    if not genbank_record.strip():

        st.error(
            "The selected GenBank record "
            "was empty."
        )

        st.stop()

    st.success(
        f"Successfully retrieved "
        f"{selected_record['accession']}."
    )

    # -----------------------------
    # GenBank record
    # -----------------------------

    with st.expander(
        "View GenBank record"
    ):

        st.code(
            genbank_record,
            language="text"
        )

    # -----------------------------
    # Step 7: Extract CDS
    # -----------------------------

    st.subheader(
        "3️⃣ Coding Sequence"
    )

    cds_sequence, cds_start, cds_end = extract_cds(
        genbank_record
    )

    if cds_sequence is None:

        st.error(
            "Could not extract the CDS "
            "from the GenBank record."
        )

        st.stop()

    st.success(
        "CDS extracted successfully."
    )

    st.write(
        f"**CDS coordinates:** "
        f"{cds_start}..{cds_end}"
    )

    st.write(
        f"**CDS length:** "
        f"{len(cds_sequence)} nt"
    )

    # -----------------------------
    # Nucleotide FASTA
    # -----------------------------

    nucleotide_fasta = (
        f">{locus_id}_CDS\n"
        f"{cds_sequence}"
    )

    st.write(
        "**Nucleotide FASTA**"
    )

    st.code(
        nucleotide_fasta,
        language="text"
    )

import time
import requests
import pandas as pd
import streamlit as st
from Bio.Seq import Seq


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="Bio API Tester",
    page_icon="🧬",
    layout="wide"
)


# ============================================================
# API CONFIGURATION
# ============================================================

NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Replace this with your real email address
NCBI_EMAIL = "your_email@example.com"

NCBI_TOOL = "BioAPITester"

UNIPROT_SEARCH_URL = (
    "https://rest.uniprot.org/uniprotkb/search"
)

INTERPRO_BASE_URL = (
    "https://www.ebi.ac.uk/interpro/api"
)


# ============================================================
# NCBI REQUEST HELPER
# ============================================================

def ncbi_get(endpoint, params, timeout=30):

    params = params.copy()

    params["tool"] = NCBI_TOOL
    params["email"] = NCBI_EMAIL

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


# ============================================================
# NCBI GENE SEARCH
# ============================================================

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


# ============================================================
# NCBI GENE SUMMARY
# ============================================================

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


# ============================================================
# GENE → NUCLEOTIDE LINKS
# ============================================================

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


# ============================================================
# NUCLEOTIDE SUMMARIES
# ============================================================

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


# ============================================================
# FETCH ONE GENBANK RECORD
# ============================================================

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


# ============================================================
# EXTRACT CDS FROM GENBANK
# ============================================================

def extract_cds(genbank_record):

    lines = genbank_record.splitlines()

    cds_start = None
    cds_end = None

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

    origin_index = None

    for i, line in enumerate(lines):

        if line.startswith("ORIGIN"):

            origin_index = i

            break

    if origin_index is None:

        return None, None, None

    sequence_parts = []

    for line in lines[origin_index + 1:]:

        if line.startswith("//"):

            break

        parts = line.split()

        if parts:

            sequence_parts.extend(
                parts[1:]
            )

    nucleotide_sequence = "".join(
        sequence_parts
    ).upper()

    cds_sequence = nucleotide_sequence[
        cds_start - 1:cds_end
    ]

    return (
        cds_sequence,
        cds_start,
        cds_end
    )


# ============================================================
# TRANSLATE CDS
# ============================================================

def translate_cds(cds_sequence):

    protein_sequence = str(
        Seq(cds_sequence).translate()
    )

    if protein_sequence.endswith("*"):

        protein_sequence = protein_sequence[:-1]

    return protein_sequence


# ============================================================
# UNIPROT SEARCH
# ============================================================

@st.cache_data(ttl=3600)
def search_uniprot(gene_name, taxonomy_id):

    query = (
        f"gene:{gene_name} "
        f"AND organism_id:{taxonomy_id}"
    )

    params = {
        "query": query,
        "format": "json",
        "size": 10
    }

    response = requests.get(
        UNIPROT_SEARCH_URL,
        params=params,
        timeout=30
    )

    if response.status_code != 200:

        return None, response.status_code

    data = response.json()

    return (
        data.get("results", []),
        response.status_code
    )


# ============================================================
# GET UNIPROT RECORD
# ============================================================

@st.cache_data(ttl=3600)
def get_uniprot_record(uniprot_id):

    url = (
        "https://rest.uniprot.org/"
        f"uniprotkb/{uniprot_id}.json"
    )

    response = requests.get(
        url,
        timeout=30
    )

    if response.status_code != 200:

        return None, response.status_code

    return (
        response.json(),
        response.status_code
    )


# ============================================================
# INTERPRO
# ============================================================

@st.cache_data(ttl=3600)
def get_interpro_annotations(uniprot_id):

    url = (
        f"{INTERPRO_BASE_URL}/"
        "entry/interpro/protein/uniprot/"
        f"{uniprot_id}"
    )

    response = requests.get(
        url,
        timeout=30
    )

    if response.status_code != 200:

        return None, response.status_code

    return (
        response.json(),
        response.status_code
    )


# ============================================================
# HEADER
# ============================================================

st.title("🧬 Bio API Tester")

st.write(
    "Retrieve, analyze, and annotate biological "
    "sequences using NCBI, UniProt, and InterPro."
)

st.markdown(
    """
    **Workflow**

    `Locus ID`
    →
    `NCBI Gene`
    →
    `CDS`
    →
    `Protein`
    →
    `UniProt`
    →
    `InterPro`
    →
    `Results`
    """
)

st.divider()


# ============================================================
# USER INPUT
# ============================================================

locus_id = st.text_input(
    "Locus ID",
    placeholder="Example: At1g01010"
)


# ============================================================
# MAIN WORKFLOW
# ============================================================

if locus_id:

    locus_id = locus_id.strip()

    if not locus_id:

        st.warning(
            "Please enter a Locus ID."
        )

        st.stop()


    # ========================================================
    # STEP 1 — GENE SEARCH
    # ========================================================

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
            "Multiple NCBI Gene records were found."
        )

        selected_gene_id = st.selectbox(
            "Select Gene ID",
            gene_ids
        )

    else:

        selected_gene_id = gene_ids[0]


    # ========================================================
    # STEP 2 — GENE INFORMATION
    # ========================================================

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


    # ========================================================
    # STEP 3 — NUCLEOTIDE RECORDS
    # ========================================================

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


    # ========================================================
    # STEP 4 — NUCLEOTIDE SUMMARIES
    # ========================================================

    with st.spinner(
        "Identifying the RefSeq transcript..."
    ):

        summaries = get_nucleotide_summaries(
            tuple(nuccore_ids)
        )

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

        st.stop()


    # ========================================================
    # STEP 5 — SELECT TRANSCRIPT
    # ========================================================

    if len(refseq_candidates) == 1:

        selected_record = refseq_candidates[0]

    else:

        st.warning(
            "Multiple RefSeq mRNA records were found."
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


    # ========================================================
    # STEP 6 — FETCH GENBANK
    # ========================================================

    with st.spinner(
        f"Retrieving "
        f"{selected_record['accession']}..."
    ):

        genbank_record = fetch_one_genbank(
            selected_record["uid"]
        )

    if not genbank_record.strip():

        st.error(
            "The selected GenBank record was empty."
        )

        st.stop()

    st.success(
        f"Successfully retrieved "
        f"{selected_record['accession']}."
    )

    with st.expander(
        "View GenBank record"
    ):

        st.code(
            genbank_record,
            language="text"
        )


    # ========================================================
    # STEP 7 — CODING SEQUENCE
    # ========================================================

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


    # ========================================================
    # NUCLEOTIDE FASTA
    # ========================================================

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


    # ========================================================
    # STEP 8 — PROTEIN TRANSLATION
    # ========================================================

    st.subheader(
        "4️⃣ Protein Sequence"
    )

    if len(cds_sequence) % 3 != 0:

        st.error(
            "CDS length is not divisible by 3. "
            "Translation cannot be performed safely."
        )

        st.stop()

    protein_sequence = translate_cds(
        cds_sequence
    )

    st.success(
        "CDS translated successfully."
    )

    st.write(
        f"**Protein length:** "
        f"{len(protein_sequence)} aa"
    )


    # ========================================================
    # PROTEIN FASTA
    # ========================================================

    protein_fasta = (
        f">{locus_id}_protein\n"
        f"{protein_sequence}"
    )

    st.write(
        "**Protein FASTA**"
    )

    st.code(
        protein_fasta,
        language="text"
    )


    # ========================================================
    # STEP 9 — UNIPROT
    # ========================================================

    st.subheader(
        "5️⃣ UniProt"
    )

    with st.spinner(
        "Searching UniProt..."
    ):

        uniprot_results, uniprot_status = search_uniprot(
            gene_name,
            taxonomy_id
        )

    if uniprot_results is None:

        st.error(
            f"UniProt search failed. "
            f"HTTP status: {uniprot_status}"
        )

        st.stop()

    if not uniprot_results:

        st.warning(
            "No UniProt record was found "
            "for this gene and organism."
        )

        st.stop()


    # ========================================================
    # SELECT UNIPROT RESULT
    # ========================================================

    if len(uniprot_results) == 1:

        selected_uniprot = uniprot_results[0]

    else:

        st.warning(
            f"Found {len(uniprot_results)} "
            "UniProt results."
        )

        uniprot_options = [
            result.get(
                "primaryAccession",
                "Unknown"
            )
            for result in uniprot_results
        ]

        selected_uniprot_id = st.selectbox(
            "Select UniProt entry",
            uniprot_options
        )

        selected_uniprot = next(
            result
            for result in uniprot_results
            if result.get("primaryAccession")
            == selected_uniprot_id
        )


    uniprot_id = selected_uniprot.get(
        "primaryAccession",
        "Not available"
    )

    uniprot_entry_name = selected_uniprot.get(
        "uniProtkbId",
        "Not available"
    )


    # ========================================================
    # GET FULL UNIPROT RECORD
    # ========================================================

    with st.spinner(
        f"Retrieving UniProt record "
        f"{uniprot_id}..."
    ):

        uniprot_record, record_status = get_uniprot_record(
            uniprot_id
        )

    if uniprot_record is None:

        st.error(
            f"Could not retrieve UniProt record. "
            f"HTTP status: {record_status}"
        )

        st.stop()


    uniprot_sequence = (
        uniprot_record
        .get("sequence", {})
        .get("value", "")
    )


    # ========================================================
    # UNIPROT INFORMATION
    # ========================================================

    col1, col2, col3 = st.columns(3)

    with col1:

        st.metric(
            "UniProt ID",
            uniprot_id
        )

    with col2:

        st.metric(
            "Entry Name",
            uniprot_entry_name
        )

    with col3:

        st.metric(
            "Protein Length",
            f"{len(uniprot_sequence)} aa"
        )


    # ========================================================
    # SEQUENCE COMPARISON
    # ========================================================

    sequences_match = (
        protein_sequence
        == uniprot_sequence
    )

    st.write(
        f"**Our translated protein:** "
        f"{len(protein_sequence)} aa"
    )

    st.write(
        f"**UniProt protein:** "
        f"{len(uniprot_sequence)} aa"
    )

    if sequences_match:

        st.success(
            "Protein sequences are identical."
        )

    else:

        st.warning(
            "Protein sequences are not identical."
        )

        st.write(
            "The UniProt sequence and translated "
            "CDS should be reviewed before treating "
            "them as the same protein."
        )


    # ========================================================
    # UNIPROT SEQUENCE
    # ========================================================

    with st.expander(
        "View UniProt protein sequence"
    ):

        st.code(
            uniprot_sequence,
            language="text"
        )


    # ========================================================
    # ANALYSIS SUMMARY TABLE
    # ========================================================

    st.divider()

    st.subheader(
        "📊 Analysis Summary"
    )

    summary_table = pd.DataFrame(
        [
            {
                "Gene": gene_name,
                "Organism": organism_name,
                "NCBI Gene ID": selected_gene_id,
                "RefSeq": selected_record["accession"],
                "CDS Length": f"{len(cds_sequence)} nt",
                "Protein Length": f"{len(protein_sequence)} aa",
                "UniProt": uniprot_id
            }
        ]
    )

    st.dataframe(
        summary_table,
        use_container_width=True,
        hide_index=True
    )


    # ========================================================
    # STEP 10 — INTERPRO
    # ========================================================

    st.subheader(
        "6️⃣ InterPro"
    )

    with st.spinner(
        "Retrieving InterPro annotations..."
    ):

        interpro_data, interpro_status = (
            get_interpro_annotations(
                uniprot_id
            )
        )

    if interpro_data is None:

        st.error(
            f"InterPro request failed. "
            f"HTTP status: {interpro_status}"
        )

        interpro_results = []

    else:

        interpro_results = interpro_data.get(
            "results",
            []
        )

        if not interpro_results:

            st.warning(
                "No InterPro annotations were found "
                "for this UniProt protein."
            )

        else:

            st.success(
                f"Found {len(interpro_results)} "
                "InterPro annotation(s)."
            )

            for result in interpro_results:

                metadata = result.get(
                    "metadata",
                    {}
                )

                accession = metadata.get(
                    "accession",
                    "Not available"
                )

                name = metadata.get(
                    "name",
                    "Not available"
                )

                entry_type = metadata.get(
                    "type",
                    "Not available"
                )

                st.write(
                    f"**InterPro ID:** {accession}"
                )

                st.write(
                    f"**Name:** {name}"
                )

                st.write(
                    f"**Type:** {entry_type}"
                )

                st.divider()


    # ========================================================
    # STEP 11 — FINAL ANNOTATION TABLE
    # ========================================================

    st.subheader(
        "7️⃣ Final Annotation Table"
    )

    annotation_rows = []


    # NCBI
    annotation_rows.append(
        {
            "Source": "NCBI",
            "Record ID": selected_gene_id,
            "Annotation": description,
            "Organism": organism_name,
            "Details": "Gene"
        }
    )


    # RefSeq
    annotation_rows.append(
        {
            "Source": "RefSeq",
            "Record ID": selected_record["accession"],
            "Annotation": selected_record["title"],
            "Organism": organism_name,
            "Details": (
                f"CDS: {cds_start}..{cds_end}; "
                f"Length: {len(cds_sequence)} nt"
            )
        }
    )


    # UniProt
    annotation_rows.append(
        {
            "Source": "UniProt",
            "Record ID": uniprot_id,
            "Annotation": uniprot_entry_name,
            "Organism": organism_name,
            "Details": (
                f"Protein length: "
                f"{len(uniprot_sequence)} aa; "
                f"Sequence match: "
                f"{'Yes' if sequences_match else 'No'}"
            )
        }
    )


    # InterPro
    for result in interpro_results:

        metadata = result.get(
            "metadata",
            {}
        )

        accession = metadata.get(
            "accession",
            "Not available"
        )

        name = metadata.get(
            "name",
            "Not available"
        )

        entry_type = metadata.get(
            "type",
            "Not available"
        )

        annotation_rows.append(
            {
                "Source": "InterPro",
                "Record ID": accession,
                "Annotation": name,
                "Organism": organism_name,
                "Details": entry_type
            }
        )


    annotation_table = pd.DataFrame(
        annotation_rows
    )


    st.dataframe(
        annotation_table,
        use_container_width=True,
        hide_index=True
    )


    # ========================================================
    # STEP 12 — DOWNLOAD RESULTS
    # ========================================================

    st.subheader(
        "8️⃣ Download Results"
    )


    col1, col2 = st.columns(2)

    with col1:

        st.download_button(
            label="⬇️ Download Nucleotide FASTA",
            data=nucleotide_fasta,
            file_name=f"{locus_id}_CDS.fasta",
            mime="text/plain",
            use_container_width=True
        )

    with col2:

        st.download_button(
            label="⬇️ Download Protein FASTA",
            data=protein_fasta,
            file_name=f"{locus_id}_protein.fasta",
            mime="text/plain",
            use_container_width=True
        )


    annotation_csv = annotation_table.to_csv(
        index=False
    )

    annotation_json = annotation_table.to_json(
        orient="records",
        indent=2
    )


    col3, col4 = st.columns(2)

    with col3:

        st.download_button(
            label="⬇️ Download Annotation CSV",
            data=annotation_csv,
            file_name="annotation_results.csv",
            mime="text/csv",
            use_container_width=True
        )

    with col4:

        st.download_button(
            label="⬇️ Download Annotation JSON",
            data=annotation_json,
            file_name="annotation_results.json",
            mime="application/json",
            use_container_width=True
        )

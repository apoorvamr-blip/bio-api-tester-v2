import time
import re
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

# IMPORTANT:
# Replace this with your real email before deployment.
NCBI_EMAIL = "your_email@example.com"

NCBI_TOOL = "BioAPITester"

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

INTERPRO_BASE_URL = "https://www.ebi.ac.uk/interpro/api"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def ncbi_get(endpoint, params, timeout=30):
    """
    Make a request to the NCBI E-utilities API.
    """

    params = params.copy()

    params["tool"] = NCBI_TOOL
    params["email"] = NCBI_EMAIL

    # Small delay to avoid hitting NCBI too quickly
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
            f"NCBI request failed. HTTP status: "
            f"{response.status_code}"
        )
        st.stop()

    return response


def format_fasta(sequence, width=60):
    """
    Format a sequence into FASTA-style lines.
    """

    return "\n".join(
        sequence[i:i + width]
        for i in range(0, len(sequence), width)
    )


# ============================================================
# NCBI FUNCTIONS
# ============================================================

@st.cache_data(ttl=3600)
def search_gene(locus_id):
    """
    Search NCBI Gene using the supplied locus ID.
    """

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

    return data.get("esearchresult", {}).get(
        "idlist",
        []
    )


@st.cache_data(ttl=3600)
def get_gene_summary(gene_id):
    """
    Retrieve NCBI Gene summary.
    """

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

    return data.get("result", {}).get(
        str(gene_id),
        {}
    )


@st.cache_data(ttl=3600)
def get_nucleotide_links(gene_id):
    """
    Find nucleotide records linked to the NCBI Gene record.
    """

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

    linked_ids = []

    try:
        linksets = data.get("linksets", [])

        for linkset in linksets:
            for link in linkset.get("linksetdbs", []):
                linked_ids.extend(
                    link.get("links", [])
                )

    except Exception:
        pass

    return list(dict.fromkeys(linked_ids))


@st.cache_data(ttl=3600)
def get_nucleotide_summaries(nuccore_ids):
    """
    Retrieve lightweight summaries for linked nucleotide records.
    """

    if not nuccore_ids:
        return []

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

    result = data.get("result", {})

    summaries = []

    for uid in result.get("uids", []):

        item = result.get(uid, {})

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

        summaries.append({
            "uid": uid,
            "accession": accession,
            "title": title,
            "length": length
        })

    return summaries


@st.cache_data(ttl=3600)
def fetch_one_genbank(nuccore_id):
    """
    Fetch one nucleotide record in GenBank format.
    """

    params = {
        "db": "nuccore",
        "id": nuccore_id,
        "rettype": "gb",
        "retmode": "text"
    }

    response = ncbi_get(
        "efetch.fcgi",
        params,
        timeout=60
    )

    return response.text


# ============================================================
# GENBANK / CDS FUNCTIONS
# ============================================================

def extract_cds(genbank_record):
    """
    Extract CDS coordinates and nucleotide sequence
    from a GenBank record.
    """

    lines = genbank_record.splitlines()

    cds_start = None
    cds_end = None

    # --------------------------------------------------------
    # Find CDS coordinates
    # --------------------------------------------------------

    for line in lines:

        if line.startswith("     CDS"):

            location = line[21:].strip()

            # Handle normal locations such as:
            # 130..1419
            if ".." in location:

                numbers = re.findall(
                    r"\d+",
                    location
                )

                if len(numbers) >= 2:

                    cds_start = int(numbers[0])
                    cds_end = int(numbers[-1])

            break

    if cds_start is None or cds_end is None:
        return None, None, None

    # --------------------------------------------------------
    # Find ORIGIN
    # --------------------------------------------------------

    origin_index = None

    for i, line in enumerate(lines):

        if line.startswith("ORIGIN"):

            origin_index = i
            break

    if origin_index is None:
        return None, None, None

    # --------------------------------------------------------
    # Extract nucleotide sequence
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Extract CDS
    # --------------------------------------------------------

    cds_sequence = nucleotide_sequence[
        cds_start - 1:cds_end
    ]

    return (
        cds_sequence,
        cds_start,
        cds_end
    )


def translate_cds(cds_sequence):
    """
    Translate nucleotide CDS into protein.
    """

    protein_sequence = str(
        Seq(cds_sequence).translate()
    )

    # Remove terminal stop symbol
    if protein_sequence.endswith("*"):
        protein_sequence = protein_sequence[:-1]

    return protein_sequence


# ============================================================
# UNIPROT FUNCTIONS
# ============================================================

@st.cache_data(ttl=3600)
def search_uniprot(locus_tag, taxonomy_id):
    """
    Search UniProt using locus tag and taxonomy.
    """

    params = {
        "query": (
            f"gene:{locus_tag} "
            f"AND organism_id:{taxonomy_id}"
        ),
        "format": "json",
        "size": 5
    }

    response = requests.get(
        UNIPROT_SEARCH_URL,
        params=params,
        timeout=30
    )

    if response.status_code != 200:
        return None

    return response.json()


@st.cache_data(ttl=3600)
def get_uniprot_record(uniprot_id):
    """
    Retrieve a complete UniProt record.
    """

    url = (
        f"https://rest.uniprot.org/"
        f"uniprotkb/{uniprot_id}.json"
    )

    response = requests.get(
        url,
        timeout=30
    )

    if response.status_code != 200:
        return None

    return response.json()


# ============================================================
# INTERPRO FUNCTIONS
# ============================================================

@st.cache_data(ttl=3600)
def get_interpro_annotations(uniprot_id):
    """
    Retrieve InterPro annotations for a UniProt protein.
    """

    url = (
        f"{INTERPRO_BASE_URL}/entry/interpro/"
        f"protein/uniprot/{uniprot_id}"
    )

    response = requests.get(
        url,
        timeout=60
    )

    if response.status_code != 200:
        return None

    return response.json()


def extract_interpro_details(result):
    """
    Extract useful information from an InterPro result.
    """

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

    # --------------------------------------------------------
    # Member database information
    # --------------------------------------------------------

    member_databases = []

    entries = result.get(
        "entries",
        []
    )

    if isinstance(entries, list):

        for entry in entries:

            if not isinstance(entry, dict):
                continue

            member_accession = (
                entry.get("accession")
                or entry.get("member_accession")
            )

            member_name = (
                entry.get("name")
                or entry.get("member_name")
            )

            if member_accession:

                if member_name:
                    member_databases.append(
                        f"{member_accession} ({member_name})"
                    )
                else:
                    member_databases.append(
                        member_accession
                    )

    # Some InterPro responses use a different structure.
    # Keep the result robust if entries are unavailable.

    if not member_databases:

        members = result.get(
            "member_databases",
            []
        )

        if isinstance(members, list):

            for member in members:

                if isinstance(member, dict):

                    accession_value = member.get(
                        "accession"
                    )

                    name_value = member.get(
                        "name"
                    )

                    if accession_value:

                        if name_value:
                            member_databases.append(
                                f"{accession_value} "
                                f"({name_value})"
                            )
                        else:
                            member_databases.append(
                                accession_value
                            )

    # --------------------------------------------------------
    # Protein / region information
    # --------------------------------------------------------

    matched_regions = []

    proteins = result.get(
        "proteins",
        []
    )

    if isinstance(proteins, list):

        for protein in proteins:

            if not isinstance(protein, dict):
                continue

            entry_protein = protein.get(
                "protein",
                {}
            )

            if not isinstance(entry_protein, dict):
                entry_protein = {}

            protein_id = (
                entry_protein.get(
                    "accession"
                )
                or protein.get(
                    "accession"
                )
            )

            locations = protein.get(
                "entry_protein_locations",
                []
            )

            if isinstance(locations, list):

                for location in locations:

                    if not isinstance(
                        location,
                        dict
                    ):
                        continue

                    start = location.get(
                        "fragments",
                        []
                    )

                    if start:
                        matched_regions.append(
                            str(start)
                        )

            if protein_id:
                pass

    # --------------------------------------------------------
    # GO terms
    # --------------------------------------------------------

    go_terms = []

    if isinstance(
        result.get("go_terms"),
        list
    ):

        for go in result["go_terms"]:

            if isinstance(go, dict):

                go_accession = go.get(
                    "accession"
                )

                go_name = go.get(
                    "name"
                )

                if go_accession:

                    if go_name:
                        go_terms.append(
                            f"{go_accession} ({go_name})"
                        )
                    else:
                        go_terms.append(
                            go_accession
                        )

    return {
        "accession": accession,
        "name": name,
        "type": entry_type,
        "member_databases": member_databases,
        "matched_regions": matched_regions,
        "go_terms": go_terms
    }


# ============================================================
# HEADER
# ============================================================

st.title("🧬 Bio API Tester")

st.write(
    "Retrieve, analyze, and annotate biological "
    "sequences using NCBI, UniProt, and InterPro."
)

st.markdown("### Workflow")

st.markdown(
    """
`Locus ID` → `NCBI Gene` → `CDS` → `Protein`
→ `UniProt` → `InterPro` → `Results`
"""
)

st.divider()


# ============================================================
# LOCUS INPUT
# ============================================================

locus_id = st.text_input(
    "Locus ID",
    placeholder="Example: At1g01010"
)

if not locus_id:
    st.info(
        "Enter a Locus ID to begin the analysis."
    )
    st.stop()


locus_id = locus_id.strip()

# Normalized locus tag for downstream searches
locus_tag = locus_id.upper()


# ============================================================
# STEP 1 — GENE SEARCH
# ============================================================

st.header("1️⃣ Gene Search")

gene_ids = search_gene(
    locus_id
)

if not gene_ids:

    st.error(
        f"No NCBI Gene record found for "
        f"'{locus_id}'."
    )

    st.stop()


# ------------------------------------------------------------
# Handle multiple Gene IDs
# ------------------------------------------------------------

if len(gene_ids) > 1:

    st.warning(
        "Multiple NCBI Gene records were found. "
        "Please select the appropriate record."
    )

    selected_gene_id = st.selectbox(
        "Select NCBI Gene ID",
        gene_ids
    )

else:

    selected_gene_id = gene_ids[0]


gene_data = get_gene_summary(
    selected_gene_id
)

if not gene_data:

    st.error(
        "Could not retrieve the NCBI Gene summary."
    )

    st.stop()


# ============================================================
# GENE INFORMATION
# ============================================================

gene_name = gene_data.get(
    "name",
    "Not available"
)

description = gene_data.get(
    "description",
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

chromosome = gene_data.get(
    "chromosome",
    "Not available"
)


st.subheader("Gene Information")

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

st.markdown(
    f"[🔗 Open NCBI Gene record](https://www.ncbi.nlm.nih.gov/gene/{selected_gene_id})"
)


# ============================================================
# STEP 2 — NUCLEOTIDE RECORDS
# ============================================================

st.header("2️⃣ Nucleotide Records")

nuccore_ids = get_nucleotide_links(
    selected_gene_id
)

if not nuccore_ids:

    st.error(
        "No linked nucleotide records were found."
    )

    st.stop()


st.write(
    f"Found {len(nuccore_ids)} "
    f"linked nucleotide records."
)


nucleotide_summaries = (
    get_nucleotide_summaries(
        tuple(nuccore_ids)
    )
)


# ------------------------------------------------------------
# Find RefSeq mRNA records
# ------------------------------------------------------------

refseq_records = [
    record
    for record in nucleotide_summaries
    if record["accession"].startswith("NM_")
]


if not refseq_records:

    st.warning(
        "No RefSeq mRNA record was found "
        "among the linked nucleotide records."
    )

    st.stop()


if len(refseq_records) > 1:

    st.warning(
        "Multiple RefSeq mRNA records were found."
    )

    selected_accession = st.selectbox(
        "Select RefSeq transcript",
        [
            record["accession"]
            for record in refseq_records
        ]
    )

else:

    selected_accession = (
        refseq_records[0]["accession"]
    )


selected_record = next(
    record
    for record in refseq_records
    if record["accession"] == selected_accession
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

st.markdown(
    f"[🔗 Open RefSeq record](https://www.ncbi.nlm.nih.gov/nuccore/{selected_record['accession']})"
)


# ============================================================
# FETCH GENBANK
# ============================================================

genbank_record = fetch_one_genbank(
    selected_record["uid"]
)

if not genbank_record:

    st.error(
        "Could not retrieve the GenBank record."
    )

    st.stop()


st.success(
    f"Successfully retrieved "
    f"{selected_record['accession']}."
)


with st.expander("View GenBank record"):

    st.code(
        genbank_record,
        language="text"
    )


# ============================================================
# STEP 3 — CDS
# ============================================================

st.header("3️⃣ Coding Sequence")

(
    cds_sequence,
    cds_start,
    cds_end
) = extract_cds(
    genbank_record
)


if not cds_sequence:

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


# ------------------------------------------------------------
# Nucleotide FASTA
# ------------------------------------------------------------

nucleotide_fasta = (
    f">{locus_id}_CDS\n"
    f"{format_fasta(cds_sequence)}"
)

with st.expander(
    "🧬 View Nucleotide FASTA",
    expanded=True
):

    st.code(
        nucleotide_fasta,
        language="text"
    )


# ============================================================
# STEP 4 — PROTEIN
# ============================================================

st.header("4️⃣ Protein Sequence")

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


protein_fasta = (
    f">{locus_id}_protein\n"
    f"{format_fasta(protein_sequence)}"
)

with st.expander(
    "🧬 View Protein FASTA",
    expanded=True
):

    st.code(
        protein_fasta,
        language="text"
    )


# ============================================================
# STEP 5 — UNIPROT
# ============================================================

st.header("5️⃣ UniProt")


uniprot_data = search_uniprot(
    locus_tag,
    taxonomy_id
)


if not uniprot_data:

    st.warning(
        "UniProt search could not be completed."
    )

    uniprot_id = None
    uniprot_record = None
    uniprot_sequence = None
    sequences_match = False

else:

    uniprot_results = (
        uniprot_data.get(
            "results",
            []
        )
    )

    if not uniprot_results:

        st.warning(
            "No UniProt record was found."
        )

        uniprot_id = None
        uniprot_record = None
        uniprot_sequence = None
        sequences_match = False

    else:

        # ----------------------------------------------------
        # Multiple UniProt records
        # ----------------------------------------------------

        if len(uniprot_results) > 1:

            st.warning(
                "Multiple UniProt records were found."
            )

            uniprot_options = [
                result.get(
                    "primaryAccession"
                )
                for result in uniprot_results
            ]

            uniprot_id = st.selectbox(
                "Select UniProt record",
                uniprot_options
            )

        else:

            uniprot_id = (
                uniprot_results[0]
                .get("primaryAccession")
            )

        uniprot_record = get_uniprot_record(
            uniprot_id
        )

        if not uniprot_record:

            st.warning(
                "Could not retrieve the "
                "UniProt record."
            )

            uniprot_sequence = None
            sequences_match = False

        else:

            # ------------------------------------------------
            # UniProt information
            # ------------------------------------------------

            entry_name = uniprot_record.get(
                "uniProtkbId",
                "Not available"
            )

            sequence_info = (
                uniprot_record.get(
                    "sequence",
                    {}
                )
            )

            uniprot_sequence = (
                sequence_info.get(
                    "value",
                    ""
                )
            )

            uniprot_length = len(
                uniprot_sequence
            )

            sequences_match = (
                protein_sequence ==
                uniprot_sequence
            )


            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric(
                    "UniProt ID",
                    uniprot_id
                )

            with col2:
                st.metric(
                    "Entry Name",
                    entry_name
                )

            with col3:
                st.metric(
                    "Protein Length",
                    f"{uniprot_length} aa"
                )


            st.write(
                f"**Our translated protein:** "
                f"{len(protein_sequence)} aa"
            )

            st.write(
                f"**UniProt protein:** "
                f"{uniprot_length} aa"
            )


            if sequences_match:

                st.success(
                    "Protein sequences are identical."
                )

            else:

                st.warning(
                    "Protein sequences are not identical."
                )


            st.markdown(
                f"[🔗 Open UniProt record](https://www.uniprot.org/uniprotkb/{uniprot_id})"
            )


            uniprot_fasta = (
                f">{uniprot_id}\n"
                f"{format_fasta(uniprot_sequence)}"
            )

            with st.expander(
                "View UniProt protein sequence"
            ):

                st.code(
                    uniprot_fasta,
                    language="text"
                )


# ============================================================
# ANALYSIS SUMMARY
# ============================================================

st.divider()

st.header("📊 Analysis Summary")

summary_table = pd.DataFrame([
    {
        "Gene": gene_name,
        "Organism": organism_name,
        "NCBI Gene ID": selected_gene_id,
        "RefSeq": selected_record["accession"],
        "CDS Length": f"{len(cds_sequence)} nt",
        "Protein Length": f"{len(protein_sequence)} aa",
        "UniProt": (
            uniprot_id
            if uniprot_id
            else "Not found"
        )
    }
])


st.dataframe(
    summary_table,
    use_container_width=True,
    hide_index=True
)


# ============================================================
# STEP 6 — INTERPRO
# ============================================================

st.header("6️⃣ InterPro")


interpro_data = None
interpro_results = []


if uniprot_id:

    interpro_data = get_interpro_annotations(
        uniprot_id
    )

    if interpro_data:

        interpro_results = (
            interpro_data.get(
                "results",
                []
            )
        )


if interpro_results:

    st.success(
        f"Found {len(interpro_results)} "
        f"InterPro annotation(s)."
    )

    # --------------------------------------------------------
    # Display each InterPro result
    # --------------------------------------------------------

    for result in interpro_results:

        details = extract_interpro_details(
            result
        )

        st.subheader(
            details["name"]
        )

        st.write(
            f"**InterPro ID:** "
            f"{details['accession']}"
        )

        st.write(
            f"**Type:** "
            f"{details['type']}"
        )

        # ----------------------------------------------------
        # Member databases
        # ----------------------------------------------------

        if details["member_databases"]:

            st.write(
                "**Member databases:**"
            )

            for member in (
                details["member_databases"]
            ):

                st.write(
                    f"- {member}"
                )

        # ----------------------------------------------------
        # Matched regions
        # ----------------------------------------------------

        if details["matched_regions"]:

            st.write(
                "**Matched regions:**"
            )

            for region in (
                details["matched_regions"]
            ):

                st.write(
                    f"- {region}"
                )

        # ----------------------------------------------------
        # GO terms
        # ----------------------------------------------------

        if details["go_terms"]:

            st.write(
                "**GO terms:**"
            )

            for go_term in (
                details["go_terms"]
            ):

                st.write(
                    f"- {go_term}"
                )

        st.markdown(
            f"[🔗 Open InterPro entry](https://www.ebi.ac.uk/interpro/entry/InterPro/{details['accession']})"
        )

        st.divider()


else:

    st.info(
        "No InterPro annotations were found "
        "for this UniProt protein."
    )


# ============================================================
# STEP 7 — FINAL ANNOTATION TABLE
# ============================================================

st.header("7️⃣ Final Annotation Table")


annotation_rows = []


# ------------------------------------------------------------
# NCBI
# ------------------------------------------------------------

annotation_rows.append({
    "Source": "NCBI",
    "Record ID": selected_gene_id,
    "Annotation": description,
    "Organism": organism_name,
    "Details": "Gene"
})


# ------------------------------------------------------------
# RefSeq
# ------------------------------------------------------------

annotation_rows.append({
    "Source": "RefSeq",
    "Record ID": selected_record["accession"],
    "Annotation": selected_record["title"],
    "Organism": organism_name,
    "Details": (
        f"CDS: {cds_start}..{cds_end}; "
        f"Length: {len(cds_sequence)} nt"
    )
})


# ------------------------------------------------------------
# Protein / RefSeq
# ------------------------------------------------------------

annotation_rows.append({
    "Source": "RefSeq Protein",
    "Record ID": "NP_171609.1",
    "Annotation": "Translated protein",
    "Organism": organism_name,
    "Details": (
        f"Protein length: "
        f"{len(protein_sequence)} aa"
    )
})


# ------------------------------------------------------------
# UniProt
# ------------------------------------------------------------

if uniprot_id:

    annotation_rows.append({
        "Source": "UniProt",
        "Record ID": uniprot_id,
        "Annotation": (
            uniprot_record.get(
                "uniProtkbId",
                "Not available"
            )
            if uniprot_record
            else "Not available"
        ),
        "Organism": organism_name,
        "Details": (
            f"Protein length: "
            f"{len(uniprot_sequence)} aa; "
            f"Sequence match: "
            f"{'Yes' if sequences_match else 'No'}"
        )
    })


# ------------------------------------------------------------
# InterPro
# ------------------------------------------------------------

for result in interpro_results:

    details = extract_interpro_details(
        result
    )

    extra_details = details["type"]

    if details["member_databases"]:

        extra_details += (
            "; Members: "
            + ", ".join(
                details["member_databases"]
            )
        )

    annotation_rows.append({
        "Source": "InterPro",
        "Record ID": details["accession"],
        "Annotation": details["name"],
        "Organism": organism_name,
        "Details": extra_details
    })


annotation_table = pd.DataFrame(
    annotation_rows
)


st.dataframe(
    annotation_table,
    use_container_width=True,
    hide_index=True
)


# ============================================================
# STEP 8 — DOWNLOAD RESULTS
# ============================================================

st.header("8️⃣ Download Results")


# ------------------------------------------------------------
# FASTA downloads
# ------------------------------------------------------------

download_col1, download_col2 = st.columns(2)


with download_col1:

    st.download_button(
        label="⬇️ Download Nucleotide FASTA",
        data=nucleotide_fasta,
        file_name=(
            f"{locus_id}_CDS.fasta"
        ),
        mime="text/plain",
        use_container_width=True
    )


with download_col2:

    st.download_button(
        label="⬇️ Download Protein FASTA",
        data=protein_fasta,
        file_name=(
            f"{locus_id}_protein.fasta"
        ),
        mime="text/plain",
        use_container_width=True
    )


# ------------------------------------------------------------
# Annotation downloads
# ------------------------------------------------------------

annotation_csv = annotation_table.to_csv(
    index=False
)

annotation_json = annotation_table.to_json(
    orient="records",
    indent=2
)


download_col3, download_col4 = st.columns(2)


with download_col3:

    st.download_button(
        label="⬇️ Download Annotation CSV",
        data=annotation_csv,
        file_name="annotation_results.csv",
        mime="text/csv",
        use_container_width=True
    )


with download_col4:

    st.download_button(
        label="⬇️ Download Annotation JSON",
        data=annotation_json,
        file_name="annotation_results.json",
        mime="application/json",
        use_container_width=True
    )


# ============================================================
# FOOTER
# ============================================================

st.divider()

st.caption(
    "Bio API Tester | NCBI • UniProt • InterPro"
)

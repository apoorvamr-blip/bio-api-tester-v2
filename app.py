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
# Replace with your real email before deployment.
NCBI_EMAIL = "your_email@example.com"

NCBI_TOOL = "BioAPITester"

UNIPROT_SEARCH_URL = (
    "https://rest.uniprot.org/uniprotkb/search"
)

INTERPRO_BASE_URL = (
    "https://www.ebi.ac.uk/interpro/api"
)


# ============================================================
# GENERAL HELPERS
# ============================================================

def format_fasta(sequence, width=60):
    """Format sequence into readable FASTA lines."""

    return "\n".join(
        sequence[i:i + width]
        for i in range(0, len(sequence), width)
    )


def ncbi_get(endpoint, params, timeout=30):
    """Make a request to NCBI E-utilities."""

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
# NCBI GENE
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

    return data.get(
        "esearchresult",
        {}
    ).get(
        "idlist",
        []
    )


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

    return data.get(
        "result",
        {}
    ).get(
        str(gene_id),
        {}
    )


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

    linked_ids = []

    for linkset in data.get("linksets", []):

        for link in linkset.get(
            "linksetdbs",
            []
        ):

            linked_ids.extend(
                link.get(
                    "links",
                    []
                )
            )

    return list(
        dict.fromkeys(
            linked_ids
        )
    )


@st.cache_data(ttl=3600)
def get_nucleotide_summaries(
    nuccore_ids
):

    if not nuccore_ids:
        return []

    params = {
        "db": "nuccore",
        "id": ",".join(
            nuccore_ids
        ),
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
            item.get(
                "accessionversion"
            )
            or item.get(
                "caption"
            )
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
def fetch_one_genbank(
    nuccore_id
):

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
# GENBANK PARSING
# ============================================================

def extract_cds(
    genbank_record
):
    """
    Extract CDS sequence, coordinates,
    protein ID and product from GenBank.
    """

    lines = genbank_record.splitlines()

    cds_start = None
    cds_end = None

    protein_id = None
    product = None

    # --------------------------------------------------------
    # Find CDS location
    # --------------------------------------------------------

    for line in lines:

        if line.startswith("     CDS"):

            location = line[21:].strip()

            numbers = re.findall(
                r"\d+",
                location
            )

            if len(numbers) >= 2:

                cds_start = int(
                    numbers[0]
                )

                cds_end = int(
                    numbers[-1]
                )

            break

    if cds_start is None:
        return None, None, None, None, None

    # --------------------------------------------------------
    # Find CDS qualifiers
    # --------------------------------------------------------

    in_cds = False

    for line in lines:

        if line.startswith("     CDS"):

            in_cds = True
            continue

        if in_cds:

            # End of CDS feature
            if (
                line.startswith("     ")
                and not line.startswith(
                    "                     "
                )
            ):
                break

            stripped = line.strip()

            if stripped.startswith(
                "/protein_id="
            ):

                protein_id = (
                    stripped
                    .split(
                        "=",
                        1
                    )[1]
                    .strip()
                    .strip('"')
                )

            elif stripped.startswith(
                "/product="
            ):

                product = (
                    stripped
                    .split(
                        "=",
                        1
                    )[1]
                    .strip()
                    .strip('"')
                )

    # --------------------------------------------------------
    # Find ORIGIN
    # --------------------------------------------------------

    origin_index = None

    for i, line in enumerate(lines):

        if line.startswith("ORIGIN"):

            origin_index = i
            break

    if origin_index is None:

        return (
            None,
            None,
            None,
            protein_id,
            product
        )

    # --------------------------------------------------------
    # Extract nucleotide sequence
    # --------------------------------------------------------

    sequence_parts = []

    for line in lines[
        origin_index + 1:
    ]:

        if line.startswith("//"):
            break

        parts = line.split()

        if parts:

            sequence_parts.extend(
                parts[1:]
            )

    nucleotide_sequence = (
        "".join(
            sequence_parts
        ).upper()
    )

    cds_sequence = (
        nucleotide_sequence[
            cds_start - 1:cds_end
        ]
    )

    return (
        cds_sequence,
        cds_start,
        cds_end,
        protein_id,
        product
    )


def extract_locus_tag(
    genbank_record
):

    for line in genbank_record.splitlines():

        stripped = line.strip()

        if stripped.startswith(
            "/locus_tag="
        ):

            return (
                stripped
                .split(
                    "=",
                    1
                )[1]
                .strip()
                .strip('"')
            )

    return None


def translate_cds(
    cds_sequence
):

    protein_sequence = str(
        Seq(
            cds_sequence
        ).translate()
    )

    if protein_sequence.endswith("*"):

        protein_sequence = (
            protein_sequence[:-1]
        )

    return protein_sequence


# ============================================================
# UNIPROT
# ============================================================

@st.cache_data(ttl=3600)
def search_uniprot(
    locus_tag,
    taxonomy_id
):

    queries = [
        (
            f'gene:{locus_tag} '
            f'AND organism_id:{taxonomy_id}'
        ),
        (
            f'"{locus_tag}" '
            f'AND organism_id:{taxonomy_id}'
        )
    ]

    for query in queries:

        params = {
            "query": query,
            "format": "json",
            "size": 5
        }

        response = requests.get(
            UNIPROT_SEARCH_URL,
            params=params,
            timeout=30
        )

        if response.status_code != 200:
            continue

        data = response.json()

        if data.get(
            "results"
        ):

            return data

    return None


@st.cache_data(ttl=3600)
def search_uniprot_by_protein(
    protein_id
):

    if not protein_id:
        return None

    params = {
        "query": (
            f"xref:RefSeq_{protein_id}"
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
def get_uniprot_record(
    uniprot_id
):

    url = (
        "https://rest.uniprot.org/"
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
# INTERPRO
# ============================================================

@st.cache_data(ttl=3600)
def get_interpro_annotations(
    uniprot_id
):

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


# ============================================================
# UNIPROT / INTERPRO ANNOTATION PARSERS
# ============================================================

def get_uniprot_annotation_details(record):
    """Extract UniProt's own annotation/evidence fields."""

    if not record:
        return {
            "annotation_score": "Not available",
            "protein_name": "Not available",
            "function": "Not available",
            "go_terms": [],
            "evidence": []
        }

    annotation_score = record.get(
        "annotationScore",
        "Not available"
    )

    protein_name = "Not available"
    description = record.get(
        "proteinDescription",
        {}
    )

    recommended = description.get(
        "recommendedName",
        {}
    )

    full_name = recommended.get(
        "fullName",
        {}
    )

    if isinstance(full_name, dict):
        protein_name = full_name.get(
            "value",
            "Not available"
        )
    elif full_name:
        protein_name = str(full_name)

    function_texts = []
    evidence_texts = []

    for comment in record.get("comments", []):
        if not isinstance(comment, dict):
            continue

        comment_type = comment.get("commentType", "")

        if comment_type == "FUNCTION":
            texts = comment.get("texts", [])
            for item in texts:
                if isinstance(item, dict):
                    value = item.get("value")
                    if value:
                        function_texts.append(value)

        elif comment_type in {
            "CATALYTIC ACTIVITY",
            "PATHWAY",
            "SUBUNIT",
            "SUBCELLULAR LOCATION"
        }:
            texts = comment.get("texts", [])
            for item in texts:
                if isinstance(item, dict):
                    value = item.get("value")
                    if value:
                        evidence_texts.append(
                            f"{comment_type}: {value}"
                        )

    go_terms = []

    for xref in record.get(
        "uniProtKBCrossReferences",
        []
    ):
        if not isinstance(xref, dict):
            continue
        if xref.get("database") != "GO":
            continue

        properties = xref.get("properties", [])
        term = xref.get("id", "")
        for prop in properties:
            if isinstance(prop, dict):
                key = prop.get("key", "")
                value = prop.get("value", "")
                if key and value:
                    term = f"{term} ({value})"
                    break
        if term:
            go_terms.append(term)

    return {
        "annotation_score": annotation_score,
        "protein_name": protein_name,
        "function": " ".join(function_texts) if function_texts else "Not available",
        "go_terms": go_terms,
        "evidence": evidence_texts
    }


def _recursive_find_values(obj, target_keys):
    """Find values for selected keys anywhere in a nested API response."""

    found = []

    if isinstance(obj, dict):
        for key, value in obj.items():
            if key.lower() in target_keys:
                found.append(value)
            found.extend(
                _recursive_find_values(
                    value,
                    target_keys
                )
            )

    elif isinstance(obj, list):
        for item in obj:
            found.extend(
                _recursive_find_values(
                    item,
                    target_keys
                )
            )

    return found


def _format_interpro_locations(result):
    """Extract InterPro matched amino-acid regions when present."""

    locations = _recursive_find_values(
        result,
        {"location", "locations", "fragments"}
    )

    formatted = []

    def walk(value):
        if isinstance(value, dict):
            start = value.get("start")
            end = value.get("end")
            if start is not None and end is not None:
                formatted.append(f"{start}–{end} aa")
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)

    for location in locations:
        walk(location)

    return list(dict.fromkeys(formatted))


def _format_interpro_scores(result):
    """Extract score/E-value fields only when InterPro actually supplies them."""

    score_keys = {
        "score",
        "evalue",
        "e_value",
        "e-value",
        "bit_score",
        "bitscore"
    }

    values = []

    def walk(obj):
        if isinstance(obj, dict):
            for key, value in obj.items():
                if key.lower() in score_keys and value not in (None, ""):
                    values.append(
                        f"{key}: {value}"
                    )
                else:
                    walk(value)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(result)
    return list(dict.fromkeys(values))


def extract_interpro_details(result):
    """Extract InterPro entry, member signatures, regions and native evidence."""

    metadata = result.get("metadata", {})

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

    member_databases = []

    entries = result.get("entries", [])

    if isinstance(entries, list):
        for entry in entries:
            if not isinstance(entry, dict):
                continue

            accession_value = (
                entry.get("accession")
                or entry.get("member_accession")
            )

            name_value = (
                entry.get("name")
                or entry.get("member_name")
            )

            if accession_value:
                if name_value:
                    member_databases.append(
                        f"{accession_value} ({name_value})"
                    )
                else:
                    member_databases.append(
                        str(accession_value)
                    )

    scores = _format_interpro_scores(result)
    locations = _format_interpro_locations(result)

    return {
        "accession": accession,
        "name": name,
        "type": entry_type,
        "member_databases": list(dict.fromkeys(member_databases)),
        "scores": scores,
        "locations": locations
    }


def build_final_annotation(
    description,
    uniprot_details,
    interpro_details_list
):
    """Create a conservative putative name/function from database evidence."""

    original_text = (description or "").lower()
    uniprot_name = uniprot_details.get(
        "protein_name",
        "Not available"
    )

    hypothetical_terms = {
        "hypothetical",
        "uncharacterized",
        "unknown function",
        "putative protein"
    }

    has_specific_uniprot_name = (
        uniprot_name
        and uniprot_name != "Not available"
        and not any(
            term in uniprot_name.lower()
            for term in hypothetical_terms
        )
    )

    if has_specific_uniprot_name:
        final_name = uniprot_name
        name_basis = "UniProt protein annotation"
    else:
        domain_name = None
        for item in interpro_details_list:
            candidate = item.get("name", "")
            if candidate and candidate != "Not available":
                domain_name = candidate
                break

        if domain_name:
            final_name = f"Putative {domain_name}"
            name_basis = "InterPro domain/family evidence"
        elif "hypothetical" in original_text or "uncharacterized" in original_text:
            final_name = "Putative uncharacterized protein"
            name_basis = "No specific domain name available"
        else:
            final_name = description or "Protein annotation unavailable"
            name_basis = "NCBI annotation"

    function = uniprot_details.get(
        "function",
        "Not available"
    )

    if function == "Not available":
        go_terms = uniprot_details.get("go_terms", [])
        if go_terms:
            function = (
                "Putative function supported by UniProt Gene Ontology "
                "annotations: " + "; ".join(go_terms[:5])
            )
        elif interpro_details_list:
            domain_names = [
                item.get("name")
                for item in interpro_details_list
                if item.get("name")
                and item.get("name") != "Not available"
            ]
            if domain_names:
                function = (
                    "Putative function associated with the detected "
                    "InterPro domain/family: "
                    + "; ".join(dict.fromkeys(domain_names[:3]))
                )
            else:
                function = "Function could not be assigned from available database evidence."
        else:
            function = "Function could not be assigned from available database evidence."

    return {
        "final_name": final_name,
        "function": function,
        "basis": name_basis
    }


# ============================================================
# HEADER
# ============================================================

st.title(
    "🧬 Bio API Tester"
)

st.write(
    "Retrieve, analyze, and annotate biological "
    "sequences using NCBI, UniProt, and InterPro."
)

st.markdown(
    "### Workflow"
)

st.markdown(
    """
`Locus ID` → `NCBI Gene` → `RefSeq` → `CDS`
→ `Protein` → `UniProt` → `InterPro` → `Results`
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


# ============================================================
# STEP 1 — NCBI GENE
# ============================================================

st.header(
    "1️⃣ Gene Search"
)

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
# Multiple NCBI records
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

st.markdown(
    f"[🔗 Open NCBI Gene record]"
    f"(https://www.ncbi.nlm.nih.gov/gene/"
    f"{selected_gene_id})"
)


# ============================================================
# STEP 2 — LINKED NUCLEOTIDE RECORDS
# ============================================================

st.header(
    "2️⃣ Nucleotide Records"
)

nuccore_ids = get_nucleotide_links(
    selected_gene_id
)

nucleotide_summaries = (
    get_nucleotide_summaries(
        tuple(nuccore_ids)
    )
)


st.write(
    f"Found {len(nucleotide_summaries)} "
    f"linked nucleotide records."
)


# ------------------------------------------------------------
# Accept BOTH NM_ and XM_
# ------------------------------------------------------------

refseq_records = [
    record
    for record in nucleotide_summaries
    if (
        record["accession"].startswith(
            "NM_"
        )
        or
        record["accession"].startswith(
            "XM_"
        )
    )
]


selected_record = None


if refseq_records:

    if len(refseq_records) > 1:

        st.warning(
            "Multiple RefSeq transcript records "
            "were found. Please select one."
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
        if record["accession"]
        == selected_accession
    )

    st.success(
        f"Selected transcript: "
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

else:

    st.warning(
        "No NM_ or XM_ RefSeq transcript "
        "was found among the linked records."
    )

    st.info(
        "The application will try to continue "
        "using other available records."
    )


# ============================================================
# REFSEQ PROTEIN DISCOVERY
# ============================================================

protein_records = [
    record
    for record in nucleotide_summaries
    if (
        record["accession"].startswith(
            "NP_"
        )
        or
        record["accession"].startswith(
            "XP_"
        )
    )
]


selected_protein_record = None


if protein_records:

    if len(protein_records) > 1:

        protein_accession = st.selectbox(
            "Select RefSeq protein",
            [
                record["accession"]
                for record in protein_records
            ]
        )

    else:

        protein_accession = (
            protein_records[0]["accession"]
        )

    selected_protein_record = next(
        record
        for record in protein_records
        if record["accession"]
        == protein_accession
    )

    st.write(
        f"**RefSeq protein:** "
        f"{selected_protein_record['accession']}"
    )


# ============================================================
# STEP 3 — GENBANK / CDS
# ============================================================

st.header(
    "3️⃣ Coding Sequence"
)

cds_sequence = None
cds_start = None
cds_end = None
protein_id = None
product = None
genbank_record = None


if selected_record:

    genbank_record = fetch_one_genbank(
        selected_record["uid"]
    )

    if genbank_record:

        (
            cds_sequence,
            cds_start,
            cds_end,
            protein_id,
            product
        ) = extract_cds(
            genbank_record
        )

        if cds_sequence:

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

            if protein_id:

                st.write(
                    f"**Protein ID:** "
                    f"{protein_id}"
                )

        else:

            st.warning(
                "The transcript was found, "
                "but a CDS could not be extracted."
            )

        with st.expander(
            "View GenBank record"
        ):

            st.code(
                genbank_record,
                language="text"
            )

    else:

        st.warning(
            "Could not retrieve the GenBank "
            "transcript record."
        )


# ============================================================
# NUCLEOTIDE FASTA
# ============================================================

nucleotide_fasta = None


if cds_sequence:

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

st.header(
    "4️⃣ Protein Sequence"
)

protein_sequence = None


# ------------------------------------------------------------
# Preferred route: translate CDS
# ------------------------------------------------------------

if cds_sequence:

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


# ------------------------------------------------------------
# Protein FASTA
# ------------------------------------------------------------

protein_fasta = None


if protein_sequence:

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

elif selected_protein_record:

    st.info(
        "A RefSeq protein record was found, "
        "but a transcript CDS was not available "
        "for direct translation."
    )

    st.markdown(
        f"[🔗 Open RefSeq protein record]"
        f"(https://www.ncbi.nlm.nih.gov/protein/"
        f"{selected_protein_record['accession']})"
    )

else:

    st.warning(
        "No protein sequence is currently available "
        "from the selected NCBI records."
    )


# ============================================================
# STEP 5 — UNIPROT
# ============================================================

st.header(
    "5️⃣ UniProt"
)

uniprot_data = None
uniprot_id = None
uniprot_record = None
uniprot_sequence = None
sequences_match = False
uniprot_details = {
    "annotation_score": "Not available",
    "protein_name": "Not available",
    "function": "Not available",
    "go_terms": [],
    "evidence": []
}


# ------------------------------------------------------------
# Search by locus tag first
# ------------------------------------------------------------

uniprot_data = search_uniprot(
    locus_id,
    taxonomy_id
)


# ------------------------------------------------------------
# If no result, try RefSeq protein accession
# ------------------------------------------------------------

if (
    not uniprot_data
    and protein_id
):

    uniprot_data = (
        search_uniprot_by_protein(
            protein_id
        )
    )


if uniprot_data:

    uniprot_results = (
        uniprot_data.get(
            "results",
            []
        )
    )

    if uniprot_results:

        if len(
            uniprot_results
        ) > 1:

            st.warning(
                "Multiple UniProt records "
                "were found."
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
                .get(
                    "primaryAccession"
                )
            )

        uniprot_record = (
            get_uniprot_record(
                uniprot_id
            )
        )

        if uniprot_record:

            uniprot_details = get_uniprot_annotation_details(
                uniprot_record
            )

            entry_name = (
                uniprot_record.get(
                    "uniProtkbId",
                    "Not available"
                )
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

            if protein_sequence:

                sequences_match = (
                    protein_sequence
                    == uniprot_sequence
                )

            col1, col2, col3 = (
                st.columns(3)
            )

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
                    f"{len(uniprot_sequence)} aa"
                )

            if protein_sequence:

                if sequences_match:

                    st.success(
                        "Protein sequences are identical."
                    )

                else:

                    st.warning(
                        "Protein sequences are not identical."
                    )

            st.write(
                f"**UniProt Annotation Score:** "
                f"{uniprot_details['annotation_score']} / 5"
                if uniprot_details["annotation_score"] != "Not available"
                else "**UniProt Annotation Score:** Not available"
            )

            st.write(
                f"**UniProt Protein Name:** "
                f"{uniprot_details['protein_name']}"
            )

            if uniprot_details["function"] != "Not available":
                st.write(
                    f"**UniProt Function:** "
                    f"{uniprot_details['function']}"
                )

            if uniprot_details["evidence"]:
                with st.expander("View UniProt functional evidence"):
                    for item in uniprot_details["evidence"]:
                        st.write(f"- {item}")

            if uniprot_details["go_terms"]:
                with st.expander("View UniProt GO annotations"):
                    for term in uniprot_details["go_terms"]:
                        st.write(f"- {term}")

            st.markdown(
                f"[🔗 Open UniProt record]"
                f"(https://www.uniprot.org/uniprotkb/"
                f"{uniprot_id})"
            )

        else:

            st.warning(
                "UniProt record could not be retrieved."
            )

    else:

        st.info(
            "No UniProt record was found."
        )

else:

    st.info(
        "No UniProt record was found for this gene."
    )


# ============================================================
# ANALYSIS SUMMARY
# ============================================================

st.divider()

st.header(
    "📊 Analysis Summary"
)


summary_table = pd.DataFrame([
    {
        "Gene": gene_name,
        "Organism": organism_name,
        "NCBI Gene ID": selected_gene_id,
        "RefSeq": (
            selected_record["accession"]
            if selected_record
            else "Not available"
        ),
        "Protein ID": (
            protein_id
            if protein_id
            else (
                selected_protein_record[
                    "accession"
                ]
                if selected_protein_record
                else "Not available"
            )
        ),
        "CDS Length": (
            f"{len(cds_sequence)} nt"
            if cds_sequence
            else "Not available"
        ),
        "Protein Length": (
            f"{len(protein_sequence)} aa"
            if protein_sequence
            else (
                "Not available"
            )
        ),
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

st.header(
    "6️⃣ InterPro"
)

interpro_data = None
interpro_results = []
parsed_interpro_details = []


if uniprot_id:

    interpro_data = (
        get_interpro_annotations(
            uniprot_id
        )
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
        f"Found {len(interpro_results)} InterPro annotation(s)."
    )

    parsed_interpro_details = []

    for result in interpro_results:

        details = extract_interpro_details(result)
        parsed_interpro_details.append(details)

    interpro_display_rows = []

    for details in parsed_interpro_details:
        interpro_display_rows.append({
            "InterPro ID": details["accession"],
            "Domain / Family": details["name"],
            "Type": details["type"],
            "InterPro Match Evidence": (
                "; ".join(details["scores"])
                if details["scores"]
                else "Not provided by InterPro API"
            ),
            "Matched Region": (
                "; ".join(details["locations"])
                if details["locations"]
                else "Not provided"
            ),
            "Member Database": (
                "; ".join(details["member_databases"])
                if details["member_databases"]
                else "Not available"
            )
        })

    interpro_display_table = pd.DataFrame(
        interpro_display_rows
    )

    st.dataframe(
        interpro_display_table,
        use_container_width=True,
        hide_index=True
    )

    for details in parsed_interpro_details:
        with st.expander(
            f"{details['accession']} — {details['name']}"
        ):
            st.write(
                f"**Type:** {details['type']}"
            )

            st.write(
                "**InterPro match evidence:** "
                + (
                    "; ".join(details["scores"])
                    if details["scores"]
                    else "Not provided by InterPro API"
                )
            )

            st.write(
                "**Matched region:** "
                + (
                    "; ".join(details["locations"])
                    if details["locations"]
                    else "Not provided"
                )
            )

            if details["member_databases"]:
                st.write("**Member databases:**")
                for member in details["member_databases"]:
                    st.write(f"- {member}")

            st.markdown(
                f"[🔗 Open InterPro entry]"
                f"(https://www.ebi.ac.uk/interpro/"
                f"entry/InterPro/"
                f"{details['accession']})"
            )

else:

    st.info(
        "No InterPro annotations were found "
        "for the available UniProt protein."
    )


# ============================================================
# STEP 7 — FINAL ANNOTATION TABLE
# ============================================================

st.header(
    "7️⃣ Final Annotation Table"
)

# ------------------------------------------------------------
# Final functional interpretation
# ------------------------------------------------------------

final_annotation = build_final_annotation(
    description,
    uniprot_details,
    parsed_interpro_details
)

st.subheader("🧬 Final Functional Annotation")

final_annotation_table = pd.DataFrame([
    {
        "Field": "Final Protein Name",
        "Result": final_annotation["final_name"]
    },
    {
        "Field": "Predicted / Reported Function",
        "Result": final_annotation["function"]
    },
    {
        "Field": "Annotation Basis",
        "Result": final_annotation["basis"]
    },
    {
        "Field": "UniProt Annotation Score",
        "Result": (
            f"{uniprot_details['annotation_score']} / 5"
            if uniprot_details["annotation_score"] != "Not available"
            else "Not available"
        )
    },
    {
        "Field": "InterPro Domains / Families",
        "Result": (
            "; ".join(
                item["name"]
                for item in parsed_interpro_details
                if item.get("name")
                and item["name"] != "Not available"
            )
            if parsed_interpro_details
            else "None found"
        )
    }
])

st.dataframe(
    final_annotation_table,
    use_container_width=True,
    hide_index=True
)

# ------------------------------------------------------------
# Evidence table combining NCBI, UniProt and InterPro
# ------------------------------------------------------------

annotation_rows = []

annotation_rows.append({
    "Source": "NCBI",
    "Record ID": selected_gene_id,
    "Annotation": description,
    "Organism": organism_name,
    "Native Match / Score": "Not applicable",
    "Matched Region": "—",
    "Function / Evidence": "Gene-level annotation"
})

if selected_record:
    annotation_rows.append({
        "Source": "RefSeq Transcript",
        "Record ID": selected_record["accession"],
        "Annotation": selected_record["title"],
        "Organism": organism_name,
        "Native Match / Score": "Not applicable",
        "Matched Region": (
            f"{cds_start}..{cds_end}"
            if cds_sequence
            else "Not available"
        ),
        "Function / Evidence": "RefSeq transcript / CDS"
    })

if final_protein_id:
    annotation_rows.append({
        "Source": "RefSeq Protein",
        "Record ID": final_protein_id,
        "Annotation": product or "RefSeq protein",
        "Organism": organism_name,
        "Native Match / Score": "Not applicable",
        "Matched Region": "Full protein",
        "Function / Evidence": (
            f"Protein length: {len(protein_sequence)} aa"
            if protein_sequence
            else "Protein sequence available"
        )
    })

if uniprot_id:
    annotation_rows.append({
        "Source": "UniProt",
        "Record ID": uniprot_id,
        "Annotation": uniprot_details["protein_name"],
        "Organism": organism_name,
        "Native Match / Score": (
            f"Annotation Score: {uniprot_details['annotation_score']}/5"
            if uniprot_details["annotation_score"] != "Not available"
            else "Annotation Score: Not available"
        ),
        "Matched Region": "Full protein",
        "Function / Evidence": uniprot_details["function"]
    })
else:
    annotation_rows.append({
        "Source": "UniProt",
        "Record ID": "Not found",
        "Annotation": "Not available",
        "Organism": organism_name,
        "Native Match / Score": "Not available",
        "Matched Region": "—",
        "Function / Evidence": "No UniProt record found"
    })

for details in parsed_interpro_details:
    annotation_rows.append({
        "Source": "InterPro",
        "Record ID": details["accession"],
        "Annotation": details["name"],
        "Organism": organism_name,
        "Native Match / Score": (
            "; ".join(details["scores"])
            if details["scores"]
            else "Not provided by InterPro API"
        ),
        "Matched Region": (
            "; ".join(details["locations"])
            if details["locations"]
            else "Not provided"
        ),
        "Function / Evidence": (
            "; ".join(details["member_databases"])
            if details["member_databases"]
            else details["type"]
        )
    })

annotation_table = pd.DataFrame(
    annotation_rows
)

st.subheader("📋 Database Evidence Summary")
st.dataframe(
    annotation_table,
    use_container_width=True,
    hide_index=True
)

# STEP 8 — DOWNLOADS
# ============================================================

st.header(
    "8️⃣ Download Results"
)


download_col1, download_col2 = (
    st.columns(2)
)


with download_col1:

    if nucleotide_fasta:

        st.download_button(
            label="⬇️ Download Nucleotide FASTA",
            data=nucleotide_fasta,
            file_name=(
                f"{locus_id}_CDS.fasta"
            ),
            mime="text/plain",
            use_container_width=True
        )

    else:

        st.button(
            "⬇️ Nucleotide FASTA unavailable",
            disabled=True,
            use_container_width=True
        )


with download_col2:

    if protein_fasta:

        st.download_button(
            label="⬇️ Download Protein FASTA",
            data=protein_fasta,
            file_name=(
                f"{locus_id}_protein.fasta"
            ),
            mime="text/plain",
            use_container_width=True
        )

    else:

        st.button(
            "⬇️ Protein FASTA unavailable",
            disabled=True,
            use_container_width=True
        )


annotation_csv = (
    annotation_table.to_csv(
        index=False
    )
)

annotation_json = (
    annotation_table.to_json(
        orient="records",
        indent=2
    )
)


download_col3, download_col4 = (
    st.columns(2)
)


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

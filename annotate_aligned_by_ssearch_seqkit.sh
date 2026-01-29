#!/usr/bin/env bash
# annotate_aligned_by_ssearch_seqkit.sh (MULTI-DATABASE VERSION)
#
# Annotate aligned sequences by searching against reference database(s) using ssearch36.
# Produces annotated alignments with top hit names, strand, and bit scores.
#
# USAGE:
#   Single database:
#     ./annotate_aligned_by_ssearch_seqkit.sh <aligned.fa> <db.fasta> [out_prefix] [threads]
#
#   Multiple databases:
#     ./annotate_aligned_by_ssearch_seqkit.sh <aligned.fa> <db1.fa,db2.fa,db3.fa> [out_prefix] [threads]
#     Or with custom names:
#     ./annotate_aligned_by_ssearch_seqkit.sh <aligned.fa> <db1.fa:Name1,db2.fa:Name2> [out_prefix] [threads]
#
# EXAMPLES:
#   # Single database
#   ./annotate_aligned_by_ssearch_seqkit.sh sample.clw consensus_db.fasta sample.annot 8
#
#   # Multiple databases (auto-naming from filenames)
#   ./annotate_aligned_by_ssearch_seqkit.sh sample.clw repbase.fa,dfam.fa sample.annot 8
#
#   # Multiple databases (custom names)
#   ./annotate_aligned_by_ssearch_seqkit.sh sample.clw repbase.fa:RepBase,dfam.fa:Dfam sample.annot 8
#
# OUTPUTS:
#   <prefix>.aligned.annot.fasta            - All sequences (original order)
#   <prefix>.aligned.annot.sorted.fasta     - All sequences (sorted by hit name)
#   <prefix>.best_hits.tsv                  - TSV: qid, subject, bitscore, strand
#   <prefix>.<DB>.aligned.fasta             - Per-database alignments (multi-db mode)
#   <prefix>.unmatched.aligned.fasta        - Sequences with no match (multi-db mode)
#
# HEADER FORMAT:
#   Matched: >{DB1}hit_name|+|bs=450|original_id
#   Unmatched: >original_id (unchanged)
#
# Requirements: seqkit, ssearch36 (FASTA), awk, xargs, python3
set -euo pipefail

# Parse arguments
ALIGN=""
DB=""
OUT_PREFIX=""
THREADS=""

while [[ $# -gt 0 ]]; do
    if [ -z "$ALIGN" ]; then
        ALIGN="$1"
    elif [ -z "$DB" ]; then
        DB="$1"
    elif [ -z "$OUT_PREFIX" ]; then
        OUT_PREFIX="$1"
    elif [ -z "$THREADS" ]; then
        THREADS="$1"
    else
        echo "Unknown argument: $1" >&2
        exit 1
    fi
    shift
done

# Set defaults
OUT_PREFIX="${OUT_PREFIX:-${ALIGN%.*}.annot}"
THREADS="${THREADS:-$(nproc)}"

if [ -z "$ALIGN" ] || [ -z "$DB" ]; then
    echo "Usage: $0 <aligned.fa> <db.fasta|db1.fa,db2.fa> [out_prefix] [threads]" >&2
    exit 1
fi

command -v seqkit >/dev/null 2>&1 || { echo "seqkit required but not found" >&2; exit 1; }
command -v ssearch36 >/dev/null 2>&1 || { echo "ssearch36 required but not found" >&2; exit 1; }
command -v xargs >/dev/null 2>&1 || { echo "xargs required but not found" >&2; exit 1; }

WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/annotate_split_XXXX")
trap 'rm -rf "$WORKDIR"' EXIT

# Detect if multiple databases specified (comma-separated)
MULTI_DB=false
SPLIT_DBS=""
MERGED_DB=""

if [[ "$DB" == *","* ]]; then
    MULTI_DB=true
    echo "Multiple databases detected, preparing merged database..."
    MERGED_DB="${WORKDIR}/merged_db.fasta"
    
    # Parse database specifications: file.fa or file.fa:Name
    IFS=',' read -ra DB_SPECS <<< "$DB"
    DB_NAMES=()
    
    for spec in "${DB_SPECS[@]}"; do
        if [[ "$spec" == *":"* ]]; then
            # Custom name provided: file.fa:Name
            db_file="${spec%%:*}"
            db_name="${spec##*:}"
        else
            # Auto-generate name from filename
            db_file="$spec"
            db_name=$(basename "$db_file" | sed 's/\.[^.]*$//' | tr '[:lower:]' '[:upper:]')
        fi
        
        if [ ! -f "$db_file" ]; then
            echo "ERROR: Database file not found: $db_file" >&2
            exit 1
        fi
        
        echo "  Adding $db_file as {$db_name}"
        DB_NAMES+=("$db_name")
        
        # Add prefix and append to merged database
        sed "s/^>/>\\{${db_name}\\}/" "$db_file" >> "$MERGED_DB"
    done
    
    # Set SPLIT_DBS for later use
    SPLIT_DBS=$(IFS=','; echo "${DB_NAMES[*]}")
    DB="$MERGED_DB"
    
    echo "  Created merged database: $MERGED_DB"
    echo "  Database names: $SPLIT_DBS"
fi

# 1) produce aligned FASTA (preserve gaps) for final renaming step
ALIGNED_FA="${WORKDIR}/aligned_with_gaps.fasta"
seqkit seq "$ALIGN" > "$ALIGNED_FA"

# 2) produce ungapped FASTA (seqkit seq -g) and split into single-sequence files
UNALIGNED_FA="${WORKDIR}/unaligned_ungapped.fasta"
SPLITDIR="${WORKDIR}/singles"
mkdir -p "$SPLITDIR"
seqkit seq -g "$ALIGN" > "$UNALIGNED_FA"
seqkit split -s 1 -O "$SPLITDIR" "$UNALIGNED_FA" 2>/dev/null

# 3) per-file worker: run ssearch36 and write a .hit file with qid<TAB>subject<TAB>bitscore<TAB>strand
PROCESS_SCRIPT="${WORKDIR}/process_query.sh"
cat > "$PROCESS_SCRIPT" <<'SH'
#!/usr/bin/env bash
f="$1"
db="$2"
qid=$(sed -n '1s/^>//p' "$f" | awk '{print $1}')

# run ssearch36 (tabular -m 8); pick best hit by highest numeric value in last column (NF)
if ! ssearch_output=$(ssearch36 -m 8 "$f" "$db" 2>&1); then
    echo "WARNING: ssearch36 failed for $qid" >&2
    echo -e "${qid}\tNone\t0\t+" > "${f}.hit"
    exit 0
fi

# Parse output and get best hit with strand information
# In ssearch36 -m 8 format: qid, sid, %id, alnlen, mismatch, gaps, qstart, qend, sstart, send, evalue, bitscore
echo "$ssearch_output" | awk -F"\t" -v q="$qid" '
BEGIN{best_s="None";best=0;found=0;strand="+"} 
!/^#/ && NF>=12 {
    score=$(NF)+0
    # Determine strand: if sstart > send, it is reverse complement
    if(NF>=10){
        sstart=$(NF-3)
        send=$(NF-2)
        curr_strand = (sstart > send) ? "-" : "+"
    }
    if(found==0 || score>best){
        best=score
        best_s=$2
        strand=curr_strand
        found=1
    }
} 
END{ 
    if(found) 
        print q "\t" best_s "\t" int(best+0.5) "\t" strand
    else 
        print q "\tNone\t0\t+"
}' > "${f}.hit"
SH
chmod +x "$PROCESS_SCRIPT"

# 4) run all workers in parallel
NUM_SEQS=$(find "$SPLITDIR" -type f -name '*.fasta' | wc -l)
echo "Running ssearch36 in parallel on $NUM_SEQS single-sequence files ($THREADS threads)"

if [ "$NUM_SEQS" -eq 0 ]; then
    echo "ERROR: No sequences found to process" >&2
    exit 1
fi

find "$SPLITDIR" -type f -name "*.fasta" -print0 \
  | xargs -0 -n1 -P "$THREADS" -I {} bash -c '"$0" "{}" "'"$DB"'"' "$PROCESS_SCRIPT"

# 5) collect hits into single TSV (qid \t subject \t bitscore \t strand)
BEST_TSV="${OUT_PREFIX}.best_hits.tsv"

HIT_FILES=$(find "$SPLITDIR" -type f -name "*.hit" | wc -l)
if [ "$HIT_FILES" -eq 0 ]; then
    echo "ERROR: No hit files were generated" >&2
    exit 1
fi

find "$SPLITDIR" -type f -name "*.hit" -print0 | xargs -0 cat > "$BEST_TSV" || true

# Dedupe by keeping HIGHEST bitscore for each qid
awk -F"\t" '
{
    score=$3+0
    if(!($1 in best_score) || score>best_score[$1]){
        best[$1]=$0
        best_score[$1]=score
    }
} 
END{
    for(k in best) print best[k]
}' "$BEST_TSV" > "${BEST_TSV}.uniq" && mv "${BEST_TSV}.uniq" "$BEST_TSV"

# Validate that we got some hits
TOTAL_HITS=$(wc -l < "$BEST_TSV")
REAL_HITS=$(grep -v $'\tNone\t' "$BEST_TSV" | wc -l)

if [ "$REAL_HITS" -eq 0 ]; then
    echo "WARNING: No sequences matched the database! Check that database is appropriate." >&2
fi

# 6) annotate aligned FASTA headers
ANNOT_FA="${OUT_PREFIX}.aligned.annot.fasta"
ANNOT_FA_SORTED="${OUT_PREFIX}.aligned.annot.sorted.fasta"

# Pass SPLIT_DBS to Python script
python3 - "$ALIGNED_FA" "$BEST_TSV" "$ANNOT_FA" "$ANNOT_FA_SORTED" "$OUT_PREFIX" "$SPLIT_DBS" <<'PY'
import sys
import re
if len(sys.argv)!=7:
    print("usage: script aligned_fa best_tsv out_annot_fa out_sorted_fa out_prefix split_dbs", file=sys.stderr)
    sys.exit(2)

aligned_fa, best_tsv, out_fa, out_sorted_fa, out_prefix, split_dbs_arg = sys.argv[1:]
best = {}

# Read best hits with strand info
with open(best_tsv) as fh:
    for line in fh:
        line=line.rstrip("\n")
        if not line: 
            continue
        cols=line.split("\t")
        q=cols[0]
        s=cols[1] if len(cols)>1 else "None"
        b=cols[2] if len(cols)>2 else "0"
        strand=cols[3] if len(cols)>3 else "+"
        best[q]=(s,b,strand)

def clean(x):
    """Replace problematic characters that interfere with delimiter"""
    return x.replace(" ", "_").replace("|","_")

def extract_db_prefix(hit_name):
    """Extract database prefix like {DB1} from hit name"""
    match = re.match(r'\{([^}]+)\}', hit_name)
    return match.group(1) if match else None

# Parse split_dbs argument
split_dbs = []
if split_dbs_arg and split_dbs_arg != "":
    split_dbs = [db.strip() for db in split_dbs_arg.split(",")]

# Annotate aligned FASTA
annotated_count = 0
missing_count = 0
sequences = []  # Store for sorting
db_sequences = {db: [] for db in split_dbs}  # Store per-database sequences
unmatched_sequences = []  # Store sequences without hits (when split_dbs enabled)

with open(aligned_fa) as inf:
    header=None
    seq=[]
    
    for line in inf:
        if line.startswith(">"):
            # Save previous sequence if exists
            if header is not None:
                seq_str = "".join(seq)
                sequences.append((header, seq_str))
                
                # If splitting by database, categorize this sequence
                if split_dbs:
                    # Extract DB from header (it's in the hit name)
                    hit_part = header[1:].split("|")[0]  # Get first part after ">"
                    db_prefix = extract_db_prefix(hit_part)
                    if db_prefix and db_prefix in split_dbs:
                        db_sequences[db_prefix].append((header, seq_str))
                    else:
                        # No DB prefix found - this is unmatched (keeps original name)
                        unmatched_sequences.append((header, seq_str))
            
            # Process new header
            raw=line[1:].rstrip("\n")
            qid=raw.split()[0]
            
            if qid in best:
                s, b, strand = best[qid]
                annotated_count += 1
            else:
                s, b, strand = "None", "0", "+"
                missing_count += 1
                print(f"WARNING: No hit found for sequence: {qid}", file=sys.stderr)
            
            # For matched sequences, prepend hit name
            # For unmatched sequences, keep original name unchanged
            if s != "None":
                s_clean = clean(s)
                new_header = f">{s_clean}|{strand}|bs={b}|{raw}"
            else:
                # Keep original header for unmatched sequences
                new_header = f">{raw}"
            header=new_header
            seq=[]
        else:
            seq.append(line)
    
    # Save last sequence
    if header is not None:
        seq_str = "".join(seq)
        sequences.append((header, seq_str))
        
        if split_dbs:
            hit_part = header[1:].split("|")[0]
            db_prefix = extract_db_prefix(hit_part)
            if db_prefix and db_prefix in split_dbs:
                db_sequences[db_prefix].append((header, seq_str))
            else:
                unmatched_sequences.append((header, seq_str))

# Write unsorted version
with open(out_fa, "w") as out:
    for header, seq in sequences:
        out.write(header+"\n")
        out.write(seq)

# Write sorted version (by top hit name)
sequences_sorted = sorted(sequences, key=lambda x: x[0].split("|")[0].lower())
with open(out_sorted_fa, "w") as out:
    for header, seq in sequences_sorted:
        out.write(header+"\n")
        out.write(seq)

# Write per-database files if requested
if split_dbs:
    for db in split_dbs:
        if db_sequences[db]:
            db_fa = f"{out_prefix}.{db}.aligned.fasta"
            with open(db_fa, "w") as out:
                for header, seq in db_sequences[db]:
                    out.write(header+"\n")
                    out.write(seq)
            print(f"Wrote {len(db_sequences[db])} sequences to {db_fa}", file=sys.stderr)
        else:
            print(f"WARNING: No sequences found for database {db}", file=sys.stderr)
    
    # Write unmatched sequences file
    if unmatched_sequences:
        unmatched_fa = f"{out_prefix}.unmatched.aligned.fasta"
        with open(unmatched_fa, "w") as out:
            for header, seq in unmatched_sequences:
                out.write(header+"\n")
                out.write(seq)
        print(f"Wrote {len(unmatched_sequences)} unmatched sequences to {unmatched_fa}", file=sys.stderr)

print(f"Annotated {annotated_count} sequences, {missing_count} had no hits", file=sys.stderr)
PY

echo "Done.

OUTPUTS:
  $ANNOT_FA - annotated alignment (original order), $TOTAL_HITS sequences
  $ANNOT_FA_SORTED - annotated alignment (sorted by hit name), $TOTAL_HITS sequences
  $BEST_TSV - best hits table, $REAL_HITS matched ($(awk "BEGIN {printf \"%.1f\", ($REAL_HITS/$TOTAL_HITS)*100}")%)
"

# Report per-database outputs if created
if [ -n "$SPLIT_DBS" ]; then
    IFS=',' read -ra DBS <<< "$SPLIT_DBS"
    for db in "${DBS[@]}"; do
        db_file="${OUT_PREFIX}.${db}.aligned.fasta"
        if [ -f "$db_file" ]; then
            count=$(grep -c "^>" "$db_file" || true)
            echo "  $db_file - $db matches, $count sequences"
        fi
    done
    
    # Report unmatched file if it exists
    unmatched_file="${OUT_PREFIX}.unmatched.aligned.fasta"
    if [ -f "$unmatched_file" ]; then
        count=$(grep -c "^>" "$unmatched_file" || true)
        echo "  $unmatched_file - no database match, $count sequences"
    fi
    echo ""
fi

exit 0

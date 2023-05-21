#/bin/sh
T0=$(date +%s); T1=$(date +%s)
Secs_HMS() { echo "  done in $(( ${1} / 3600 ))h $(( (${1} / 60) % 60 ))m $(( ${1} % 60 ))s"; }
#
if [[ $# -eq 0 ]] ; then
    echo 'No arguments supplied!'
    exit 1
fi
Files=($(which "mafft") $(which "cons") $(which "seqkit") $(which "seqret") "$1" )
for f in "${Files[@]}" ; do
if [ ! -f "$f" ]; then
echo "File $f: not found"
exit 1
fi
done
if [ -z "$2" ]
        then BnkSz=50
        else BnkSz=$2
fi
#
echo Reordering bank and making consensus sequences
mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) --nuc --quiet --retree 0 --reorder $1 \
        | seqkit split2 -s $BnkSz -O ./ > /dev/null  2>&1 \
                && rename -f "s/stdin.part/${1%.*}/" stdin.part*.fasta \
                && rename -f "s/fasta$/bnk/" *.fasta \
        || exit 1
Secs_HMS "$(($(date +%s) - ${T1}))"; T1=$(date +%s)
#
echo Aligning consensus sequences
for b in $(find . -name "${1%.*}_*.bnk"); do
        mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) --nuc --reorder --quiet $b \
        | cons -plurality 18 -name $(basename "${b%.*}") -filter \
        | awk '!/^>/{gsub(/[Nn]/, "-")}1' > $(basename "${b%.*}").cons
done;
Secs_HMS "$(($(date +%s) - ${T1}))"; T1=$(date +%s)
#
echo Final alignment
cat ${1%.*}_*.cons > ${1%.*}.clw
#mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) --nuc --auto --reorder --quiet ${1%.*}.clw \
mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet ${1%.*}.clw \
        | seqret -filter -osformat2 msf > ${1%.*}.msf
Secs_HMS "$(($(date +%s) - ${T1}))";
#
echo Job completed
Secs_HMS "$(($(date +%s) - ${T0}))"

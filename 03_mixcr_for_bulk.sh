# Written by Etienne Becht
# Modified by Maxime Meylan
#
# get clonotype information from bulk RNAseq
#
input_dir="$1"
out_dir="$2"
mixcr_executable="$3"
files="$4"
lines=`cat $files | awk '{print $1}'`
echo $lines
for file in $lines; do
    echo "$file"
    R1="${input_dir}${file}_R1.fastq.gz"
    R2="${input_dir}${file}_R2.fastq.gz"
    #R1="${input_dir}${file}_RNA_R1_fastq.gz"
    #R2="${input_dir}${file}_RNA_R2_fastq.gz"
    echo "$R1"

    if test -f "$R1"; then
  out="${out_dir}/${file}"
  mkdir -p "$out"
  mkdir -p "$out/clones_fastq"
  cd $out

  $mixcr_executable -Xmx48g analyze shotgun \
    --species hs \
    --starting-material rna \
    --only-productive \
    --threads 24 \
    --align "-OsaveOriginalReads=true" \
    --contig-assembly \
    --impute-germline-on-export \
    --assemble "-ObadQualityThreshold=0" \
    "${R1}" "${R2}" "${file}repertoire"

    $mixcr_executable assemble \
                --threads 24 -f -a "${file}repertoire.vdjca" "${file}repertoire.clna"

  $mixcr_executable exportReadsForClones \
    -s "${file}repertoire.clna" "./clones_fastq/${file}_reads.fastq.gz"
  fi

done

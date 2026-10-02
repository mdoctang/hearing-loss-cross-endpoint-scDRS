#!/usr/bin/env bash
# chunked_download.sh URL OUT SIZE [NCHUNKS] [PAR]
# 分块并发下载大文件并拼接，带字节数校验。适用于支持 Range 请求的 HTTP(S) 源。
set -u
URL="$1"; OUT="$2"; SIZE="$3"; N="${4:-16}"; PAR="${5:-8}"
PDIR="${OUT}.parts"
export URL PDIR
mkdir -p "$PDIR"
rm -f "$PDIR"/part_*
CHUNK=$(( SIZE / N + 1 ))
for i in $(seq -w 0 $((N-1))); do
  s=$((10#$i * CHUNK))
  e=$(( (10#$i + 1) * CHUNK - 1 ))
  if [ "$e" -ge "$SIZE" ]; then e=$((SIZE-1)); fi
  echo "$i $s $e"
done | xargs -P "$PAR" -n 3 bash -c 'curl -sL --retry 5 --retry-delay 2 -r "$1-$2" -o "$PDIR/part_$0" "$URL"'
total=0
for f in "$PDIR"/part_*; do sz=$(stat -c %s "$f"); total=$((total+sz)); done
echo "chunk bytes total: $total / expected: $SIZE"
if [ "$total" -ne "$SIZE" ]; then echo "SIZE MISMATCH - aborting"; exit 1; fi
cat "$PDIR"/part_* > "$OUT"
rm -rf "$PDIR"
echo "assembled: $OUT ($(stat -c %s "$OUT") bytes)"

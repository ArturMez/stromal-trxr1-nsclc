#!/usr/bin/env bash
# Atlas stage 00: download every public input listed in atlas/sources.tsv and check it.
#
# WHAT: for every row of sources.tsv whose kind is "file" and whose "how" starts with "download",
#       fetches the url into <dest>/<local_path> (curl, resumable, three retries), then compares the
#       sha256 of the file with the sources.tsv column. The E-MTAB-6653 matrix is served as a zip
#       made on request (its bytes change with every fetch), so the zip is unpacked into
#       raw_counts/matrix_market/ and the three members are checked instead. Rows that are record
#       snapshots or local tables are skipped (they are not read by the build); the two E-MTAB-13526
#       author h5ad files (97 GiB) are never downloaded, the build reads their obs by HTTP range reads.
# FROM: atlas/sources.tsv.
# TO:   <dest>, by default raw_root/external_scRNA from config.json (the "scrna_raw" key overrides it),
#       laid out as the build expects; a report on stdout; exit code 1 when a file fails its check.
# WHY:  one command fetches and verifies every public deposit the atlas build reads (LuCA, the three
#       ArrayExpress deposits, the Single Cell Expression Atlas export, the STTT portal files, the
#       GEO series). Space: about 27 GB.
#
# RUN:  bash atlas/00_download.sh [--dest DIR] [--check-only] [--urls-only] [--only PREFIX]
#         --check-only   verify the checksums of files already on disk, download nothing
#         --urls-only    ranged one-byte GET of every url; report status and remote size
#         --only PREFIX  restrict to local paths starting with PREFIX (for example GSE325414/)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SOURCES="$HERE/sources.tsv"
DEST=""
MODE="download"
ONLY=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dest) DEST="$2"; shift 2 ;;
    --check-only) MODE="check"; shift ;;
    --urls-only) MODE="urls"; shift ;;
    --only) ONLY="$2"; shift 2 ;;
    *) echo "unknown argument $1"; exit 2 ;;
  esac
done
if [ -z "$DEST" ]; then
  DEST="$(cd "$HERE" && python3 -c 'import common; print(common.SCRNA)')" || exit 2
fi
echo "sources: $SOURCES"
echo "destination: $DEST   mode: $MODE"

sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
fails=0; done_n=0; skipped=0

# columns: 1 source_id 2 kind 3 arm 4 study 5 accession 6 repository 7 url 8 local_path 9 size_bytes
#          10 sha256 11 version_or_date 12 how 13 read_by_build 14 url_check 15 license
while IFS=$'\t' read -r sid kind arm study acc repo url path size sum ver how readb chk lic; do
  [ "$kind" = "file" ] || continue
  [ -z "$ONLY" ] || case "$path" in "$ONLY"*) ;; *) continue ;; esac
  case "$how" in
    download*) ;;
    extracted*) [ "$MODE" = "urls" ] && continue
                if [ -f "$DEST/$path" ]; then
                  got="$(sha "$DEST/$path")"
                  if [ "$got" = "$sum" ]; then echo "  ok    $path"; else echo "  FAIL  $path sha256 $got"; fails=$((fails+1)); fi
                else echo "  MISSING $path (unpack the zip first)"; fails=$((fails+1)); fi
                continue ;;
    *) skipped=$((skipped+1)); continue ;;
  esac
  if [ "$MODE" = "urls" ]; then
    hdr="$(curl -sS -L -m 60 -r 0-0 -o /dev/null -D - "$url" 2>&1)"  # discard the body; keep the headers
    code="$(printf '%s\n' "$hdr" | awk 'toupper($1) ~ /^HTTP/ {c=$2} END {print c}')"
    rsz="$(printf '%s\n' "$hdr" | awk -F/ 'tolower($0) ~ /^content-range/ {gsub(/\r/,"",$2); print $2}')"
    if [ "$code" = "200" ] || [ "$code" = "206" ]; then  # 200 = full answer, 206 = partial content (ranged request)
      echo "  ok    $code ${rsz:-no size} $path"
      [ -n "$rsz" ] && [ "$rsz" != "$size" ] && { echo "        remote size $rsz differs from $size"; fails=$((fails+1)); }
    else echo "  FAIL  ${code:-no answer} $url"; fails=$((fails+1)); fi
    done_n=$((done_n+1))
    continue
  fi
  if [ "$MODE" = "download" ]; then
    mkdir -p "$DEST/$(dirname "$path")"
    have="no"
    if [ -f "$DEST/$path" ]; then
      case "$path" in
        *.zip) [ "$(wc -c < "$DEST/$path" | tr -d ' ')" = "$size" ] && have="yes" ;;
        *) [ "$(sha "$DEST/$path")" = "$sum" ] && have="yes" ;;
      esac
    fi
    if [ "$have" = "yes" ]; then
      echo "  have  $path"
    else
      echo "  get   $path  <- $url"
      curl -sS -L -C - --retry 3 --retry-delay 5 --max-time 36000 -o "$DEST/$path" "$url" \
        || { echo "  FAIL  download $path"; fails=$((fails+1)); continue; }
    fi
  fi
  if [ ! -f "$DEST/$path" ]; then echo "  MISSING $path"; fails=$((fails+1)); continue; fi
  case "$path" in
    E-MTAB-6653/*.zip|E-MTAB-6653/*/*.zip)
      # made on request; the bytes differ from fetch to fetch, so only the size is compared here
      got="$(wc -c < "$DEST/$path" | tr -d ' ')"
      if [ "$got" = "$size" ]; then echo "  ok    $path (size; zip made on request)"; else echo "  FAIL  $path size $got"; fails=$((fails+1)); fi
      if [ "$path" = "E-MTAB-6653/raw_counts/E-MTAB-6653_raw_counts_mtx.zip" ] && [ "$MODE" = "download" ]; then
        mkdir -p "$DEST/E-MTAB-6653/raw_counts/matrix_market"
        unzip -o -q "$DEST/$path" -d "$DEST/E-MTAB-6653/raw_counts/matrix_market"
      fi ;;
    *)
      got="$(sha "$DEST/$path")"
      if [ "$got" = "$sum" ]; then echo "  ok    $path"; else echo "  FAIL  $path sha256 $got, expected $sum"; fails=$((fails+1)); fi ;;
  esac
  done_n=$((done_n+1))
done < <(tail -n +2 "$SOURCES")

echo "files handled: $done_n; rows skipped (records, local tables, range-read, constituents): $skipped; failures: $fails"
[ "$fails" -eq 0 ]

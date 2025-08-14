!/bin/bash
set -euo pipefail

if [ "$#" -eq 0 ]; then
  echo "Usage:"
  echo "  $0 /path/to/files/*.gz            # pass a glob (UNQUOTED)"
  echo "  $0 file1.gz file2.gz ...          # or list files explicitly"
  echo "  $0 /path/to/dir                   # or a directory (processes *.gz inside)"
  exit 1
fi

process_file() {
  local gz="$1"
  [[ -f "$gz" ]] || { echo "Skip (not a file): $gz"; return; }

  echo "Processing: $gz"
  gunzip -f -- "$gz"                           # decompress; leaves file without .gz
  local decompressed="${gz%.gz}"
  local new_name="${decompressed/_R2/_R1}"     # rename _R2 -> _R1
  if [[ "$decompressed" != "$new_name" ]]; then
    mv -- "$decompressed" "$new_name"
  fi
  gzip -f -- "$new_name"                        # recompress
}

# If a single directory is given, process all *.gz inside it
if [[ "$#" -eq 1 && -d "$1" ]]; then
  while IFS= read -r -d '' f; do
    process_file "$f"
  done < <(find "$1" -type f -name '*.gz' -print0)
else
  # Otherwise, treat all arguments as files (supports glob expansion)
  for f in "$@"; do
    process_file "$f"
  done
fi

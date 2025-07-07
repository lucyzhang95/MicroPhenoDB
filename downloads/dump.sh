# dump.sh: download all URLs passed as arguments
set -euo pipefail

urls=(
  "http://lilab2.sysu.edu.cn/Microbe/download/core_table.txt"
  "http://lilab2.sysu.edu.cn/Microbe/download/EFO.txt"
  "http://lilab2.sysu.edu.cn/Microbe/download/NCIT.txt"
)
for url in "${urls[@]}"; do
  echo "Downloading $url …"

  wget \
    --progress=bar:force \
    --tries=3 \
    --timeout=30 \
    --waitretry=5 \
    --continue \
    --user-agent="Mozilla/5.0 (Macintosh; Intel Mac OS X 13_0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36" \
    "$url"

  echo
done

echo "All MicroPhenoDB downloads are complete."

# Usage:
# change the directory to where the script is located
# chmod +x dump.sh
# ./dump.sh
#!/usr/bin/env bash
set -euo pipefail

# clean_src.sh
# Purpose: Remove compiled artifacts from `src/` so the directory is ready
#          for a fresh compilation.
# Usage: ./clean_src.sh [--yes]
# Options:
#  --yes   : do not prompt for confirmation; remove files immediately.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC_DIR="$ROOT_DIR/src"

echo "[clean_src] Target directory: $SRC_DIR"

cd "$SRC_DIR"

# Show common compiled artifact patterns (if present)
echo "[clean_src] Found the following artifacts (if any):"
ls -1 -- main *.o *.mod *.mod~ 2>/dev/null || true

if [ "${1:-}" != "--yes" ]; then
  read -p "Proceed to remove the listed artifacts? [y/N] " ans
  case "$ans" in
    [Yy]*) ;;
    *) echo "[clean_src] Aborted by user."; exit 0 ;;
  esac
fi

echo "[clean_src] Removing artifacts..."
rm -f -- main *.o *.mod *.mod~ || true

echo "[clean_src] Clean complete. $SRC_DIR is ready for recompilation."

exit 0

#!/usr/bin/env bash
# Validates an integration-mapping-pipeline workdir configuration.
#
# Usage:
#   bash validate.sh <workdir>
#
# Prints VALID (exit 0) or INVALID with error messages (exit 1).
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <workdir>" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${SCRIPT_DIR}/validate.py" "$1"

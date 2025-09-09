#!/usr/bin/env bash
set -e

HOST="xxx"
USER="xxx"                    # your GEO user
PASS="xxx"         # your GEO password/token

# Make sure this is the directory you actually want to push
LOCAL_DIR="/Users/djuna/Documents/abca7_dryad/counts"  
REMOTE_DIR="uploads/xxx"   # adjust to your GEO folder

lftp -u "$USER","$PASS" "$HOST" <<EOF
  set net:timeout 30
  set ftp:passive-mode true

  # -R        = reverse (upload local → remote)
  # -c        = continue partially transferred files
  # -P 4      = use 4 parallel transfers
  # --verbose = show each file as it’s transferred
  mirror -R -c --verbose "$LOCAL_DIR" "$REMOTE_DIR"

  quit
EOF

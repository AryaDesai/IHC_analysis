#!/bin/sh
# Copies code, raw data and results between this computer and the lab server.
#
#   ./sync.sh pull    get the latest code, raw data and results
#   ./sync.sh push    send this computer's results to the server
#
# Add --dry-run to either one to list what would be copied without copying.
# The raw data is mirrored from the server, so local raw files that are not on
# the server are deleted. Results are copied in both directions without deleting
# anything, and a file is replaced only when the incoming copy is newer.

set -e
SERVER="${SYNC_SERVER:-/Volumes/Arya CDB/MEK_oscillation_project}"
RESULTS=analysis_scripts/ihc_results
HERE=$(cd "$(dirname "$0")" && pwd -P)
LOCAL=$(dirname "$HERE")

# The rsync from Homebrew is used because the one that comes with macOS
# (openrsync) cannot show the overall progress of a copy.
RSYNC=/opt/homebrew/bin/rsync
if [ ! -x "$RSYNC" ]; then
    echo "Homebrew rsync was not found at $RSYNC. Install it with: brew install rsync"
    exit 1
fi
if [ ! -d "$SERVER" ]; then
    echo "Cannot reach the server folder: $SERVER (is the share mounted?)"
    exit 1
fi
if [ "$LOCAL" = "$(cd "$SERVER" && pwd -P)" ]; then
    echo "This copy is on the server. Run sync from a copy on this computer."
    exit 1
fi

DRY=""
[ "$2" = "--dry-run" ] && DRY="--dry-run --itemize-changes"
# -rt copies folders and keeps modification times. Owners and permissions are
# left out because the SMB share does not keep them. --info=progress2 prints one
# line with the amount copied, the speed and the percentage done, and
# --no-inc-recursive lists all files before copying so that the percentage is
# out of the full total.
OPTS="-rt --info=progress2 --exclude=.DS_Store $DRY"
SKIP="--exclude=/models --exclude=/runtime"

case "$1" in
pull)
    [ -z "$DRY" ] && git -C "$HERE" pull --ff-only
    echo "Raw data:"
    "$RSYNC" $OPTS --delete --exclude=/analysis_scripts --exclude='/analysis_scripts_backup_*' \
        --exclude=.venv --exclude='*.ims' "$SERVER/" "$LOCAL/"
    echo "Results:"
    mkdir -p "$LOCAL/$RESULTS"
    "$RSYNC" $OPTS --update $SKIP "$SERVER/$RESULTS/" "$LOCAL/$RESULTS/"
    ;;
push)
    echo "Results:"
    "$RSYNC" $OPTS --update $SKIP "$LOCAL/$RESULTS/" "$SERVER/$RESULTS/"
    echo "Code changes not yet committed:"
    git -C "$HERE" status --short
    ;;
*)
    echo "Usage: ./sync.sh pull | ./sync.sh push   (add --dry-run to only list changes)"
    exit 1
    ;;
esac

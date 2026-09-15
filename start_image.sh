#!/bin/bash
# Usage: ./start_image.sh [host_data_dir]
# The image bakes the whole repo into /data (see Dockerfile's `ADD . $HOME/data/`),
# so mounting a host folder over all of /data would hide the APE-Gen code itself.
# Instead this mounts host_data_dir (default: ./data next to this script) onto
# /data/intermediate_files, which is where results land by default (--dir).
# Results survive after the container exits, instead of being lost with --rm.
# If you pass --dir to New_APE-Gen.py, point it under intermediate_files/ (e.g.
# --dir intermediate_files/batch1) so it also lands on the host.
HOST_DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data}"
mkdir -p "$HOST_DIR"
docker run -it --rm -v "$HOST_DIR":/data/intermediate_files kavrakilab/apegen2.0:rc1

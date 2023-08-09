#!/bin/bash

while /home/gridsan/djuna/inotify-tools-3.14/src/inotifywait -r /home/gridsan/djuna/homer/github/ABCA7lof2/pdf_figures --format true; do
    echo "syncing to"
    rclone sync --include='*Figure*' /home/gridsan/djuna/homer/github/ABCA7lof2/pdf_figures/. dropbox:Apps/Overleaf/ABCA7lof2/pdf_figures/
done


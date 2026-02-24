md5sum -c all.md5
for f in *.gz; do tar -zxvf "$f"; done

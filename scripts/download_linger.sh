#!/bin/bash
Datadir=/home/jiehoonk/project/crisprscreen/data/adj
cd $Datadir
wget --load-cookies /tmp/cookies.txt "https://drive.usercontent.google.com/download?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.usercontent.google.com/download?id=1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b'  -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b" -O data_bulk.tar.gz && rm -rf /tmp/cookies.txt
#!/usr/bin/env bash

url_list="mouse_reference_urls.txt"
url_list=$1
echo "url_list: $url_list"

for u in `cat $url_list`; do echo "processing $u" && wget $u; done

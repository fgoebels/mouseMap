#!/bin/bash

filename="$1"
while read -r line
do
    name=$line
    echo "Downloading: $line"
    curl -O $name
done < "$filename"

#!/bin/bash


: ${num_cores:=1}

while getopts ":c:" opt; do
  case $opt in
    c) num_cores="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

mkdir /data/pfamscan_results
mkdir /data/final_results
snakemake -s /usr/local/bin/snakefile -j $num_cores



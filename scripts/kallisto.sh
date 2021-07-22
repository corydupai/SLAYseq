#!/bin/bash
while getopts i:l:r: flag
do
    case "${flag}" in
        i) index=${OPTARG};;
        l) read1=${OPTARG};;
        r) read2=${OPTARG};;
    esac
done
echo "Index: $index";
echo "Read1: $read1";
echo "Read2: $read2";
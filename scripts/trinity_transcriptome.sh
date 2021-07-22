#!/bin/bash
while getopts c:d: flag
do
    case "${flag}" in
        c) contrasts=${OPTARG};;
        d) directory=${OPTARG};;
    esac
done
#!/usr/bin/env bash

bam_file=$1
read_length=$(samtools view $bam_file | perl -lane 'print scalar(split //,$F[9]) and last if $F[5]=~/^[\dM]*$/;')
echo $read_length

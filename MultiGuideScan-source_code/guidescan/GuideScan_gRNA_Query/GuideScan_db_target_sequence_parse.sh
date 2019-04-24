#!/bin/bash

set -e

samtools=""
db=""
outfile=""

`$samtools view $db | awk '{print $1}' - | cut -c1-20 > $outfile`

#!/usr/bin/env bash
CODE='BEGIN {
    print "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT"
}

/^#/ {
    next
}

{
    cov = $7 + $8
    strand = ($6 == "+") ? "F" : "R"
    pct1 = (cov == 0) ? "0" : sprintf("%.2f", $8 * 100 / cov)
    pct2 = (cov == 0) ? "0" : sprintf("%.2f", $7 * 100 / cov)

    print $1 ":" $3 "\t" $1 "\t" $3 "\t" strand "\t" ($7 + $8) "\t" pct1 "\t" pct2
}'

if [ -z "$1" ]; then
	awk "$CODE"
elif [ -f "$1" ]; then
	if [ -z `file "$1" | grep -F gzip` ]; then
		cat "$1" | awk "$CODE"
	else
		gunzip -c "$1" | awk "$CODE"
	fi
else
	echo "Not a valid file: $1" 2>&1
fi
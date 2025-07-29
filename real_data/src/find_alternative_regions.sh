#!/bin/bash

module load impg

paf_file="$1"
bed_in="$2"
bed_out="$3"

> "$bed_out"

while read -r chr start end name_rest; do
    name=${name_rest:-}

    region="${chr}\t${start}\t${end}"

    impg_regions=$(impg query -f \
        --paf-files "$paf_file" \
        -b <(echo -e "$region") | cut -f1-3 | \
        awk -v chr="$chr" -v start="$start" -v end="$end" -F'\t' '
        {
            key = $1 ":" $2 "-" $3;
            query_key = chr ":" start "-" end;
            if (key == query_key) next;

            name = $1;
            gsub(/:[^:]*-[0-9]+$/, "", name);
            printf("%s%s:%s-%s", (NR==1 ? "" : ","), name, $2, $3);
        } END { print "" }')

    if [ -z "$impg_regions" ]; then
        if [ -z "$name" ]; then
            echo -e "$chr\t$start\t$end\tunknown" >> "$bed_out"
        else
            echo -e "$chr\t$start\t$end\t$name" >> "$bed_out"
        fi
    else
        if [ -z "$name" ]; then
            echo -e "$chr\t$start\t$end\tunknown\t$impg_regions" >> "$bed_out"
        else
            echo -e "$chr\t$start\t$end\t$name\t$impg_regions" >> "$bed_out"
        fi
    fi
done < "$bed_in"


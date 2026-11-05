root=$(pwd)
out="$root/undesired_rxns_by_enzyme.tsv"

{
    printf 'directory\tenzyme\treactions\n'

    find "$root" -type f -name 'runRBA.undesired_rxns_used.txt' -print0 |
    while IFS= read -r -d '' file; do
        dir=$(basename "$(dirname "$file")")

        awk -F '\t' -v OFS='\t' -v dir="$dir" '
        {
            field = $1
            gsub(/^"|"$/, "", field)

            # Save ordinary reaction names for matching ENZLOAD rows.
            if (field !~ /^(ENZSYN-|ENZLOAD-|PTM-)/ &&
                field != "rxn") {
                rxn[field] = 1
            }

            # ENZLOAD rows have the form:
            # ENZLOAD-<reaction>-<enzyme>
            if (field ~ /^ENZLOAD-/) {
                load = substr(field, length("ENZLOAD-") + 1)
                best_rxn = ""
                best_prefix = ""

                # Choose the longest known reaction prefix.
                for (r in rxn) {
                    prefix = r
                    sub(/^RXN-/, "", prefix)
                    prefix = prefix "-"
                    if (index(load, prefix) == 1 &&
                        length(prefix) > length(best_prefix)) {
                        best_rxn = r
                        best_prefix = prefix
                    }
                }

                if (best_rxn != "") {
                    reaction_part = best_rxn
                    sub(/^RXN-/, "", reaction_part)
                    enzyme = substr(load, length(reaction_part) + 2)
                    key = enzyme SUBSEP best_rxn

                    if (!(key in seen)) {
                        seen[key] = 1
                        if (enzyme in reactions) {
                            reactions[enzyme] = reactions[enzyme] ", " reaction_part
                        } else {
                            reactions[enzyme] = reaction_part
                            enzyme_order[++n] = enzyme
                        }
                    }
                }
            }
        }

        END {
            for (i = 1; i <= n; i++) {
                enzyme = enzyme_order[i]
                print (i == 1 ? dir : ""), enzyme, reactions[enzyme]
            }
        }
        ' "$file"
    done
} > "$out"

# run python3 singleTermPipeline.py for each term in the file ../data/go_terms_test_candidates_maxlen500_minmembers100_ordered.tsv
# Usage: ./runPipeline.sh

IFS=$'\t'

# read ../data/go_terms_test_candidates_maxlen500_minmembers100_ordered.tsv and run singleTermPipeline.py for each term
while read line; do
    # extract the first field of the line
    term=$(echo $line | cut -f1)
    # run singleTermPipeline.py for each term
    python3 singleTermPipeline.py $term
done < ../data/go_terms_test_candidates_maxlen500_minmembers100_ordered.tsv



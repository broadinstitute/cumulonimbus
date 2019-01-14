import csv
import sys

sort_header = sys.argv[1]
in_file_name = sys.argv[2]

with open(in_file_name, "r", newline="") as in_file:
    in_reader = csv.reader(in_file, delimiter='\t')
    header_row = next(in_reader)
    sort_header_index = header_row.index(sort_header)
    script_file_name = "header_sort_script.sh"
    with open(script_file_name, "w") as script_file:
        script_content = "(head -n 1 %s && tail -n +2 %s | sort -b -k %s)" \
                         %(in_file_name, in_file_name, sort_header_index)
        script_file.write(script_content)


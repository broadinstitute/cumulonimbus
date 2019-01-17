import csv
import sys

file_name = sys.argv[1]
numerical_headers = sys.argv[2:]

with open(file_name, "r", newline="") as file:
    reader = csv.reader(file, delimiter='\t')
    header_row = next(reader)
    numerical_cols = list(map(lambda header: header_row.index(header), numerical_headers))
    for row in reader:
        for col in numerical_cols:
            value = float(row[col])

import csv
import sys

file_name = sys.argv[1]
header1 = sys.argv[2]
header2 = sys.argv[3]
epsilon = float(sys.argv[4])

with open(file_name, "r", newline="") as file:
    reader = csv.reader(file, delimiter='\t')
    header_row = next(reader)
    print("\t".join(header_row))
    col1 = header_row.index(header1)
    col2 = header_row.index(header2)
    for row in reader:
        valueString1 = row[col1]
        valueString2 = row[col2]
        if row[col1] != row[col2]:
            value1 = float(valueString1)
            value2 = float(valueString2)
            denom = (abs(value1) + abs(value2))/2
            if denom != 0 and abs(value1 - value2)/denom > epsilon:
                print("\t".join(row))



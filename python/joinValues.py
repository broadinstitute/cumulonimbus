import csv
import sys

file_name1 = sys.argv[1]
id_header1 = sys.argv[2]
value_header1 = sys.argv[3]
file_name2 = sys.argv[4]
id_header2 = sys.argv[5]
value_header2 = sys.argv[6]
file_name_out = sys.argv[7]
id_header_out = sys.argv[8]
value_header1_out = sys.argv[9]
value_header2_out = sys.argv[10]

with open(file_name1, "r", newline="") as file1, \
        open(file_name2, "r", newline="") as file2, \
        open(file_name_out, "w", newline="") as file_out:
    reader1 = csv.reader(file1, delimiter='\t')
    reader2 = csv.reader(file2, delimiter='\t')
    writer = csv.writer(file_out, delimiter='\t')
    header_row1 = next(reader1)
    id_index1 = header_row1.index(id_header1)
    value_index1 = header_row1.index(value_header1)
    header_row2 = next(reader2)
    id_index2 = header_row2.index(id_header2)
    value_index2 = header_row2.index(value_header2)

    header_row_out = [id_header_out, value_header1_out, value_header2_out]
    writer.writerow(header_row_out)

    def next_or_none(iterator):
        try:
            result = next(iterator)
        except StopIteration:
            result = None
        return result


    row1 = next_or_none(reader1)
    row2 = next_or_none(reader2)
    while row1 is not None and row2 is not None:
        id1 = row1[id_index1]
        id2 = row2[id_index2]
        if id1 == id2:
            value1 = row1[value_index1]
            value2 = row2[value_index2]
            row_out = [id1, value1, value2]
            writer.writerow(row_out)
            row1 = next_or_none(reader1)
            row2 = next_or_none(reader2)
        else:
            if id1 < id2:
                row1 = next_or_none(reader1)
            else:
                row2 = next_or_none(reader2)



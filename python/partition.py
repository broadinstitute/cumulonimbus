import sys
import csv

class ArgumentException(Exception):
    def __init__(self, message):
        self.message = message

if len(sys.argv) != 6:
    raise ArgumentException("Got %s arguments, but need 5." % (len(sys.argv) - 1))

cutoff = float(sys.argv[1])
key = sys.argv[2]
in_file_name = sys.argv[3]
out_file_rare_name = sys.argv[4]
out_file_common_name = sys.argv[5]

with open(in_file_name, 'r', newline='') as in_file, \
        open(out_file_rare_name, 'w') as out_file_rare, \
        open(out_file_common_name, 'w') as out_file_common:
    in_reader = csv.reader(in_file, delimiter='\t')
    rare_writer = csv.writer(out_file_rare, delimiter='\t')
    common_writer = csv.writer(out_file_common, delimiter='\t')
    header_row = next(in_reader)
    key_index = header_row.index(key)
    rare_writer.writerow(header_row)
    common_writer.writerow(header_row)
    for row in in_reader:
        frequency = float(row[key_index])
        if frequency > 0.5:
            frequency = 1.0 - frequency
        if frequency < cutoff:
            rare_writer.writerow(row)
        else:
            common_writer.writerow(row)



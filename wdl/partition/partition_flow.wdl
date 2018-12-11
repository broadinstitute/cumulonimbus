version 1.0

workflow partition_flow {
    input {
        Float cutoff_frequency
        String frequency_column
        File variants
        String variants_rare_name
        String variants_common_name
    }

    call partition {
        input:
            cutoff_frequency = cutoff_frequency,
            frequency_column = frequency_column,
            variants = variants,
            variants_rare_name = variants_rare_name,
            variants_common_name  = variants_common_name
    }
}

task partition {
    input {
        Float cutoff_frequency
        String frequency_column
        File variants
        String variants_rare_name
        String variants_common_name
    }
    runtime {
        docker: "us.gcr.io/broad-gdr-dig-storage/metal-python:2018-08-28"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 20 HDD"
    }
    command <<<
        set -e
        python3 --version
        cat << EOF >partition.py
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
        EOF
        echo "=== BEGIN partition.py ==="
        cat partition.py
        echo "=== END partition.py ==="
        python3 partition.py ~{cutoff_frequency} ~{frequency_column} ~{variants} \
          ~{variants_rare_name} ~{variants_common_name}
    >>>
    output {
        File variants_rare = variants_rare_name
        File variants_common = variants_common_name
    }
}
version 1.0

struct InputFile {
    File file
    String dataset_name
    Map[String, String]? column_mapping
}

struct Ancestry {
    String name
    Array[InputFile] input_files
    Map[String, String]? column_mapping
}

workflow bottom_line {
    input {
        Array[Ancestry] ancestries
        Float cutoff_frequency
        String frequency_column = "FREQ"
        String marker_column = "MARKER"
        String size_column = "N"
        String output_prefix
        String output_suffix
        String scheme
        String std_err = "STDERR"
        String effect = "EFFECT"
        String p_value = "PVALUE"
        String alt_allele = "ALLELE1"
        String ref_allele = "ALLELE2"
        Map[String, String]? column_mapping = { "blub": "blub" }
    }

    scatter(ancestry in ancestries) {
        scatter(input_file in ancestry.input_files) {
            String ancestry_dataset_file_name = output_prefix + "_" + ancestry.name + "_" + input_file.dataset_name
            call partition {
                input:
                    file_column_mapping = input_file.column_mapping,
                    ancestry_column_mapping = ancestry.column_mapping,
                    global_column_mapping = column_mapping,
                    cutoff_frequency = cutoff_frequency,
                    frequency_column = frequency_column,
                    input_file = input_file.file,
                    variants_rare_name = ancestry_dataset_file_name + "_rare." + output_suffix,
                    variants_common_name  = ancestry_dataset_file_name + "_common." + output_suffix
            }
        }
        String ancestry_file_name = output_prefix + "_" + ancestry.name
        call pick_largest {
            input:
                input_files = partition.variants_rare,
                marker_column = marker_column,
                size_column = size_column,
                output_file_name = ancestry_file_name + "_rare." + output_suffix
        }
        call metal as metal_common_per_ancestry {
            input:
                input_files = partition.variants_common,
                column_counting = "LENIENT",
                overlap = true,
                marker = marker_column,
                weight_column = size_column,
                out_prefix = output_prefix,
                out_postfix = output_suffix,
                scheme = scheme,
                average_freq = true,
                min_max_freq = true,
                std_err = std_err,
                effect = effect,
                p_value = p_value,
                alt_allele = alt_allele,
                ref_allele = ref_allele
        }
    }
    call metal as metal_all_rare {
        input:
            input_files = pick_largest.output_file,
            column_counting = "LENIENT",
            overlap = false,
            marker = marker_column,
            weight_column = size_column,
            out_prefix = output_prefix,
            out_postfix = output_suffix,
            scheme = scheme,
            average_freq = true,
            min_max_freq = true,
            std_err = std_err,
            effect = effect,
            p_value = p_value,
            alt_allele = alt_allele,
            ref_allele = ref_allele
    }
    call metal as metal_all_common {
        input:
            input_files = metal_common_per_ancestry.out,
            column_counting = "LENIENT",
            overlap = false,
            marker = marker_column,
            weight_column = size_column,
            out_prefix = output_prefix,
            out_postfix = output_suffix,
            scheme = scheme,
            average_freq = true,
            min_max_freq = true,
            std_err = std_err,
            effect = effect,
            p_value = p_value,
            alt_allele = alt_allele,
            ref_allele = ref_allele
    }
    call concat {
        input:
            input_files = [metal_all_common.out, metal_all_rare.out],
            output_file_name = output_prefix + "." + output_suffix
    }
}

task partition {
    input {
        Map[String, String] file_column_mapping = { "blub": "blub" }
        Map[String, String] ancestry_column_mapping = { "blub": "blub" }
        Map[String, String] global_column_mapping = { "blub": "blub" }
        Float cutoff_frequency
        String frequency_column
        File input_file
        String variants_rare_name
        String variants_common_name
    }
    File settings_file = write_json(
        object {
            file_column_mapping: file_column_mapping,
            ancestry_column_mapping: ancestry_column_mapping,
            global_column_mapping: global_column_mapping,
            cutoff_frequency: cutoff_frequency,
            frequency_column: frequency_column,
            variants_rare_name: variants_rare_name,
            variants_common_name: variants_common_name
        }
    )
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
        import json
        import csv
        settingsFile = open("~{settings_file}", "r")
        settings = json.load(settingsFile)
        settingsFile.close()
        print("=== BEGIN settings ===")
        print(json.dumps(settings, sort_keys=True, indent=4))
        print("=== END settings ===")
        file_column_mapping = settings["file_column_mapping"]
        ancestry_column_mapping = settings["ancestry_column_mapping"]
        global_column_mapping = settings["global_column_mapping"]
        cutoff = float(settings["cutoff_frequency"])
        key = settings["frequency_column"]
        in_file_name = "~{input_file}"
        out_file_rare_name = settings["variants_rare_name"]
        out_file_common_name = settings["variants_common_name"]

        def map_column_name(name):
            name = file_column_mapping.get(name) or name
            name = ancestry_column_mapping.get(name) or name
            name = global_column_mapping.get(name) or name
            return name

        with open(in_file_name, 'r', newline='') as in_file, \
                open(out_file_rare_name, 'w') as out_file_rare, \
                open(out_file_common_name, 'w') as out_file_common:
            in_reader = csv.reader(in_file, delimiter='\t')
            rare_writer = csv.writer(out_file_rare, delimiter='\t')
            common_writer = csv.writer(out_file_common, delimiter='\t')
            header_row = list(map(map_column_name, next(in_reader)))
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
        python3 partition.py
    >>>
    output {
        File variants_rare = variants_rare_name
        File variants_common = variants_common_name
    }
}

task pick_largest {
    input {
        Array[File] input_files
        String marker_column
        String size_column
        String output_file_name
    }
    File settings_file = write_json(
        object {
            marker_column: marker_column,
            size_column: size_column,
            output_file: output_file_name
        }
    )
    runtime {
        docker: "us.gcr.io/broad-gdr-dig-storage/metal-python:2018-08-28"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 20 HDD"
    }
    command <<<
        set -e
        python3 --version
        cat << EOF > union.py
        import json
        import csv
        settingsFile = open("~{settings_file}", "r")
        settings = json.load(settingsFile)
        settingsFile.close()
        print("=== BEGIN settings ===")
        print(json.dumps(settings, sort_keys=True, indent=4))
        print("=== END settings ===")
        in_file_names = ["~{sep='\", \"' input_files}"]
        marker_column = settings["marker_column"]
        size_column = settings["size_column"]
        out_file_name = settings["output_file"]
        marker_list = []
        markers = set(marker_list)
        column_list = []
        union_data = {}
        for in_file_name in in_file_names:
            with open(in_file_name, 'r', newline='') as in_file:
                in_reader = csv.reader(in_file, delimiter='\t')
                header_row = next(in_reader)
                for column in header_row:
                    if not column in column_list:
                        column_list.append(column)
                marker_column_index = header_row.index(marker_column)
                size_column_index = header_row.index(size_column)
                for row in in_reader:
                    marker = row[marker_column_index]
                    size = float(row[marker_column_index])
                    if(not marker in markers):
                        marker_list.append(marker)
                        markers.add(marker)
                    row_dict = dict(zip(header_row, row))
                    union_entry = union_data[marker]
                    if(union_entry is None):
                        union_data[marker] = row_dict
                    else:
                        union_size = float(union_entry[size_column])
                        if(size > union_size):
                            union_data[marker] = row_dict
        with open(out_file_name, 'w') as out_file:
            out_writer = csv.writer(out_file, delimiter='\t')
            out_writer.writerow(column_list)
            for marker in marker_list:
                entry = union_data[marker]
                row = list(map(lambda column: entry[column] or "NA", column_list))
                out_writer.writerow(row)
        EOF
        echo "=== BEGIN union.py ==="
        cat union.py
        echo "=== END union.py ==="
        python3 union.py
    >>>
    output {
        File output_file = output_file_name
    }
}

task metal {
    input {
        Array[File] input_files
        String column_counting = "STRICT"
        Boolean overlap = false
        String marker = "MARKER"
        String weight_column = "N"
        String out_prefix
        String out_postfix
        String scheme = "scheme"
        Boolean average_freq = false
        Boolean min_max_freq = false
        String std_err = "STDERR"
        String effect = "EFFECT"
        String p_value = "PVALUE"
        String freq = "FREQ"
        String alt_allele = "ALLELE1"
        String ref_allele = "ALLELE2"
        Array[String] custom_variables = []
    }
    File settings_file = write_json(
        object {
            overlap: overlap,
            column_counting: column_counting,
            marker: marker,
            out_prefix: out_prefix,
            out_postfix: out_postfix,
            scheme: scheme,
            average_freq: average_freq,
            min_max_freq: min_max_freq,
            std_err: std_err,
            effect: effect,
            custom_variables: custom_variables,
            p_value: p_value,
            freq: freq,
            alt_allele: alt_allele,
            ref_allele: ref_allele
        }
    )
    runtime {
        docker: "us.gcr.io/broad-gdr-dig-storage/metal-python:2018-08-28"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 20 HDD"
    }
    command <<<
        set -e
        python3 --version
        cat << EOF >metalcast.py
        import json
        settingsFile = open("~{settings_file}", "r")
        settings = json.load(settingsFile)
        settingsFile.close()
        print("=== BEGIN settings ===")
        print(json.dumps(settings, sort_keys=True, indent=4))
        print("=== END settings ===")
        lines = []
        def addLine(lines, line):
            lines.append(line)
        def addValue(lines, prefix, dict, key):
            value = dict.get(key)
            if(value is not None):
                lines.append(prefix + " " + value)
        def addFlag(lines, prefix, dict, key):
            value = dict.get(key)
            if(value is not None):
                if(value):
                    lines.append(prefix + " ON")
                else:
                    lines.append(prefix + " OFF")
        def addTwoValues(lines, prefix, dict, key1, key2):
            value1 = dict.get(key1) or defaultValue1
            value2 = dict.get(key2) or defaultValue2
            if(value1 is not None and value2 is not None):
                lines.append(prefix + " " + value1 + " " + value2)
        def addCustomVariables(lines, dict, key):
            variables = dict.get(key)
            if(variables is not None):
                for variable in variables:
                    lines.append("CUSTOMVARIABLE " + variable)
                    lines.append("LABEL " + variable + " AS " + variable)
        addLine(lines, "SEPARATOR TAB")
        addValue(lines, "SCHEME", settings, "scheme")
        addValue(lines, "WEIGHTLABE", settings, "weight_column")
        addFlag(lines, "OVERLAP", settings, "overlap")
        addFlag(lines, "AVERAGEFREQ", settings, "average_freq")
        addFlag(lines, "MINMAXFREQ", settings, "min_max_freq")
        addValue(lines, "COLUMNCOUNTING", settings, "column_counting")
        addValue(lines, "STDERR", settings, "std_err")
        addValue(lines, "EFFECT", settings, "effect")
        addValue(lines, "MARKER", settings, "marker")
        addValue(lines, "PVALUE", settings, "p_value")
        addValue(lines, "FREQ", settings, "freq")
        addTwoValues(lines, "ALLELE", settings, "alt_allele", "ref_allele")
        addCustomVariables(lines, settings, "custom_variables")
        for input_file in ["~{sep='\", \"' input_files}"]:
            addLine(lines, "PROCESS " + input_file)
        addTwoValues(lines, "OUTFILE", globalSettings, "out_prefix", "out_postfix")
        addLine(lines, "ANALYZE")
        addLine("QUIT")
        script = reduce(lambda line1, line2: line1 + "\n" + line2, lines) + "\n"
        scriptFile = open("script.metal", "w")
        scriptFile.write(script)
        scriptFile.close()
        EOF
        echo "=== BEGIN metalcast.py ==="
        cat metalcast.py
        echo "=== END metalcast.py ==="
        python3 metalcast.py
        echo "=== BEGIN script.metal ==="
        cat script.metal
        echo "=== END script.metal ==="
        /metal script.metal
    >>>
    output {
        File out = glob(out_prefix + "*" + out_postfix)[0]
    }
}

task concat {
    input {
        Array[File] input_files
        String output_file_name
    }
    File settings_file = write_json(
        object {
            output_file_name: output_file_name
        }
    )
    runtime {
        docker: "us.gcr.io/broad-gdr-dig-storage/metal-python:2018-08-28"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 20 HDD"
    }
    command <<<
        set -e
        python3 --version
        cat << EOF >concat.py
        import json
        import csv
        settingsFile = open("~{settings_file}", "r")
        settings = json.load(settingsFile)
        settingsFile.close()
        print("=== BEGIN settings ===")
        print(json.dumps(settings, sort_keys=True, indent=4))
        print("=== END settings ===")
        in_file_names = ["~{sep='\", \"' input_files}"]
        out_file_name = settings["output_file_name"]
        with open(out_file_name, 'w') as out_file:
            out_writer = csv.writer(out_file, delimiter='\t')
            header_is_written = False
            for in_file_name in in_file_names:
                print("Now going to add " + in_file_name)
                with open(in_file_name, 'r') as in_file:
                    in_reader = csv.reader(in_file, delimiter='\t')
                    header_row = next(in_reader)
                    if not header_is_written:
                        out_writer.writerow(header_row)
                        header_is_written = True
                    for row in in_reader:
                        out_writer.writerow(row)
        print("Done!")
        EOF
        echo "=== BEGIN concat.py ==="
        cat concat.py
        echo "=== END concat.py ==="
        python3 concat.py
    >>>
    output {
        File output_file = output_file_name
    }
}
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
        String output_prefix
        String output_suffix
        String frequency_col = "FREQ"
        String marker_col = "MARKER"
        String size_col = "N"
        String stderr_col = "STDERR"
        String effect_col = "EFFECT"
        String p_value_col = "PVALUE"
        String alt_allele_col = "ALLELE1"
        String ref_allele_col = "ALLELE2"
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
                    frequency_col = frequency_col,
                    input_file = input_file.file,
                    variants_rare_name = ancestry_dataset_file_name + "_rare." + output_suffix,
                    variants_common_name  = ancestry_dataset_file_name + "_common." + output_suffix
            }
        }
        String ancestry_file_name = output_prefix + "_" + ancestry.name
        call metal as metal_samplesize_per_ancestry {
            input:
                input_files = partition.variants_common,
                column_counting = "LENIENT",
                overlap = true,
                out_prefix = ancestry_file_name + "_common_samplesize",
                out_postfix = "." + output_suffix,
                scheme = "SAMPLESIZE",
                average_freq = true,
                min_max_freq = true,
                marker_col = marker_col,
                weight_col = size_col,
                frequency_col = frequency_col,
                stderr_col = stderr_col,
                effect_col = effect_col,
                p_value_col = p_value_col,
                alt_allele_col = alt_allele_col,
                ref_allele_col = ref_allele_col
        }
        call metal as metal_stderr_per_ancestry {
            input:
                input_files = partition.variants_common,
                column_counting = "LENIENT",
                overlap = false,
                marker_col = marker_col,
                weight_col = size_col,
                frequency_col = frequency_col,
                out_prefix = ancestry_file_name + "_common_stderr",
                out_postfix = "." + output_suffix,
                scheme = "STDERR",
                average_freq = true,
                min_max_freq = true,
                stderr_col = stderr_col,
                effect_col = effect_col,
                p_value_col = p_value_col,
                alt_allele_col = alt_allele_col,
                ref_allele_col = ref_allele_col
        }
        call merge_metal_schemes as merge_metal_schemes_per_ancestry {
            input:
                metal_samplesize = metal_samplesize_per_ancestry.out,
                metal_stderr = metal_stderr_per_ancestry.out,
                out_file_name = ancestry_file_name + "_common_merged." + output_suffix
        }
        call pick_largest {
            input:
                rare_variants_files = partition.variants_rare,
                metal_file = merge_metal_schemes_per_ancestry.out,
                marker_col = marker_col,
                frequency_col = frequency_col,
                size_col = size_col,
                stderr_col = stderr_col,
                effect_col = effect_col,
                p_value_col = p_value_col,
                alt_allele_col = alt_allele_col,
                ref_allele_col = ref_allele_col,
                out_file_name = ancestry_file_name + "_picked." + output_suffix
        }
    }
    call metal as metal_samplesize {
        input:
            input_files = pick_largest.output_file,
            column_counting = "LENIENT",
            overlap = false,
            out_prefix = output_prefix + "_samplesize",
            out_postfix = "." + output_suffix,
            scheme = "SAMPLESIZE",
            average_freq = true,
            min_max_freq = true,
            marker_col = "MarkerName",
            weight_col = "Weight",
            frequency_col = "Freq1",
            stderr_col = "StdErr",
            effect_col = "Effect",
            p_value_col = "P-value",
            alt_allele_col = "Allele1",
            ref_allele_col = "Allele2"
    }
    call metal as metal_stderr {
        input:
            input_files = pick_largest.output_file,
            column_counting = "LENIENT",
            overlap = false,
            out_prefix = output_prefix + "_stderr",
            out_postfix = "." + output_suffix,
            scheme = "STDERR",
            average_freq = true,
            min_max_freq = true,
            marker_col = "MarkerName",
            weight_col = "Weight",
            frequency_col = "Freq1",
            stderr_col = "StdErr",
            effect_col = "Effect",
            p_value_col = "P-value",
            alt_allele_col = "Allele1",
            ref_allele_col = "Allele2"
    }
    call merge_metal_schemes {
        input:
            metal_samplesize = metal_samplesize.out,
            metal_stderr = metal_stderr.out,
            out_file_name = output_prefix + "." + output_suffix
    }
}

task partition {
    input {
        Map[String, String] file_column_mapping = { "blub": "blub" }
        Map[String, String] ancestry_column_mapping = { "blub": "blub" }
        Map[String, String] global_column_mapping = { "blub": "blub" }
        Float cutoff_frequency
        String frequency_col
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
            frequency_column: frequency_col,
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

task merge_metal_schemes {
    input {
        File metal_samplesize
        File metal_stderr
        String out_file_name
        String marker_col = "MarkerName"
        String stderr_col = "StdErr"
        String effect_col = "Effect"
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
        cat << EOF > merge.py
        import csv
        samplesize_file_name = "~{metal_samplesize}"
        stderr_file_name = "~{metal_stderr}"
        out_file_name = "~{out_file_name}"
        marker_header = "~{marker_col}"
        stderr_header = "~{stderr_col}"
        effect_header = "~{effect_col}"
        stderr_dict = {}
        with open(stderr_file_name, 'r', newline='') as stderr_file:
            stderr_reader = csv.reader(stderr_file, delimiter='\t')
            stderr_header_row = next(stderr_reader)
            marker_col_index = stderr_header_row.index(marker_header)
            stderr_col_index = stderr_header_row.index(stderr_header)
            effect_col_index = stderr_header_row.index(effect_header)
            for row in stderr_reader:
                marker = row[marker_col_index]
                effect = row[effect_col_index]
                stderr = row[stderr_col_index]
                stderr_dict[marker] = {stderr_header: stderr, effect_header: effect}
        with open(samplesize_file_name, 'r', newline='') as samplesize_file, \
                open(out_file_name, 'w', newline='') as out_file:
            samplesize_reader = csv.reader(samplesize_file, delimiter='\t')
            out_writer = csv.writer(out_file, delimiter='\t')
            samplesize_header_row = next(samplesize_reader)
            marker_col_index = samplesize_header_row.index(marker_header)
            out_header_row = samplesize_header_row + [stderr_header, effect_header]
            out_writer.writerow(out_header_row)
            for row in samplesize_reader:
                marker = row[marker_col_index]
                if marker in stderr_dict:
                    stderr_data = stderr_dict[marker]
                    del stderr_dict[marker]
                    stderr = stderr_data[stderr_header]
                    effect = stderr_data[effect_header]
                    out_row = row + [stderr, effect]
                    out_writer.writerow(out_row)
        EOF
        echo "=== BEGIN merge.py ==="
        cat merge.py
        echo "=== END merge.py ==="
        python3 merge.py
    >>>
    output {
        File out = out_file_name
    }
}

task pick_largest {
    input {
        Array[File] rare_variants_files
        File metal_file
        String out_file_name
        String frequency_col = "FREQ"
        String marker_col = "MARKER"
        String size_col = "N"
        String stderr_col = "STDERR"
        String effect_col = "EFFECT"
        String p_value_col = "PVALUE"
        String alt_allele_col = "ALLELE1"
        String ref_allele_col = "ALLELE2"
    }
    runtime {
        docker: "us.gcr.io/broad-gdr-dig-storage/metal-python:2018-08-28"
        cpu: 1
        memory: "5 GB"
        disks: "local-disk 20 HDD"
    }
    command <<<
        set -e
        python3 --version
        cat << EOF > union.py
        import csv
        rare_variants_file_names = ["~{sep='", "' rare_variants_files}"]
        metal_file = "~{metal_file}"
        marker_col = "MarkerName"
        frequency_col = "Freq1"
        out_col_list = [marker_col, "Weight", frequency_col, "StdErr", "Effect", "P-value", "Allele1",
                        "Allele2"]
        metal_col_list = out_col_list
        rare_col_list = ["~{marker_col}", "~{size_col}", "~{frequency_col}", "~{stderr_col}",
                         "~{effect_col}", "~{p_value_col}", "~{alt_allele_col}", "~{ref_allele_col}"]
        metal_col_dict = dict(zip(metal_col_list, out_col_list))
        rare_col_dict = dict(zip(rare_col_list, out_col_list))
        marker_list = []
        markers = set(marker_list)
        union_data = {}

        def read_file(in_file_name, file_col_dict):
            with open(in_file_name, "r", newline="") as in_file:
                in_reader = csv.reader(in_file, delimiter='\t')
                header_row = next(in_reader)
                col_to_index = {}
                for index, header in enumerate(header_row):
                    out_header = file_col_dict.get(header)
                    if out_header is not None:
                        col_to_index[out_header] = index
                for row in in_reader:
                    row_data = {col: row[index] for col, index in col_to_index.items() }
                    marker = row_data[marker_col]
                    if marker not in markers:
                        marker_list.append(marker)
                        markers.add(marker)
                    union_entry = union_data.get(marker)
                    if union_entry is None:
                        union_data[marker] = row_data
                    else:
                        row_frequency = float(row_data[frequency_col])
                        union_frequency = float(union_entry[frequency_col])
                        if row_frequency > union_frequency:
                            union_data[marker] = row_data


        for in_file_name in rare_variants_file_names:
            read_file(in_file_name, rare_col_dict)
        read_file(metal_file, metal_col_dict)
        out_file_name = "~{out_file_name}"
        with open(out_file_name, 'w') as out_file:
            out_writer = csv.writer(out_file, delimiter='\t')
            out_writer.writerow(out_col_list)
            for marker in marker_list:
                entry = union_data[marker]
                row = list(map(lambda col: entry[col] or "NA", out_col_list))
                out_writer.writerow(row)
        EOF
        echo "=== BEGIN union.py ==="
        cat union.py
        echo "=== END union.py ==="
        python3 union.py
    >>>
    output {
        File output_file = out_file_name
    }
}

task metal {
    input {
        Array[File] input_files
        String column_counting = "STRICT"
        Boolean overlap = false
        String out_prefix
        String out_postfix
        String scheme = "scheme"
        Boolean average_freq = false
        Boolean min_max_freq = false
        String marker_col = "MARKER"
        String weight_col = "N"
        String stderr_col = "STDERR"
        String effect_col = "EFFECT"
        String p_value_col = "PVALUE"
        String frequency_col = "FREQ"
        String alt_allele_col = "ALLELE1"
        String ref_allele_col = "ALLELE2"
        Array[String] custom_variables = []
    }
    File settings_file = write_json(
        object {
            overlap: overlap,
            column_counting: column_counting,
            marker: marker_col,
            weight_column: weight_col,
            out_prefix: out_prefix,
            out_postfix: out_postfix,
            scheme: scheme,
            average_freq: average_freq,
            min_max_freq: min_max_freq,
            std_err: stderr_col,
            effect: effect_col,
            custom_variables: custom_variables,
            p_value: p_value_col,
            freq: frequency_col,
            alt_allele: alt_allele_col,
            ref_allele: ref_allele_col
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
        addValue(lines, "WEIGHTLABEL", settings, "weight_column")
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
        for input_file in ["~{sep='", "' input_files}"]:
            addLine(lines, "PROCESS " + input_file)
        addTwoValues(lines, "OUTFILE", settings, "out_prefix", "out_postfix")
        addLine(lines, "ANALYZE")
        addLine(lines, "QUIT")
        script = "\n".join(lines) + "\n"
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

version 1.0

struct GlobalSettings {
    String? scheme
    Boolean? average_freq
    Boolean? min_max_freq
    String? std_err
    Array[String]? custom_variables
    String out_prefix
    String out_postfix
}

struct FileOptions {
    String? marker
    String? p_value
    String? freq
    String? alt_allele
    String? ref_allele
    String? effect
    Map[String, String]? custom_variable_map
}

struct FileSetting {
    File file
    String dataset_name
    FileOptions? options
}

struct MetalSettings {
    GlobalSettings global_settings
    FileOptions? default_file_options
    Array[FileSetting] file_settings
}

struct InputFile {
    File file
    String dataset_name
    String frequency_header
}

struct Ancestry {
    String name
    Array[InputFile] input_files
}

workflow bottom_line {
    input {
        Array[Ancestry] ancestries
        Float cutoff_frequency
        String output_prefix
        String output_suffix
    }

    scatter(ancestry in ancestries) {
        scatter(input_file in ancestry.input_files) {
            String output_base_name = output_prefix + "_" + ancestry.name + "_" + input_file.dataset_name
            call partition {
                input:
                    cutoff_frequency = cutoff_frequency,
                    frequency_column = input_file.frequency_header,
                    variants = input_file.file,
                    variants_rare_name = output_base_name + "_rare." + output_suffix,
                    variants_common_name  = output_base_name + "_common." + output_suffix
            }
        }
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

task metal {
    input {
        MetalSettings settings
    }
    GlobalSettings global_settings = settings.global_settings
    runtime {
        docker: "us.gcr.io/broad-gdr-dig-storage/metal-python:2018-08-28"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 20 HDD"
    }
    command <<<
        set -e
        cat << EOF >metalcast.py
        import json
        settingsFile = open("~{write_json(settings)}", "r")
        settings = json.load(settingsFile)
        settingsFile.close()
        print("=== BEGIN settings ===")
        print(json.dumps(settings, sort_keys=True, indent=4))
        print("=== END settings ===")
        inputs = map(lambda fileSetting: fileSetting["file"], settings["file_settings"])
        print("=== BEGIN settings ===")
        print(settings)
        print("=== END settings ===")
        globalSettings = settings["global_settings"]
        lines = []
        def addLine(lines, line):
            lines.append(line)
        def addValue(lines, prefix, dict, key, dictDefault = {}, defaultValue = None):
            value = dict.get(key) or dictDefault.get(key) or defaultValue
            if(value is not None):
                lines.append(prefix + " " + value)
        def addFlag(lines, prefix, dict, key, dictDefault = {}, defaultValue = None):
            value = dict.get(key) or dictDefault.get(key) or defaultValue
            if(value is not None):
                if(value):
                    lines.append(prefix + " ON")
                else:
                    lines.append(prefix + " OFF")
        def addTwoValues(lines, prefix, dict, key1, key2, dictDefault = {},
                         defaultValue1 = None, defaultValue2 = None):
            value1 = dict.get(key1) or dictDefault.get(key1) or defaultValue1
            value2 = dict.get(key2) or dictDefault.get(key2) or defaultValue2
            if(value1 is not None and value2 is not None):
                lines.append(prefix + " " + value1 + " " + value2)
        def addArray(lines, prefix, dict, key, dictDefault = {}, defaultValue = None):
            values = dict.get(key) or dictDefault.get(key) or defaultValue
            if(values is not None):
                for value in values:
                    lines.append(prefix + " " + value)
        def addMap(lines, prefix, dict, key, dictDefault = {}, defaultValue = None):
            values = dict.get(key) or dictDefault.get(key) or defaultValue
            if(values is not None):
                for subKey, value in values.items() :
                    lines.append(prefix + " " + subKey + " AS " + value)
        addValue(lines, "SCHEME", globalSettings, "scheme")
        addFlag(lines, "AVERAGEFREQ", globalSettings, "average_freq")
        addFlag(lines, "MINMAXFREQ", globalSettings, "min_max_freq")
        addValue(lines, "STDERR", globalSettings, "std_err")
        addValue(lines, "EFFECT", globalSettings, "effect")
        addArray(lines, "CUSTOMVARIABLE", globalSettings, "custom_variables")
        defaultFileOptions = settings.get("default_file_options") or {}
        for fileSetting in settings["file_settings"]:
            fileOptions = fileSetting.get("options") or {}
            addValue(lines, "MARKER", fileOptions, "marker", defaultFileOptions, "MARKER")
            addValue(lines, "PVALUE", fileOptions, "p_value", defaultFileOptions, "PVALUE")
            addValue(lines, "FREQ", fileOptions, "freq", defaultFileOptions, "FREQ")
            addTwoValues(lines, "ALLELE", fileOptions, "alt_allele", "ref_allele", defaultFileOptions,
                         "ALLELE1", "ALLELE2")
            addValue(lines, "EFFECT", fileOptions, "effect", defaultFileOptions, "EFFECT")
            addMap(lines, "LABEL", fileOptions, "custom_variable_map", defaultFileOptions)
            addLine(lines, "PROCESS " + fileSetting["file"])
        addTwoValues(lines, "OUTFILE", globalSettings, "out_prefix", "out_postfix")
        addLine(lines, "ANALYZE HETEROGENEITY")
        script = reduce(lambda line1, line2: line1 + "\n" + line2, lines) + "\n"
        scriptFile = open("script.metal", "w")
        scriptFile.write(script)
        scriptFile.close()
        EOF
        echo "=== BEGIN metalcast.py ==="
        cat metalcast.py
        echo "=== END metalcast.py ==="
        python metalcast.py
        echo "=== BEGIN script.metal ==="
        cat script.metal
        echo "=== END script.metal ==="
        /metal script.metal
    >>>
    output {
        File out = glob(global_settings.out_prefix + "*" + global_settings.out_postfix)[0]
    }
}
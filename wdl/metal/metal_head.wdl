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
    FileOptions? options
}

struct MetalSettings {
    GlobalSettings global_settings
    FileOptions? default_file_options
    Array[FileSetting] file_settings
}

workflow metal_head {
    input {
        MetalSettings settings
    }
    call metal {
        input:
            settings = settings
    }
    output {
        File out = metal.out
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
        def addValue(lines, prefix, dict, key, dictDefault = {}):
            value = dict.get(key, dictDefault.get(key))
            if(value is not None):
                lines.append(prefix + " " + value)
        def addFlag(lines, prefix, dict, key, dictDefault = {}):
            value = dict.get(key, dictDefault.get(key))
            if(value is not None):
                if(value):
                    lines.append(prefix + " ON")
                else:
                    lines.append(prefix + " OFF")
        def addTwoValues(lines, prefix, dict, key1, key2, dictDefault = {}):
            value1 = dict.get(key1, dictDefault.get(key1))
            value2 = dict.get(key2, dictDefault.get(key2))
            if(value1 is not None and value2 is not None):
                lines.append(prefix + " " + value1 + " " + value2)
        def addArray(lines, prefix, dict, key, dictDefault = {}):
            values = dict.get(key, dictDefault.get(key))
            if(values is not None):
                for value in values:
                    lines.append(prefix + " " + value)
        def addMap(lines, prefix, dict, key, dictDefault = {}):
            values = dict.get(key, dictDefault.get(key))
            if(values is not None):
                for subKey, value in values :
                    lines.append(prefix + " " + subKey + " AS " + value)
        addValue(lines, "SCHEME", globalSettings, "scheme")
        addFlag(lines, "AVERAGEFREQ", globalSettings, "average_freq")
        addFlag(lines, "MINMAXFREQ", globalSettings, "min_max_freq")
        addValue(lines, "STDERR", globalSettings, "std_err")
        addValue(lines, "EFFECT", globalSettings, "effect")
        addArray(lines, "CUSTOMVARIABLE", globalSettings, "custom_variables")
        if("default_file_options" in settings):
            defaultFileOptions = settings["default_file_options"]
            addValue(lines, "MARKER", defaultFileOptions, "marker")
            addValue(lines, "PVALUE", defaultFileOptions, "p_value")
            addValue(lines, "FREQ", defaultFileOptions, "freq")
            addTwoValues(lines, "ALLELE", defaultFileOptions, "alt_allele", "ref_allele")
            addValue(lines, "EFFECT", defaultFileOptions, "effect")
            addLine(lines, "LABEL TotalSampleSize as NS")
        lines += map(lambda item: "PROCESS " + item, inputs)
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
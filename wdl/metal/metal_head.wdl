version 1.0

struct GlobalSettings {
    String? scheme
    Boolean? average_freq
    Boolean? min_max_freq
    String? std_err
    String? effect
    Array[String]? custom_variables
    String? label
    String out_prefix
    String out_postfix
}

struct FileOptions {
    String? marker
    String? p_value
    String? freq
    String? alt_allele
    String? ref_allele
    Map[String, String]? custom_variable_map
}

struct FileSetting {
    File file
    FileOptions? options
}

struct MetalSettings {
    GlobalSettings global_settings
    FileOptions? global_file_options
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
    command {
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
        preamble = """
        SCHEME STDERR
        MARKER ID
        PVALUE P
        FREQ MAF
        ALLELE ALT REF
        AVERAGEFREQ ON
        MINMAXFREQ ON
        STDERR SEBETA
        EFFECT BETA
        CUSTOMVARIABLE TotalSampleSize
        LABEL TotalSampleSize as NS
        """
        processes = reduce(lambda item1, item2: item1 + item2, map(lambda item: "PROCESS " + item + "\n", inputs))
        postamble = """
        OUTFILE """ + globalSettings["out_prefix"] + " " + globalSettings["out_postfix"] + """
        ANALYZE HETEROGENEITY
        """
        script = preamble + "\n" + processes + "\n" + postamble
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
    }
    output {
        File out = glob(global_settings.out_prefix + "*" + global_settings.out_postfix)[0]
    }
}
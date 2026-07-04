"""
CINDEL Snakemake workflow.

Two modes, selected via config["mode"]:
  - "batch":  score every guide in config["input_csv"], write config["output_csv"]
  - "single": score config["single_sequence"], write results/single_score.txt

Usage:
    snakemake -n                       # dry run
    snakemake --cores 1 --use-conda    # real run
"""

configfile: "config/config.yaml"

MODE = config.get("mode", "batch")
VERBOSE_FLAG = "-v" if config.get("verbose", False) else ""


def target_output():
    if MODE == "batch":
        return [config["output_csv"]]
    elif MODE == "single":
        return ["results/single_score.txt"]
    else:
        raise ValueError(f"config['mode'] must be 'batch' or 'single', got {MODE!r}")


rule all:
    input:
        target_output()


rule score_batch:
    input:
        csv=config["input_csv"],
    output:
        config["output_csv"],
    conda:
        "envs/environment.yaml"
    shell:
        "python CINDEL.py -a {input.csv} -o {output} {VERBOSE_FLAG}"


rule score_single:
    output:
        "results/single_score.txt",
    params:
        seq=config["single_sequence"],
    conda:
        "envs/environment.yaml"
    shell:
        "python CINDEL.py -b {params.seq} {VERBOSE_FLAG} > {output}"

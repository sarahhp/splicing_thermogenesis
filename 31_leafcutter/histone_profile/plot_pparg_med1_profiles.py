"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2024-01-25"

CHIP_IN = "/projects/imb-pkbphil/moab/Susanne_Mandrup/Browning_adipocytes_2014/SRA/bamCompare_ratio"


def find_chip(mod, cond):
    sample = list(set(config["modification"][mod]) &
                  set (config["condition"][cond]))[0]
    print(mod, cond, sample)
    fn = "_".join([sample, "hMADS_D19_ratio.bw"])
    filepath = os.path.join(CHIP_IN, fn)
    return filepath

def find_wab_chips(wildcards):
    sample1 = find_chip(wildcards.mod, wildcards.cond1)
    sample2 = find_chip(wildcards.mod, wildcards.cond2)
    return [sample1, sample2]


rule profiles:
    input:
        expand("{mod}/{mod}_{cond1}_{cond2}-{binsize}.profile.pdf",
                mod=["PPARG","MED1"],
                cond1 = "white",
                cond2 = "beige",
                binsize=10)


def choose_window_size(wildcards):
    if wildcards.mod in ["PPARG"]:
        return 2000
    elif  wildcards.mod in ["MED1"]:
        return 4000

rule compute_matrix:
    input:
        histones = find_wab_chips,
        tsses = expand("{dir}/lc_3db_{type}_TSSes.bed",
                        dir=".",
                        type=["beige","white","not_sig"]),
    threads: 4
    params:
        window = choose_window_size,
        binsize = lambda wc: int(wc.binsize),
        samples = lambda wc: expand("{mod}_{cond}",
                        mod = wc.mod,
                        cond = ["White_Adi.","Beige_Adi"])
    output:
        matrix = "{mod}/{mod}_{cond1}_{cond2}-{binsize}.gz",
    shell:
        "computeMatrix reference-point -R {input.tsses} -S {input.histones} "
            "-b {params.window} -a {params.window} "
            "-bs {params.binsize} -p {threads} "
            "--samplesLabel {params.samples} "
            "-o {output.matrix} "

rule plot_profile:
    input:
        matrix = "{mod}/{mod}_{cond1}_{cond2}-{binsize}.gz"
    params:
        title = lambda wc: wc.mod,
        regions = "Beige_TSSs White_TSSs Non-sig_TSSs",
        samples = "White_Adi. Beige_Adi."
    output:
        pdf = "{mod}/{mod}_{cond1}_{cond2}-{binsize}.profile.pdf",
        info = "{mod}/{mod}_{cond1}_{cond2}-{binsize}.profile.tab"
    shell:
        "plotProfile --regionsLabel {params.regions} "
            "--samplesLabel {params.samples} --plotTitle {params.title} "
            "--colors orange c black "
            "-m {input.matrix} -o {output.pdf} --outFileNameData {output.info}"

rule all_counts:
    input:
        expand("{mod}/window{window}/{type}.{window}_white_beige.tab",
                mod="PPARG",
                window = "250:250",
                type = ["beige","white","not_sig"]),
        expand("{mod}/window{window}/{type}.{window}_white_beige.tab",
                mod="MED1",
                window = "500:500",
                type = ["beige","white","not_sig"]),

rule window_TSS:
    input:
        bed = "lc_3db_{type}_TSSes.bed",
        chr_size = "/projects/imb-pkbphil/sp/annotations/ensembl95/hg38.chrom.sizes"
    params:
        up = lambda wc: wc.window.split(":")[0],
        down = lambda wc: wc.window.split(":")[1]
    output:
        "{mod}/window{window}/{type}.{window}.bed"
    shell:
        "bedtools slop -s -l {params.up} -r {params.down} -i {input.bed} "
            "-g {input.chr_size} > {output} "

rule count_enrichment:
    input:
        bed = "{mod}/window{window}/{type}.{window}.bed",
        bws = find_wab_chips
    output:
        counts = "{mod}/window{window}/{type}.{window}_{cond1}_{cond2}.tab",
        npz = "{mod}/window{window}/{type}.{window}_{cond1}_{cond2}.npz"
    shell:
        "multiBigwigSummary BED-file --BED {input.bed} "
            "--bwfiles {input.bws} --labels White_Adi. Beige_Adi. "
            "-o {output.npz} --outRawCounts {output.counts}"













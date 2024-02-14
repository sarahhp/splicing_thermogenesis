"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2024-01-30"


TSS_SET = ["lc_3db_beige_promoters.uniq",
                "lc_3db_white_promoters.uniq",
                "white_and_beige_promoters.merged",
                "lc_3db_not_sig_promoters.merged"]

PPARG_DIR = "PPARG_peaks"
PPARG_SET = ["PPARgamma_hMADS_brite",
        "PPARgamma_hMADS_white"]

rule all:
    input:
        expand("intersects/{tss_set}--{pparg_set}.bed",
                tss_set = TSS_SET,
                pparg_set = PPARG_SET)


rule intersect:
    input:
        tss = "{tss_set}.bed",
        ppar = PPARG_DIR + "/{pparg_set}_peaks.narrowPeak"
    output:
        "intersects/{tss_set}--{pparg_set}.bed"
    shell:
        "bedtools intersect -u -a {input.tss} -b {input.ppar} > {output} "


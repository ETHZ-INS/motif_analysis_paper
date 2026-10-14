path = "/mnt/plger/jwang/data/dat"
path_domcke = "/mnt/plger/jwang/data/domcke"
ttl = expand(path+"/01-total/total-{mot}.rds", mot=MOT)
wgt = expand(path+"/01-weight/weight-{mot}.rds", mot=MOT)
ins = expand(path+"/01-insertion/{mot}.rds", mot=MOT) #+ expand(path_domcke+"/01-insertion/{mot}.rds", mot=MOT2)
bb = expand(path + "/01-bbChromVAR/{mot}.rds", mot = MOT)
data = ttl + wgt + ins 
domcke_peak = "/mnt/plger/jwang/Domcke_scATAC_screen_60TFs/peaks/merged_peaks.narrowPeak"

rule get_origin_counts:
    priority: 99
    input:
        "workflow/code/01-dat_total.R",
        frag = path + "/00-frg/{mot}.rds"
    output: path+"/01-total/total-{mot}.rds"
    log: "logs/get_origin_counts/{mot}.log"
    params:
        peak = lambda wc: (
            glob.glob(f"/mnt/plger/plger/DTFAB/fullFrags/{wc.mot}/peaks/*")[0]),
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} res={output} frag={input[1]} peak={params.peak} genome={params.genome}" {input[0]} {log}'''

rule get_weight_counts:
    priority: 99
    input:
        "workflow/code/01-dat_weight.R",
        frag = path + "/00-frg/{mot}.rds"
    output: path+"/01-weight/weight-{mot}.rds"
    log: "logs/get_weight_counts/{mot}.log"
    params:
        peak = lambda wc: (
            glob.glob(f"/mnt/plger/plger/DTFAB/fullFrags/{wc.mot}/peaks/*")[0]),
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} res={output} frag={input[1]} peak={params.peak} genome={params.genome}" {input[0]} {log}'''

# rule get_weight_insertions:
#     priority: 99
#     input:
#         "workflow/code/01-dat_inserts.R",
#         frag = path + "/00-frg/{mot}.rds",
#         motifRanges = "/mnt/plger/plger/DTFAB/fullFrags/{mot}/runATAC_results/others/pmoi.rds",
#     output:
#         path + "/01-insertion/{mot}.rds"
#     log: "logs/get_weight_insertions/{mot}.log"
#     params:
#         peak = lambda wc: (
#             glob.glob(f"/mnt/plger/plger/DTFAB/fullFrags/{wc.mot}/peaks/*")[0]),
#         genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus"
#     shell: '''
#         {R} CMD BATCH --no-restore --no-save "--args\
#         wcs={wildcards} res={output} frag={input[1]} motfrange={input[2]} peak={params.peak} genome={params.genome}" {input[0]} {log}'''

rule get_weight_insertions:
    priority: 99
    input:
        "workflow/code/01-dat_inserts.R",
        motifRanges = "/mnt/plger/plger/DTFAB/fullFrags/{mot}/runATAC_results/others/pmoi.rds",
    output:
        path + "/01-insertion/{mot}.rds"
    log: "logs/get_weight_insertions/{mot}.log"
    params:
        peak = lambda wc: (
            glob.glob(f"/mnt/plger/plger/DTFAB/fullFrags/{wc.mot}/peaks/*")[0]),
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} res={output} motfrange={input[1]} peak={params.peak} genome={params.genome}" {input[0]} {log}'''


rule get_weight_insertions_domcke:
    priority: 99
    input:
        "workflow/code/01-dat_inserts.R",
        frag = path_domcke + "/00-frg/{mot}.rds",
        motifRanges = path_domcke + "/motifRanges_all.rds",
    output:
        path_domcke + "/01-insertion/{mot}.rds"
    log: "logs/get_weight_insertions/{mot}.log"
    params:
        peak = "/mnt/plger/jwang/Domcke_scATAC_screen_60TFs/peaks/merged_peaks.narrowPeak",
        genome = "Mmusculus"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} res={output} frag={input[1]} motfrange={input[2]} peak={params.peak} genome={params.genome}" {input[0]} {log}'''

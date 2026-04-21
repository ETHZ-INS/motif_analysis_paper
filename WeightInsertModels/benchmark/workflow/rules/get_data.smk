path = "/mnt/plger/jwang/data/dat"
ttl = expand(path+"/01-total/total-{mot}.rds", mot=MOT)
wgt = expand(path+"/01-weight/weight-{mot}.rds", mot=MOT)
data = ttl + wgt

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

# rule dat_ttl:
#     priority: 98
#     input:  "code/01-dat_total.R",
#             path+"/00-frg/{mot}.rds"
#     output: path+"/01-total/total-{mot}.rds"
#     params:
#         genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus",
#         peak = lambda wc: (
#             glob.glob(f"/mnt/plger/plger/DTFAB/fullFrags/{wc.mot}/peaks/*")[0]
#             if data == "Benchmark" else peak2
#         )
#     log:    "logs/01-total-{mot}.Rout"
#     shell: '''
#         {R} CMD BATCH --no-restore --no-save "--args\
#         wcs={wildcards} res={output} frg={input[1]} peak={params.peak} genome={params.genome}" {input[0]} {log}'''
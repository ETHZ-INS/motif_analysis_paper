path = "/mnt/plger/jwang/data/dat"
path_domcke = "/mnt/plger/jwang/data/domcke"
dif_ttl = expand("outs/dat/dif-total-{mot},{dif}.rds", mot = MOT, dif = ["chromVAR.limma", "limma-voom", 
    "betterChromVAR", "betterChromVAR.FL", "chromVAR", "bbChromVAR", "bbChromVAR_masking", "bbChromVAR_unweighted", "bbChromVAR_noFL"]) 
dif_wgt = expand("outs/dat/dif-weight-{mot},{dif}.rds", mot = MOT, dif = ["limma-voom"])
dif_ins = expand("outs/dat/dif-inserts-{mot},{dif}.rds", mot = MOT, dif = ["limma"]) 
#dif_ins =  expand("outs/domcke/dif-inserts-{mot},{dif}.rds", mot = MOT2, dif = ["limma"])
#dif_ttl = expand("outs/domcke/dif-total-{mot},{dif}.rds", mot = MOT2, dif=["chromVAR.limma"])
#diff = dif_ins + dif_ttl
diff = dif_ttl + dif_wgt + dif_ins
#diff = dif_ins

rule differential_test_origin:
    priority: 97
    input:  "workflow/code/02-dif.R",
            "workflow/code/02-dif-{dif}.R",
            path+"/01-total/total-{mot}.rds"
    output: "outs/dat/dif-total-{mot},{dif}.rds"
    params:
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus",
    log:    "logs/differential_test_origin/{mot},{dif}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} fun={input[1]} res={output} dat={input[2]} genome={params.genome}" {input[0]} {log}'''

rule differential_test_weight:
    priority: 97
    input:  "workflow/code/02-dif.R",
            "workflow/code/02-dif-{dif}.R",
            path+"/01-weight/weight-{mot}.rds"
    output: "outs/dat/dif-weight-{mot},{dif}.rds"
    params:
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus",
    log:    "logs/differential_test_weights/{mot},{dif}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} fun={input[1]} res={output} dat={input[2]} genome={params.genome}" {input[0]} {log}'''

rule differential_test_inserts:
    priority: 97
    input:  "workflow/code/02-dif.R",
            "workflow/code/02-dif-{dif}.R",
            path+"/01-insertion/{mot}.rds"
    output: "outs/dat/dif-inserts-{mot},{dif}.rds"
    params:
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus",
    log:    "logs/differential_test_inserts{mot},{dif}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} fun={input[1]} res={output} dat={input[2]} genome={params.genome}" {input[0]} {log}'''

rule differential_test_origin_domcke:
    priority: 97
    input:  "workflow/code/02-dif.R",
            "workflow/code/02-dif-{dif}.R",
            path_domcke+"/01-total/total-{mot}.rds"
    output: "outs/domcke/dif-total-{mot},{dif}.rds"
    params:
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus",
    log:    "logs/differential_test_origin/{mot},{dif}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} fun={input[1]} res={output} dat={input[2]} genome={params.genome}" {input[0]} {log}'''

rule differential_test_inserts_domcke:
    priority: 97
    input:  "workflow/code/02-dif.R",
            "workflow/code/02-dif-{dif}.R",
            path_domcke+"/01-insertion/{mot}.rds"
    output: "outs/domcke/dif-inserts-{mot},{dif}.rds"
    params:
        genome = lambda wildcards: "Hsapiens" if wildcards.mot in human else "Mmusculus",
    log:    "logs/differential_test_inserts{mot},{dif}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} fun={input[1]} res={output} dat={input[2]} genome={params.genome}" {input[0]} {log}'''


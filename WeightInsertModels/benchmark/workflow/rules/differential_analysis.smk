path = "/mnt/plger/jwang/data/dat"
dif_ttl = expand("outs/dat/dif-total-{mot},{dif}.rds", mot = MOT, dif = ["chromVAR.limma", "limma-voom", "betterChromVAR", "betterChromVAR.FL"])
dif_wgt = expand("outs/dat/dif-weight-{mot},{dif}.rds", mot = MOT, dif = ["limma-voom"])
diff = dif_ttl + dif_wgt 

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
    log:    "logs/differential_test_origin/{mot},{dif}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} fun={input[1]} res={output} dat={input[2]} genome={params.genome}" {input[0]} {log}'''

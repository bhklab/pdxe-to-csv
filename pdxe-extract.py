import json
import re
import requests
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()

# Get the JSON from the API
response = requests.get('https://www.orcestra.ca/api/xevaset/available')
objs = json.loads(response.text)

# Find PDXE dataset (Note: You can change "PDXE" to any other xevaset when they become available)
target = next((item for item in objs if item["name"] == "PDXE"), None)
if not target:
    raise ValueError("PDXE not found in available datasets")
download_link = target["downloadLink"]

# Extract the base file name from the downloadLink (This is optional)
match = re.search(r'files/([^?]+)\.rds', download_link)
if match:
    base_filename = match.group(1)
    rds_filename = f"{base_filename}.rds"
    expression_filename = f"{base_filename}_expression.csv"
    model_csv = f"{base_filename}_model.csv"
    expDesign_csv = f"{base_filename}_expDesign.csv"
    experiment_csv = f"{base_filename}_experiment.csv"
else:
    rds_filename = "downloaded.rds"
    expression_filename = "downloaded_expression.csv"
    model_csv = "downloaded_model.csv"
    expDesign_csv = "downloaded_expDesign.csv"
    experiment_csv = "downloaded_experiment.csv"

# Download the .rds file from https://orcestra.ca/api/xevaset/available
rds_response = requests.get(download_link)
with open(rds_filename, "wb") as f:
    f.write(rds_response.content)

# Load required R packages and the RDS file
robjects.r('library(Xeva)')
robjects.r('library(Biobase)')

readRDS = robjects.r['readRDS']
pdxe = readRDS(rds_filename)
robjects.globalenv['pdxe'] = pdxe

# Extract expression matrix, convert, export
exprs_df_r = robjects.r('as.data.frame(exprs(pdxe@molecularProfiles$RNASeq))')
exprs_df = pandas2ri.rpy2py(exprs_df_r)
exprs_df.to_csv(expression_filename, index=True)

# Extract model matrix, convert, export
model_df_r = robjects.r('as.data.frame(pdxe@model)')
model_df = pandas2ri.rpy2py(model_df_r)
model_df.to_csv(model_csv, index=True)

# Extract experiment design, convert, export
expDesign_r = robjects.r('''
    expDesign_list <- lapply(names(pdxe@expDesign), function(nm) {
        x <- pdxe@expDesign[[nm]]
        data.frame(
          experiment = nm,
          name = if (is.null(x$name) || length(x$name)==0) NA else as.character(x$name)[1],
          control = if (is.null(x$control) || length(x$control)==0) NA else as.character(x$control)[1],
          treatment = if (is.null(x$treatment) || length(x$treatment)==0) NA else as.character(x$treatment)[1],
          stringsAsFactors = FALSE)
    })
    expDesign_df <- do.call(rbind, expDesign_list)
    expDesign_df
''')
expDesign_df = pandas2ri.rpy2py(expDesign_r)
expDesign_df.to_csv(expDesign_csv, index=False)

# Extract experiment, convert, export
experiment_r = robjects.r('''
    experiment_list <- lapply(names(pdxe@experiment), function(nm) {
        x <- pdxe@experiment[[nm]]
        if (!is.null(x@data)) {
            df <- as.data.frame(x@data)
            df$experiment <- nm
            df$model_id <- as.character(x@`model.id`)
            df$drug <- as.character(x@drug)
            df
        } else {
            NULL
        }
    })
    experiment_df <- do.call(rbind, experiment_list)
    experiment_df
''')
experiment_df = pandas2ri.rpy2py(experiment_r)
experiment_df.to_csv(experiment_csv, index=False)

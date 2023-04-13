library(IOHanalyzer)


for (fname in list.dirs("D:/Many_Affine/Many_affine2/", full.names = T, recursive = F)) {
  tryCatch({
    dsl1 <- DataSetList(fname)

    fid <- as.integer(stringi::stri_split(basename(fname), fixed='_')[[1]][[1]])

    saveRDS(dsl1, paste0("D:/Many_Affine/rds/", basename(fname), ".rds"))


    dt1 <- generate_data.AUC(dsl1, get_ECDF_targets(dsl1, 'bbob')
                             , multiple_x = T, scale_log = T)

    dt2 <- get_FV_sample(dsl1, c(100,250,500,1000,2000,5000, 10000, 50000, 100000, 200000), output='long')

    dt3 <- get_RT_sample(dsl1, c(1e6, 1e4, 1e2, 1e0, 1e-2, 1e-4, 1e-6, 1e-8), output='long')

    dt2$funcId <- fid
    dt3$funcId <- fid
    dt1$funcId <- fid

    write.csv(dt1, paste0("D:/Many_Affine/csv/auc/", basename(fname), ".csv"))
    write.csv(dt2, paste0("D:/Many_Affine/csv/fv/", basename(fname), ".csv"))
    write.csv(dt3, paste0("D:/Many_Affine/csv/ert/", basename(fname), ".csv"))

  }, error = function(e) { print(fname)})
}


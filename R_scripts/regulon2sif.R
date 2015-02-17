# usage: Rscript regulon2sif.R regulon.RData salida.sif
# regulon.RData es un rdata que contiene un regulon, salida de marina

args <- commandArgs(trailingOnly = TRUE);

rdata = args[1];
outfile = args[2]
print(outfile);

load(rdata);
regulon = get(ls()[2]);


for (name in names(regulon)) {
    source = rep(name, length(regulon[[name]]$tfmod));
    sif    = cbind(source, names(regulon[[name]]$tfmod), regulon[[name]]$tfmod[[1]]);
    write.table(sif, file = outfile, append = TRUE, quote = FALSE, sep = "\t", col.names=FALSE, row.names=FALSE);
}




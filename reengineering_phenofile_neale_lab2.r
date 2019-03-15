library("data.table")

# Neale Lab application.
bd2 <- fread("ukb11214.csv", header=TRUE, sep=',', na.strings=c("NA", ""))
bd2 <- as.data.frame(bd2)

# First get rid of things that are > 4.
colnames(bd2) <- gsub("-", "\\.", colnames(bd2))
visit_number <- as.integer(unlist(lapply(strsplit(names(bd2), "\\."), "[", 2)))[-1]
fields <- unlist(lapply(strsplit(names(bd2), "\\."), "[", 1))[-1]

to_include <- which(visit_number <= 4)
# Get rid of the fields that have up to 31 visits - small collection of Cancer variables,
# need to look at those separately.

# Cancer variables 
cancer_fields <- unique(fields[which(visit_number > 4)])
to_include <- setdiff(to_include, which(fields %in% cancer_fields))

bd2 <- bd2[,c(1,to_include+1)]

# Now, want to go through these column names and determine, if there's a 0, include it, if not, look to 1 etc. 
fields <- unlist(lapply(strsplit(names(bd2), "\\."), "[", 1))[-1]
visit_number <- as.integer(unlist(lapply(strsplit(names(bd2), "\\."), "[", 2))[-1])
unique_fields <- unique(fields)
fields_to_use <- c()

for (i in unique_fields) {
    matches <- which(fields == i)
    match_min <- which(visit_number[matches] == min(visit_number[matches]))
    match <- matches[match_min]
    fields_to_use <- c(fields_to_use, match)
}

bd2 <- bd2[,c(1,fields_to_use+1)]

c1 <- gsub("^", "x", colnames(bd2))
c2 <- gsub("\\.", "_", c1)
c3 <- c("userId", c2[2:length(c2)])
colnames(bd2) <- c3

# Hack to ensure that the age and sex are in as columns.
if (any(colnames(bd2) == "x21022_0_0") == FALSE) {
    bd2 <- cbind(rep(1, nrow(bd2)), bd2)
    colnames(bd2)[1] <- "x21022_0_0"
}

if (any(colnames(bd2) == "x31_0_0") == FALSE) {
    bd2 <- cbind(rep(1, nrow(bd2)), bd2)
    colnames(bd2)[1] <- "x31_0_0"
}

fwrite(bd2, file="neale_lab_parsed.tsv", row.names=FALSE, quote=FALSE, sep='\t')
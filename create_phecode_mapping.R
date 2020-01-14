library(dplyr)
library(stringr)
library(purrr)
library(data.table)

# Inputs are from the interim outputs of https://github.com/umich-cphds/createUKBphenome.
icd9key =
  fread("./data/UKB_PHENOME_ICD9_PHECODE_MAP_20200109.txt", colClasses = "character", data.table = F) %>%
  rename(icd_code = ICD9) %>%
  mutate(icd_version = "icd9")

icd10key =
  fread("./data/UKB_PHENOME_ICD10_PHECODE_MAP_20200109.txt", colClasses = "character", data.table = F) %>%
  rename(icd_code = ICD10) %>%
  mutate(icd_version = "icd10")

icd_all = bind_rows(icd10key, icd9key) %>%
  group_by(phecode, sex, description, group) %>%
  summarize(
    icd_codes = str_c(str_remove(icd_code, "\\."), collapse = ',')
  ) %>%
  ungroup %>%
  mutate(
    sex = recode(sex, Both = "both_sexes", Male = "males", Female = "females")
  )

if (length(unique(icd_all$phecode)) != nrow(icd_all)) {
  stop("Inconsistent sex def.")
}

pheinfo = fread("./data/PHECODE_v1.2b1_INFO_20200109.txt", colClasses = "character", data.table = F)
pheinfo2 = subset(pheinfo, phecode %in% icd_all$phecode)

# source: https://github.com/umich-cphds/createUKBphenome/blob/master/scripts/function.expandPhecodes.r
expandPhecodes <- function(x,addIntegers=T){
    if(is.na(x) | x == "") return(NA)
    if(grepl("\\-",x)){
        # split range
        # character prefix
        i1 <- strsplit(x,"-")[[1]]

        # numeric length of digits before "."
        nprefix <- max(nchar(gsub("\\..+","",i1)))
        # numbers of digits
        ndigits <- max(c(nchar(gsub("^[0-9]+[\\.]{0,1}","",i1)),0))
        # add "." to length of formatted number if present
        addDot <- max(as.numeric(grepl("\\.",i1)))
        # create sequence of ICD codes
        seq1 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^ndigits))
        # format sequence to match intput
        seq1 <- formatC(seq1, format='f', digits=ndigits,width=nprefix+ndigits+addDot,flag=0)
        # add integers if within range
        if(addIntegers) seq1 <- unique(sort(c(seq1,gsub("\\..+","",seq1[which(round(as.numeric(seq1)) == as.numeric(seq1))]))))

        if(ndigits == 2){
            seq2 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^(ndigits-1)))
            seq2 <- formatC(seq2, format='f', digits=ndigits-1,width=nprefix+ndigits+addDot-1,flag=0)
            seq1 <- unique(sort(c(seq1,seq2)))
        }
        return(seq1)
    } else {
        return(x)
    }
}

icd_all$exclude_phecodes =
  map_chr(1:nrow(pheinfo2), function(p) {
  phecode_remove <- ""
  phecode <- pheinfo2$phecode[p]

  # collect phecodes to include from controls
  exclude_phecodes <- phecode
  if(pheinfo2$phecode_exclude_range[p] != ""){
    phecode_remove <- unlist(strsplit(gsub(" ","",pheinfo2$phecode_exclude_range[p]),",")[[1]])
    exclude_phecodes <- c(exclude_phecodes,unlist(sapply(phecode_remove,function(x) expandPhecodes(x,T),USE.NAMES=F)))
  }
  exclude_phecodes <- unique(exclude_phecodes[which(exclude_phecodes %in% pheinfo2$phecode)])
  return(str_c(sort(exclude_phecodes), collapse=","))
})

write.table(icd_all, "./data/UKB_Phecode_v1.2b1_ICD_Mapping.txt", quote = F, row.names = F, sep = "\t")


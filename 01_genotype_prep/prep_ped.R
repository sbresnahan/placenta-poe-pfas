# =============================================================================
# Stage 01 | prep_ped.R
# Infant-centred trio pedigree from the three PLINK .fam files
# (FID IID PID MID SEX PHENO; "0" = parent not genotyped).
# In : TRIOS_DIR/<PLINK_PREFIXES>.fam      Out: TRIO_PED
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
library(data.table)

##### Prepare pedigree (.ped) file #####

# Read in the fam files (PLINK_PREFIXES order: infant, mother, father)
fam_of <- function(i) file.path(P("TRIOS_DIR"), paste0(words("PLINK_PREFIXES")[i], ".fam"))
infant_fam <- fread(fam_of(1), header = FALSE, colClasses = "character")
mother_fam <- fread(fam_of(2), header = FALSE, colClasses = "character")
father_fam <- fread(fam_of(3), header = FALSE, colClasses = "character")

# Assign column names
colnames(infant_fam) <- colnames(father_fam) <- colnames(mother_fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENOTYPE")

# Merge infant with mother and father IIDs based on Family ID (FID)
# Left joins ensure that all infants are retained even if parent is missing
infant_mothers <- merge(infant_fam[, .(FID, IID, SEX)], 
                        mother_fam[, .(FID, motherID = IID)], 
                        by = "FID", all.x = TRUE)

trio_fam <- merge(infant_mothers, 
                  father_fam[, .(FID, fatherID = IID)], 
                  by = "FID", all.x = TRUE)

# Replace missing IDs with "0"
trio_fam[is.na(fatherID), fatherID := "0"]
trio_fam[is.na(motherID), motherID := "0"]

# Create new .fam data: FID, IID, PID, MID, SEX, PHENOTYPE
trio_fam_final <- trio_fam[, .(FID,
                               IID,
                               PID = fatherID,
                               MID = motherID,
                               SEX,
                               PHENOTYPE = -9)]

# Write to file
fwrite(trio_fam_final, P("TRIO_PED"),
       sep = " ", col.names = FALSE, quote = FALSE)
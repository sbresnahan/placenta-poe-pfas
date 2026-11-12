# =============================================================================
# Stage 00 | imprintDB.R
# Curated imprinted-gene symbols and aliases -> Ensembl gene IDs (biomaRt), plus
# Ensembl ID -> symbol for every gene in the placenta long-read transcriptome.
# In : TX_CLASSIF
# Out: IMPRINT_DB (geneimprint_db_ensembl.txt), PLACENTA_DB (placenta_db_ensembl.txt)
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))

library(biomaRt)

# imprint_db
# Vector of all gene names and aliases
genes_aliases <- c(
  "DIRAS3","ARHI","NOEY2",
  "RNU5D-1","U5DL","U5DS","RNU5D",
  "TP73","P73",
  "LRRTM1",
  "GPR1",
  "ZDBF2","Slx9",
  "NAP1L5","DRLM",
  "RHOBTB3",
  "ERAP2","LRAP","L-RAP",
  "VTRNA2-1","CBL3","CBL-3","hvg-5","nc886","MIR886","VTRNA2","MIRN886","svtRNA2-1a","hsa-mir-886",
  "NR3C1","GR","GCR","GRL","GCCR","GCRST",
  "ADTRP","AIG1L","C6orf105","dJ413H6.1",
  "PXDC1","C6orf145",
  "FAM50B","X5L","D6S2654E",
  "LIN28B","CSDD2",
  "AIM1","ST4","CRYBG1",
  "PLAGL1","ZAC","LOT1","ZAC1","MGC126275","MGC126276","DKFZp781P1017",
  "HYMAI","NCRNA00020",
  "SLC22A2","OCT2","MGC32628",
  "SLC22A3","EMT","EMTH","OCT3",
  "GRB10","RSS","IRBP","MEG1","GRB-IR","Grb-10","KIAA0207",
  "DDC","AADC",
  "GLI3","PHS","ACLS","GCPS","PAPA","PAPB","PAP-A","PAPA1","PPDIV",
  "HECW1","NEDL1",
  "HOXA4","HOX1","HOX1D",
  "PEG10","EDR","HB-1","Mar2","MEF3L","Mart2","RGAG3",
  "MAGI2","AIP1","AIP-1","ARIP1","SSCAM","MAGI-2","ACVRIP1",
  "SGCE","ESG","DYT11",
  "PPP1R9A","NRB1","NRBI","FLJ20068","KIAA1222","Neurabin-I",
  "TFPI2","PP5","REF1","TFPI-2","FLJ21164",
  "CCDC71L","C7orf74",
  "CPA4","CPA3",
  "COPG2IT1","CIT1","COPG2AS","FLJ41646","NCRNA00170","DKFZP761N09121",
  "MEST","PEG1","MGC8703","MGC111102","DKFZp686L18234",
  "MESTIT1","MEST-IT","PEG1-AS","MEST-AS1","MEST-IT1","NCRNA00040",
  "KLF14","BTEB5",
  "SVOPL",
  "PTPRN2","IAR","ICAAR","PTPRP","IA-2beta","R-PTP-N2",
  "DLGAP2","DAP2","SAPAP2",
  "KCNK9","KT3.2","TASK3","K2p9.1","TASK-3","MGC138268","MGC138270",
  "PEG13",
  "MYC","MRTL","MYCC","c-Myc","bHLHe39",
  "ZFAT","AITD3","ZFAT1","ZNF406",
  "ZFAT-AS1","ZFATAS","ZFAT-AS","SAS-ZFAT","NCRNA00070",
  "GLIS3","ZNF515",
  "VIM",
  "WT1-AS","WIT1","WIT-1","WT1AS","WT1-AS1",
  "WT1","GUD","AWT1","WAGR","WT33","NPHS4","WIT-2","EWS-WT1",
  "KCNQ1DN","BWRT","HSA404617",
  "OSBPL5","ORP5","OBPH1","FLJ42929",
  "PHLDA2","IPL","BRW1C","BWR1C","HLDA2","TSSC3",
  "INS","ILPR","IRDN",
  "H19","ASM","BWS","ASM1","MGC4485","PRO2605","D11S813E",
  "SLC22A18AS","BWR1B","BWSCR1B","ORCTL2S","SLC22A1LS","p27-BWR1B",
  "CDKN1C","BWS","WBS","p57","BWCR","KIP2",
  "IGF2AS","PEG8","MGC168198",
  "KCNQ1","LQT","RWS","WRS","LQT1","SQT2","ATFB1","ATFB3","JLNS1","KCNA8","KCNA9","Kv1.9","Kv7.1","KVLQT1","FLJ26167",
  "KCNQ1OT1","LIT1","Kncq1","KvDMR1","KCNQ10T1","KCNQ1-AS2","KvLQT1-AS","NCRNA00012",
  "IGF2","INSIGF","pp9974","C11orf43","FLJ22066","FLJ44734",
  "SLC22A18","HET","ITM","BWR1A","IMPT1","TSSC5","ORCTL2","BWSCR1A","SLC22A1L","p45-BWR1A","DKFZp667A184",
  "ANO1","DOG1","TAOS2","ORAOV2","TMEM16A",
  "ZC3H12C","MCPIP3",
  "NTM","HNT","NTRI","CEPU-1","IGLON2",
  "ST8SIA1","GD3S","SIAT8","SIAT8A","SIAT8-A","ST8SiaI",
  "RBP5","CRBP3","CRBPIII","CRBP-III",
  "ATP5F1EP2","ATP5EP2",
  "RB1","RB","pRb","OSRC","pp110","p105-Rb","PPP1R130","p110-RB1",
  "SMOC1","OAS",
  "DIO3","D3","5DIII","TXDI3","DIOIII",
  "MEG3","GTL2","FP504","prebp1","PRO0518","PRO2160","FLJ31163","FLJ42589",
  "DLK1","DLK","FA1","ZOG","pG2","DLK-1","PREF1","Delta1","Pref-1",
  "MEG8","Bsr","Irm","Rian","SNHG23","SNHG24","lnc-MGC","LINC00024","NCRNA00024",
  "SNORD114-1","14q(II-1)",
  "SNORD113-1","14q(I-1)",
  "DIO3OS","DIO3-OS","DIO3-AS1","C14orf134","NCRNA00041",
  "RTL1","MART1","PEG11","LOC388015",
  "MAGEL2","nM15","NDNL1",
  "UBE3A","AS","ANCR","E6-AP","HPVE6A","EPVE6AP","FLJ26981",
  "MKRN3","D15S9","RNF63","ZFP127","ZNF127","MGC88288",
  "SNORD115-48","HBII-52-48",
  "SNORD115","HBII-52",
  "SNORD116","PET1","PWCR1","HBII-85",
  "ATP10A","ATPVA","ATPVC","ATP10C","KIAA0566",
  "SNORD109A","HBII-438A",
  "SNORD107","HBII-436","HBII-436 C/D box snoRNA",
  "SNORD108","HBII-437","HBII-437 C/D box snoRNA",
  "SNRPN","SMN","PWCR","SM-D","RT-LI","HCERN3","SNRNP-N","FLJ33569","FLJ36996","FLJ39265","MGC29886","SNURF-SNRPN","DKFZp762N022","DKFZp686C0927","DKFZp761I1912","DKFZp686M12165",
  "SNORD109B","HBII-438B","HBII-438B C/D box snoRNA",
  "PWAR6","HBT8","PAR-6",
  "PWCR1","PET1","non-coding RNA in the Prader-Willi critical region",
  "NPAP1","C15orf2",
  "NDN","HsT16328",
  "SNORD64","HBII-13","HBII-13 snoRNA",
  "SNURF",
  "RASGRF1","GNRP","GRF1","CDC25","GRF55","CDC25L","H-GRF55","PP13187",
  "IRAIN","IGF1R-AS",
  "ZNF597","HIT-4",
  "PRR25","gs64",
  "NAA60","HAT4","NatF","IBGC9","NAT15","hNaa60",
  "CMTM1","CKLFH","CKLFH1","CKLFSF1",
  "ZFP90","FIK","NK10","ZNF756","zfp-90",
  "ZNF396","ZSCAN14",
  "TCEB3C","HsT829","TCEB3L2","Elongin A3",
  "PARD6G","PAR-6G","PAR6gamma",
  "PARD6G-AS1",
  "ADNP2","ZNF508",
  "DNMT1","AIM","DNMT","MCMT","CXXC9","HSN1E","ADCADN",
  "AXL","ARK","UFO","JTK11","Tyro7",
  "ZIM2","ZNF656",
  "PEG3","PW1","ZNF904","ZSCAN24",
  "MIR371A","C19MC","MIR371","MIRN371","hsa-mir-371","hsa-mir-371a",
  "NLRP2","NBS1","PAN1","NALP2","PYPAF2","CLR19.9",
  "PEG3-AS1","APEG3","PEG3AS","PEG3-AS","NCRNA00155",
  "MIMT1","MIM1","LINC00067","NCRNA00067",
  "BLCAP","BC10",
  "NNAT","Peg5",
  "MCTS2",
  "GDAP1L1","dJ881L22.1","dJ995J12.1.1",
  "SGK2","H-SGK2","dJ138B7.2",
  "GNAS","AHO","GSA","GSP","POH","GPSA","NESP","GNAS1","PHP1A","PHP1B","C20orf45","MGC33735","dJ309F20.1.1","dJ806M20.3.3",
  "L3MBTL1","L3MBTL","ZC2HC3","H-L(3)MBT","dJ138B7.3",
  "MIR298","MIRN298","hsa-mir-298",
  "GNAS-AS1","SANG","GNASAS","NESPAS","GNAS-AS","GNAS1AS","NESP-AS","NCRNA00075",
  "MIR296","MIRN296","miRNA296",
  "MIR125B2","MIRN125B2","mir-125b-2",
  "DSCAM","CHD2","CHD2-42","CHD2-52",
  "DGCR6L",
  "DGCR6"
)

# Deduplicate and clean
imprint_db <- unique(gsub("\\*", "", genes_aliases))

imprint_db.cap <- toupper(imprint_db)

imprint_db <- unique(c(imprint_db,imprint_db.cap))

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query Ensembl using gene symbols
ensembl_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "hgnc_symbol",
  values     = imprint_db,
  mart       = mart
)

write.table(ensembl_map,P("IMPRINT_DB"),col.names=T,row.names=F,sep="\t",quote=F)



# Query Ensembl using gene symbols for all placenta genes
placenta.genes <- read.table(P("TX_CLASSIF"),header=T,sep="\t")$associated_gene
placenta.genes <- sub("\\..*$", "", placenta.genes)

ensembl_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = placenta.genes,
  mart       = mart
)

write.table(ensembl_map,P("PLACENTA_DB"),col.names=T,row.names=F,sep="\t",quote=F)

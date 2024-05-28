set.seed(3112)

#Load TF
TF <- c("Mcr-002592",
        "Mcr-009840",
        "Mcr-001329",
        "Mcr-009574",
        "Mcr-017531",
        "Mcr-007661",
        "Mcr-000904",
        "Mcr-004424",
        "Mcr-006861",
        "Mcr-012616",
        "Mcr-001431")

#Load matrix
Matrix <- read.csv("RNASeq_Matrix.csv", row.names=1)

# Identify complete cases
complete_cases <- complete.cases(Matrix)

# Subset the matrix to keep only complete cases
Matrix <- Matrix[complete_cases, ]
Matrix <- as.matrix(Matrix)

# Not all the TFs are present in the dataset -> new vector for regulators specific for this dataset. 

new.regulators <- TF[match(rownames(Matrix),TF)]
new.regulators <- new.regulators[!is.na(new.regulators)] #remove NAs

#Run GENIE3
weightmat_mes <- GENIE3(Matrix, regulators = new.regulators, nCores = 15)
linklist_mes <- getLinkList(weightmat_mes)
write.csv(linklist_mes, file = "Genie3.csv")

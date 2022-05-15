make_methyl <- function(bed_path,sample_path) {
#read in bed file with the probe locations 
    
    ed <- read.table(bed_path,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    
    chr <- ed$V1
    start <- ed$V2
    stop <-  ed$V3
    junk1 <- ed$V4
    junk2 <- ed$V5
    
    #sort and get rid of duplicates
    uchr <- sort(unique(chr))
    
    filepath <- sample_path
    #read in cgmap file that has the C and T counts at each position
    
    cgmap <- read.table(filepath,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    V1 <- tools::file_path_sans_ext(sample_path)
    chrp <- cgmap$V1
    nuc <- cgmap$V2
    pos<- cgmap$V3
    CG<- cgmap$V4
    CGP<- cgmap$V5
    meth<- cgmap$V6
    countsC<- cgmap$V7
    countsTot<- cgmap$V8
    
    #compare uchr with chrp and find indices where they are the same
    #looks for each uchr{1} (i.e. chr1) and finds all indices of it in chrp
    idx <- list()
    for (i in seq(1:length(uchr))) {
        idx[[i]] <- which(chrp == uchr[i])
    }
    
    idx2 <- list()
    methmat <- matrix(0,nrow=length(chr),ncol=1)  
    countsmat <- matrix(0,nrow=length(chr),ncol=1)  
    maxcountsmat <-matrix(0,nrow=length(chr),ncol=1)  
    for (i in seq(1:length(chr))) {
        idx2 <- (uchr == chr[i])
        idxp <- unlist(idx[idx2])
        distbp <- pos[idxp] - start[i] + 60
        idxprobe <- which(abs(distbp) < 100)
        
        sumC = sum(countsC[idxp[idxprobe]])
        sumTot = sum(countsTot[idxp[idxprobe]])
        maxTot = max(countsTot[idxp[idxprobe]])
        
        
        if (sumTot > 0) {
            methmat[i,1] = sumC/sumTot;
            countsmat[i,1] = sumTot;
            maxcountsmat[i,1] = maxTot;
        }
        else {
            methmat[i,1] = NaN;
            countsmat[i,1] = NaN;
            maxcountsmat[i,1] = NaN;
        }
        
    
    }
methmat <- rbind(V1, methmat)
return(as.data.frame(t(methmat)))
}


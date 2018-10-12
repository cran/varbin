#########################
### Numerical binning ###
#########################

varbin <- function(df, x, y, p=0.05, custom_vec=NA){

  if (!is.data.frame(df)) {
    return("Data not a data.frame")
  }
  else if (is.numeric(y) | is.numeric(x)) {
    return("Column name not string")
  }
  else if (grepl("[.]", y) | grepl("[.]", x)) {
    return("Column name with a dot [.]")
  }
  else i = which(names(df) == y)
  j = which(names(df) == x)
  if (!is.numeric(df[, i])) {
    return("Target (y) not found or it is not numeric")
  }
  else if (max(df[, i], na.rm = T) != 1) {
    return("Maximum not 1")
  }
  else if (min(df[, i], na.rm = T) != 0) {
    return("Minimum not 0")
  }
  else if (p <= 0 | p > 0.5) {
    return("p must be greater than 0 and lower than 0.5 (50%)")
  }
  else if (!is.numeric(df[, j])) {
    return("Characteristic (x) not found or it is not a number")
  }else{
    df <- data.frame(df, stringsAsFactors=T)
    df_out <- data.frame(NULL, stringsAsFactors=F)

    if (is.na(custom_vec[1])){
      f <- as.formula(paste0(y," ~ ", x))
      tree <- rpart(f, data=df, na.action=na.omit,
                    control=rpart.control(cp=0.0001, maxdepth=20, minbucket=round(nrow(df)*p)))
      bin_vec <- sort(as.numeric(tree$splits[,4]))
      bin_vec <- c(bin_vec, bin_vec[length(bin_vec)])
      if (length(bin_vec)<=1){
        return("No split possible")
      }
    }else{
      bin_vec <- sort(custom_vec)
      bin_vec <- c(bin_vec, bin_vec[length(bin_vec)])
    }

    for (i in 1:length(bin_vec)){
      if (i==1){
        y_sub <- df[((df[x]<bin_vec[i]) & (!is.na(df[x]))), y]
        cp <- paste0("<", bin_vec[i])
      }else if (i==length(bin_vec)){
        y_sub <- df[((df[x]>=bin_vec[i]) & (!is.na(df[x]))), y]
        cp <- paste0(">=", bin_vec[i])
      }else{
        y_sub <- df[((df[x]>=bin_vec[i-1]) & (df[x]<bin_vec[i]) & (!is.na(df[x]))), y]
        cp <- paste0("<", bin_vec[i])
      }
      n   <- length(y_sub)
      n1 <- sum(y_sub)
      n0 <- n-n1
      df_out <- rbind(df_out, data.frame(Cutpoint=cp, CntRec=n, CntNoEvent=n0, CntEvent=n1))
    }

    n_na <- length(df[is.na(df[x]), y])
    n1_na <- sum(df[is.na(df[x]), y])
    n0_na <- n_na-n1_na
    cp_na <- "Missing"
    df_out <- rbind(df_out, data.frame(Cutpoint=cp_na, CntRec=n_na, CntNoEvent=n0_na, CntEvent=n1_na))

    n_total <- nrow(df)
    n1_total <- sum(df[y])
    n0_total <- n_total-n1_total

    df_out$CntCumRec <- cumsum(df_out$CntRec)
    df_out$CntCumNoEvent <- cumsum(df_out$CntNoEvent)
    df_out$CntCumEvent  <- cumsum(df_out$CntEvent)
    df_out$PctRec <- df_out$CntRec/n_total
    df_out$NoEventRate <- df_out$CntNoEvent/df_out$CntRec
    df_out$EventRate <- df_out$CntEvent/df_out$CntRec
    df_out$Odds <- df_out$EventRate/df_out$NoEventRate
    df_out$LnOdds <- log(df_out$EventRate/df_out$NoEventRate)
    df_out$WoE <- log((df_out$CntNoEvent/n0_total)/(df_out$CntEvent/n1_total))*100
    df_out$IV <- ((df_out$CntNoEvent/n0_total)-(df_out$CntEvent/n1_total))*df_out$WoE

    cp_total <- "Total"
    df_out <- rbind(df_out, data.frame(Cutpoint=cp_total, CntRec=n_total, CntNoEvent=n0_total,
                                       CntEvent=n1_total, CntCumRec=NA, CntCumNoEvent=NA, CntCumEvent=NA,
                                       PctRec=sum(df_out$PctRec), NoEventRate=n0_total/n_total,
                                       EventRate=n1_total/n_total,
                                       Odds=(n1_total/n_total)/(n0_total/n_total),
                                       LnOdds=log((n1_total/n_total)/(n0_total/n_total)), WoE=0.0,
                                       IV=sum(df_out$IV, na.rm=T)))

    return(df_out)
  }
}
######################
### Factor binning ###
######################
varbin.factor <- function(df, x, y, custom_vec=NA){
  if (!is.data.frame(df)) {
    return("Data not a data.frame")
  }
  else if (is.numeric(y) | is.numeric(x)) {
    return("Column name not string")
  }
  else if (grepl("[.]", y) | grepl("[.]", x)) {
    return("Column name with a dot [.]")
  }
  else i = which(names(df) == y)
  j = which(names(df) == x)
  if (!is.numeric(df[, i])) {
    return("Target (y) not found or it is not numeric")
  }
  else if (max(df[, i], na.rm = T) != 1) {
    return("Maximum not 1")
  }
  else if (any(grepl(",", df[, j]))) {
    return("Values contain comma")
  }
  else if (tolower(y) == "default") {
    return("Field name 'default' not allowed")
  }
  else if (min(df[, i], na.rm = T) != 0) {
    return("Minimum not 0")
  }
  else if (!is.factor(df[, j])) {
    return("Characteristic (x) not found or it is not a factor")
  }
  else if (length(unique(df[, j])) <= 1) {
    return("Characteristic (x) requires at leats 2 uniques categories")
  }else{
    df <- data.frame(df, stringsAsFactors=T)
    df_out <- data.frame(NULL, stringsAsFactors=F)

    if (is.na(custom_vec[1])){
      for (bin in unique(df[!is.na(df[x]), x])){
        y_sub <- df[((df[x]==bin) & (!is.na(df[x]))), y]
        cp <- paste0("=", bin)
        n <- length(y_sub)
        n1 <- sum(y_sub)
        n0 <- n-n1
        df_out <- rbind(df_out, data.frame(Cutpoint=cp, CntRec=n, CntNoEvent=n0, CntEvent=n1))
      }
    }else{
      for (bin in custom_vec){
        bin_clean <- strsplit(bin, split = ",")[[1]]
        y_sub <- df[((df[, x] %in% bin_clean) & (!is.na(df[x]))), y]
        cp <- paste0("=", bin)
        n <- length(y_sub)
        n1 <- sum(y_sub)
        n0 <- n-n1
        df_out <- rbind(df_out, data.frame(Cutpoint=cp, CntRec=n, CntNoEvent=n0, CntEvent=n1))
      }
    }

    n_na <- length(df[is.na(df[x]), y])
    n1_na <- sum(df[is.na(df[x]), y])
    n0_na <- n_na-n1_na
    cp_na <- "Missing"
    df_out <- rbind(df_out, data.frame(Cutpoint=cp_na, CntRec=n_na, CntNoEvent=n0_na, CntEvent=n1_na))

    n_total <- nrow(df)
    n1_total <- sum(df[y])
    n0_total <- n_total-n1_total

    df_out$CntCumRec <- cumsum(df_out$CntRec)
    df_out$CntCumNoEvent <- cumsum(df_out$CntNoEvent)
    df_out$CntCumEvent  <- cumsum(df_out$CntEvent)
    df_out$PctRec <- df_out$CntRec/n_total
    df_out$NoEventRate <- df_out$CntNoEvent/df_out$CntRec
    df_out$EventRate <- df_out$CntEvent/df_out$CntRec
    df_out$Odds <- df_out$EventRate/df_out$NoEventRate
    df_out$LnOdds <- log(df_out$EventRate/df_out$NoEventRate)
    df_out$WoE <- log((df_out$CntNoEvent/n0_total)/(df_out$CntEvent/n1_total))*100
    df_out$IV <- ((df_out$CntNoEvent/n0_total)-(df_out$CntEvent/n1_total))*df_out$WoE

    cp_total <- "Total"
    df_out <- rbind(df_out, data.frame(Cutpoint=cp_total, CntRec=n_total, CntNoEvent=n0_total,
                                       CntEvent=n1_total, CntCumRec=NA, CntCumNoEvent=NA, CntCumEvent=NA,
                                       PctRec=sum(df_out$PctRec), NoEventRate=n0_total/n_total,
                                       EventRate=n1_total/n_total,
                                       Odds=(n1_total/n_total)/(n0_total/n_total),
                                       LnOdds=log((n1_total/n_total)/(n0_total/n_total)), WoE=0.0,
                                       IV=sum(df_out$IV, na.rm=T)))
    return(df_out)
  }
}
###############################################
### Monotonically in- or decreasing binning ###
###############################################
varbin.monotonic <- function(df, x, y, p=0.05, increase=F, decrease=F, auto=T){
  if (!is.data.frame(df)) {
    return("Data not a data.frame")
  }
  else if (is.numeric(y) | is.numeric(x)) {
    return("Column name not string")
  }
  else if (grepl("[.]", y) | grepl("[.]", x)) {
    return("Column name with a dot [.]")
  }
  else i = which(names(df) == y)
  j = which(names(df) == x)
  if (!is.numeric(df[, i])) {
    return("Target (y) not found or it is not numeric")
  }
  else if (max(df[, i], na.rm = T) != 1) {
    return("Maximum not 1")
  }
  else if (min(df[, i], na.rm = T) != 0) {
    return("Minimum not 0")
  }
  else if (p <= 0 | p > 0.5) {
    return("p must be greater than 0 and lower than 0.5 (50%)")
  }
  else if (!is.numeric(df[, j])) {
    return("Characteristic (x) not found or it is not a number")
  }else{
    df <- data.frame(df, stringsAsFactors=T)
    # Slice data
    d1  <- df[c(y, x)]
    # Ensure no Nas for regression
    d2  <- d1[!is.na(d1[x]), ]
    # Estimate rank correlation to determine sign of slope
    c   <- cor(d2[, 2], d2[, 1], method = "spearman", use = "complete.obs")
    # Perform isotonic regression
    if (increase==F && decrease==F && auto==T){
      reg <- isoreg(d2[, 2], c / abs(c) * d2[, 1])
    }else if (increase==T && decrease==F && auto==F){
      reg <- isoreg(d2[, 2], d2[, 1])
    }else if (increase==F && decrease==T && auto==F){
      reg <- isoreg(d2[, 2], -d2[, 1])
    }
    # Get x-axis step-values
    k   <- unique(knots(as.stepfun(reg)))+min(diff(sort(unique(d1[, x]), na.last=NA)))/10
    # Perform binning on raw x-values
    sm1 <- varbin(d1, x, y, p, custom_vec=k)
    # Filter out non-informative bins and bins with too few obs.
    CntNoEvent <- "im a dummy"
    CntEvent <- "im a dummy"
    CntRec <- "im a dummy"
    Cutpoint <- "im a dummy"
    c1  <- subset(sm1, subset=CntNoEvent * CntEvent > 0 & CntRec>(nrow(df)*p), select=Cutpoint)
    # c1  <- sm1[(((sm1[,"CntNoEvent"]*sm1[,"CntEvent"])>0) &
    #             (sm1[,"CntRec"]>(nrow(df)*p))), "Cutpoint"]
    # Extract the resulting cutpoints and convert to numerics
    c2  <- apply(c1, 1, function(x) as.numeric(gsub("[^0-9\\.]", "", x)))
    # Remove Nas
    c3  <- c2[!is.na(c2)]
    if (length(c3)<=1){
      return("No Bins")
    }else{
      return(varbin(d1, x, y, p, custom_vec=c3[-length(c3)]))
    }
  }
}
#####################################
### Global/local extremum binning ###
#####################################
varbin.kink <- function(df, x, y, p=0.05){
  if (!is.data.frame(df)) {
    return("Data not a data.frame")
  }
  else if (is.numeric(y) | is.numeric(x)) {
    return("Column name not string")
  }
  else if (grepl("[.]", y) | grepl("[.]", x)) {
    return("Column name with a dot [.]")
  }
  else i = which(names(df) == y)
  j = which(names(df) == x)
  if (!is.numeric(df[, i])) {
    return("Target (y) not found or it is not numeric")
  }
  else if (max(df[, i], na.rm = T) != 1) {
    return("Maximum not 1")
  }
  else if (min(df[, i], na.rm = T) != 0) {
    return("Minimum not 0")
  }
  else if (p <= 0 | p > 0.5) {
    return("p must be greater than 0 and lower than 0.5 (50%)")
  }
  else if (!is.numeric(df[, j])) {
    return("Characteristic (x) not found or it is not a number")
  }else{
    df <- data.frame(df, stringsAsFactors=T)
    WoE_monotonic_pos <- varbin.monotonic(df, x, y, p, increase=T, decrease=F, auto=F)
    WoE_monotonic_neg <- varbin.monotonic(df, x, y, p, increase=F, decrease=T, auto=F)

    if (class(WoE_monotonic_pos)=="character"){
      return("Monotonically increase not possible")
    }else if (class(WoE_monotonic_neg)=="character"){
      return("Monotonically decrease not possible")
    }else {
      pos_cuts <- WoE_monotonic_pos$Cutpoint
      neg_cuts <- WoE_monotonic_neg$Cutpoint
      pos_cuts <- sapply(pos_cuts, function(x) as.numeric(gsub("[^0-9\\.]", "", x)))
      neg_cuts <- sapply(neg_cuts, function(x) as.numeric(gsub("[^0-9\\.]", "", x)))
      pos_cuts <- pos_cuts[!is.na(pos_cuts)]
      neg_cuts <- neg_cuts[!is.na(neg_cuts)]
      pos_cuts <- pos_cuts[-length(pos_cuts)]
      neg_cuts <- neg_cuts[-length(neg_cuts)]

      if (length(neg_cuts)>=length(pos_cuts)){
        med_cut   <- median(neg_cuts)
        start_cut <- pos_cuts[1]
        end_cut   <- pos_cuts[length(pos_cuts)]

        if(start_cut<med_cut & end_cut>med_cut){
          t1 <- sum(pos_cuts<med_cut)
          t2 <- sum(pos_cuts>med_cut)

          if(t1>=t2){
            idx1 <- which(pos_cuts<med_cut)
            idx2 <- which(neg_cuts>=med_cut)
            cutpoints <- c(pos_cuts[idx1], neg_cuts[idx2])
          }else{
            idx1 <- which(pos_cuts>med_cut)
            idx2 <- which(neg_cuts<=med_cut)
            cutpoints <- c(neg_cuts[idx2], pos_cuts[idx1])
          }

        }else if (end_cut<med_cut){
          idx <- which(end_cut<neg_cuts)
          cutpoints <- c(pos_cuts, neg_cuts[idx])
        }else if (start_cut>med_cut){
          idx <- which(start_cut>neg_cuts)
          cutpoints <- c(neg_cuts[idx], pos_cuts)
        }else{
          NA
        }
      }else{
        med_cut   <- median(pos_cuts)
        start_cut <- neg_cuts[1]
        end_cut   <- neg_cuts[length(neg_cuts)]

        if(start_cut<med_cut & end_cut>med_cut){
          t1 <- sum(neg_cuts<med_cut)
          t2 <- sum(neg_cuts>med_cut)

          if(t1>=t2){
            idx1 <- which(neg_cuts<med_cut)
            idx2 <- which(pos_cuts>=med_cut)
            cutpoints <- c(neg_cuts[idx1], pos_cuts[idx2])
          }else{
            idx1 <- which(neg_cuts>med_cut)
            idx2 <- which(pos_cuts<=med_cut)
            cutpoints <- c(pos_cuts[idx2], neg_cuts[idx1])
          }

        }else if (end_cut<med_cut){
          idx <- which(end_cut<pos_cuts)
          cutpoints <- c(neg_cuts, pos_cuts[idx])
        }else if (start_cut>med_cut){
          idx <- which(start_cut>pos_cuts)
          cutpoints <- c(pos_cuts[idx], neg_cuts)
        }else{
          NA
        }
      }
      return(varbin(df, x, y, p, custom_vec=cutpoints))
    }
  }
}
######################
### Bin conversion ###
######################
varbin.convert <- function(df, ivTable, x){
  if (!is.data.frame(df)) {
    return("Data not a data.frame")
  }else if (!is.data.frame(ivTable)){
    return("Table not a data.frame")
  }else{
    #Initialize values
    df <- data.frame(df, stringsAsFactors=T)
    df$tmp <- NA
    ivTable <- ivTable[, c("Cutpoint", "WoE")]
    ivTable <- ivTable[ivTable[,1]!="Total", ]
    ivTable[, 2] <- as.numeric(ivTable[, 2])

    if (class(df[, x])=='factor'){
      ivTable[, 1] <- sub("=", "", ivTable[, 1])
      # Handle merged bin values
      s <- strsplit(as.character(ivTable[, 1]), split=",")
      ivTable <- data.frame(Cutpoint=unlist(s), WoE=rep(ivTable[, 2], sapply(s, length)),
                            stringsAsFactors=F)
      ivTable <- ivTable[!is.na(ivTable[, 2]), ]

      # Assign WoE values to tmp variable
      for (bin in ivTable[, 1]){
        if (bin=="Missing"){
          df$tmp[is.na(df[x])] <- ivTable[ivTable[, 1]==bin, 2]
        }else{
          df$tmp[df[x]==bin] <- ivTable[ivTable[, 1]==bin, 2]
        }
      }
    }else{
      # Assign WoE values to tmp variable
      i <- 0
      for (bin in ivTable[, 1]){
        i <- i+1
        if (regexpr('>=', bin)[1]>-1){
          bin_clean <- as.numeric(sub(">=", "", bin))
          df$tmp[((df[x]>=bin_clean) & (!is.na(df[x])))] <- ivTable[ivTable[, 1]==bin, 2]

        }else if((regexpr('<', bin)[1]>-1) & (i==1)){
          bin_clean <- as.numeric(sub("<", "", bin))
          df$tmp[((df[x]<bin_clean) & (!is.na(df[x])))] <- ivTable[ivTable[, 1]==bin, 2]

        }else if((regexpr('<', bin)[1]>-1) & (i>=2)){
          bin_clean_upper <- as.numeric(sub("<", "", bin))
          bin_clean_lower <- as.numeric(sub("<", "", ivTable[, 1][i-1]))
          df$tmp[((df[x]>=bin_clean_lower) & (df[x]<bin_clean_upper) & (!is.na(df[x])))] <- ivTable[ivTable[, 1]==bin, 2]

        }else if(bin=='Missing'){
          df$tmp[is.na(df[x])] <- ivTable[ivTable[, 1]==bin, 2]
        }
      }
    }
    new_name <- paste0('WoE_', x)
    if (new_name %in% colnames(df)){
      df <- df[, !(names(df) %in% new_name)]
    }
    colnames(df)[colnames(df) == 'tmp'] <- new_name
    return(df)
  }
}
################
### Bin plot ###
################
varbin.plot <- function(ivTable){
  if (!is.data.frame(ivTable)){
    return("Table is not a data.frame")
  }
  x <- ivTable[, 'Cutpoint']
  y <- ivTable[, 'WoE']

  if (is.na(ivTable[ivTable$Cutpoint=='Missing', "WoE"])){
    x <- x[-c(length(x)-1, length(x))]
    y <- y[-c(length(y)-1, length(y))]
    col <- c(rep("black",length(x)))
    pch <- c(rep(19,length(x)))
  }else{
    x <- x[-length(x)]
    y <- y[-length(y)]
    col <- c(rep("black",length(x)-1))
    pch <- c(rep(19,length(x)-1),1)
  }
  par(mar=c(2,2,2,2))
  plot(1:length(x), y, xlab="Cutpoint", ylab="WoE", xaxt="n", cex=1, pch=pch, cex.lab=1, col=col)
  abline(h=0, lty=2, col='black')
  title(main="Weight of Evidence", line=0.25)
  axis(1, at=1:length(x), labels=x)
}

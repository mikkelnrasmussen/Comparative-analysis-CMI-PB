# Load ensembl table for mapping between gene symbols/synonyms and Ensembl IDs
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl")
gene.map <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                  mart = ensembl)
gene.synonym <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                                     'external_synonym'),  mart = ensembl)

###############################################################################

### Helper function for selecting the correct genes in the expression matrix
### using gene symbols and gene synonyms

###############################################################################

selectGeneExp <- function(exp.mat, gene.set, verbose=FALSE){
   
   # Identify the genes map onto the gene expression matrix using the gene 
   # symbol
   module.exp <- merge(exp.mat, gene.map[,c(1,2)],
                       by.x=0, by.y="ensembl_gene_id") %>%
      distinct(Row.names, .keep_all = TRUE) %>% 
      tibble::column_to_rownames('Row.names') %>% 
      filter(external_gene_name %in% gene.set) %>% 
      distinct(external_gene_name, .keep_all = TRUE)
   
   if(verbose){
      # Print messages regarding the genes in the gene set
      print(paste0("Number of genes in total: ", length(unique(gene.set))))
      print("Genes in module: ")
      print(sort(gene.set))  
   }
   
   if(!all(gene.set %in% module.exp$external_gene_name)){
      
      # Identify the missing gene(s)
      missing.genes <- gene.set[!(gene.set %in% module.exp$external_gene_name)]
      
      if(verbose){
         # Print message regarding missing genes
         print("Missing gene(s): ")
         print(missing.genes)
      }
      
      # Identify the missing genes in the synonym mapping dataframe
      missing.id <- gene.synonym[gene.synonym$external_synonym 
                                 %in% missing.genes, ]
      
      # Get the gene expression for the missing gene(s) 
      missing.exp <- merge(exp.mat, missing.id[, c(1, 3)],
                           by.x=0, by.y="ensembl_gene_id") %>%
         dplyr::rename(external_gene_name = external_synonym) %>%
         distinct(Row.names, .keep_all = TRUE) %>% 
         tibble::column_to_rownames('Row.names') %>% 
         distinct(external_gene_name, .keep_all = TRUE)
      
      # Add to the existing expression matrix
      module.exp <- rbind(module.exp, missing.exp)
      
      if(verbose){
         # Display message regarding the identification of genes
         print(paste0("Number of identified genes: ", 
                      length((module.exp$external_gene_name))))
         print("Identified genes:")
         print(sort(module.exp$external_gene_name))  
      }
      
   } else{
      if(verbose){
         # Display message regarding the identification of genes
         print(paste0("Number of identified genes: ", 
                      length((module.exp$external_gene_name))))
         print("Identified genes:")
         print(sort(module.exp$external_gene_name))  
      }
   }
   
   # Drop the gene symbols
   module.exp <- module.exp %>% 
      dplyr::select(-one_of("external_gene_name"))
   
   return(module.exp)
}

###############################################################################

### Initialize functions for calculating the gene signature score and BioAge 
### (gene module score)

###############################################################################

# arithMean_comb_validate function from Model 1: Multicohort analysis reveals 
# baseline transcriptional predictors of influenza vaccination responses (PMID: 
# 28842433) - modified to work with CMI-PB data structures and workflow 
# - (M. N. Rasmussen, Feb. 28th 2022)
# Source: 
# https://github.com/RGLab/ImmuneSignatures/blob/main/R/Single_Gene_Analysis.R
arithMean_comb_validate <- function(genes_p = NULL,
                                    genes_n = NULL,
                                    scale   = 1){
   df.gene.up <- as.data.frame(genes_p)
   df.gene.down <- as.data.frame(genes_n)
   
   # Compute score with mean
   score <- rep(0, length(c(colnames(genes_p), colnames(genes_n))))
   if( !is.null(genes_p) ){
      score <- score + colMeans(df.gene.up, na.rm=T)
   }
   if( !is.null(genes_n) ){
      score <- score - colMeans(df.gene.down, na.rm=T)
   }
   if( !is.null(genes_p) || !is.null(genes_n) ){
      if( scale == 1 ){
         score <- scale(score);
      }
   }
   
   # Transform shape of data
   score <- as.data.frame(t(score))
   df.score <- score %>%
      pivot_longer(cols = everything(), names_to = "subject_id", 
                   values_to = "predictor")
   df.score$subject_id <- as.numeric(df.score$subject_id)
   
   return(df.score)
}

# get_score function from Model 3: Broad immune activation underlies shared set 
# point signatures for vaccine responsiveness in healthy individuals and disease 
# activity in patients with lupus (32094927)
# Source:
# https://github.com/niaid/baseline/tree/master/R/functions 

# Calculate score from genes x samples matrix as average z-score
get_score <- function(x) {
    # Calculate score
    x = t(scale(t(x)))
    score <- as.data.frame(t(colMeans(x, na.rm=T)))
    
    # Transform shape of data
    df.score <- score %>%
        pivot_longer(cols = everything(), names_to = "subject_id", 
                     values_to = "predictor")
    df.score$subject_id <- as.numeric(df.score$subject_id)
    
    return(df.score)
}


# Function for calculating the gene signature score for studies, where no 
# calculation was presented
calculateGeneSigScore <- function(gene.ex.list, gencode, genes.up=NULL, 
                                  genes.down=NULL, ensembl.flag=FALSE){
    
    # Convert to dataframe
    df.gene.ex <- as.data.frame(gene.ex.list)
    
    # Get gene names and ids
    if(ensembl.flag){
        gene.ids <- rownames(gene.ex.list)
        gene.names <- gene.ids
    } else {
        gene.names <- unique(gencode$external_gene_name[gencode$ensembl_gene_id 
                                                        %in% rownames(df.gene.ex)])
        gene.ids <- unique(gencode$ensembl_gene_id[gencode$external_gene_name 
                                                   %in% gene.names])
    }
    
    # If only the gene signature score is provided, then set the up-regulated 
    # genes equal to the signature genes
    if(is.null(genes.up) & is.null(genes.down)){
        genes.up <- gene.names
    }
    
    # For each sample, sum normalized TPM counts for upregulated genes (Sum Up) 
    # and downregulated genes (Sum Down)
    if( !(is.null(genes.up)) ){
        genes.up <- unique(genes.up[genes.up %in% gene.names])
        genes.up.ids <- gene.ids[gene.names %in% genes.up]
        
        sum.up <- df.gene.ex %>%
            dplyr::filter(rownames(df.gene.ex) %in% genes.up.ids) %>% 
            dplyr::summarise_all(sum)
    } 
    
    if( !(is.null(genes.down)) ){
        genes.down <- unique(genes.down[genes.down %in% gene.names])
        genes.down.ids <- gene.ids[gene.names %in% genes.down]
        
        sum.down <- df.gene.ex %>%
            dplyr::filter(rownames(df.gene.ex) %in% genes.down.ids) %>% 
            dplyr::summarise_all(sum)
        
    }
    
    # For each sample, subtract Down to Up: X= (Sum Up) - (Sum Down)
    if(!is.null(genes.up) & !is.null(genes.down)){
        X <- sum.up - sum.down
    } else {
        X <- sum.up
    }
    
    # Calculate the average Av, and standard deviation Stdev for X across all 
    # your samples
    Av <- rowMeans(X)
    Stdev <- sd(X)
    
    # Calculate the square root of n, the total number of samples in your 
    # dataset: sqrt (n)
    squareRoot <- sqrt(length(X))
    
    # Calculate a z-score for each sample: (X - Av) / (Stdev/sqrt(n))
    score <- (X - Av)/(Stdev/squareRoot)
    df.score <- score %>%
        pivot_longer(cols = everything(), names_to = "subject_id", 
                     values_to = "predictor")
    df.score$subject_id <- as.numeric(df.score$subject_id)
    
    return(df.score)
}

# Function for the calculation of the BioAge score comes from the study by Fourati 
# et al. (Pre-vaccination inflammation and B-cell signalling predict age-related 
# hyporesponse to hepatitis B vaccination - PMID: 26742691) - modified to work 
# with CMI-PB data structures and workflow - (M. N. Rasmussen, Mar. 2nd 2022)
# Source: 
# https://storage.googleapis.com/em131_20120924/Fig3ab/20141114_HBV.Fig3ab.code.pdf

calculateBioAge <- function(ex.mat, gene.map, gene.synonym, 
                            predictive.modules=NULL, verbose=FALSE){
    
    ### Read and calculate the BioAge score on the CMI-PB dataset
    # Download Supplementary Table 3a (BioAge)
    supTable3aPath <- file.path("https://storage.googleapis.com/em131_20120924",
                                "Fig2bc/20141118_HBV.SupplementaryTable3a.csv")
    supTable3aBin <- getBinaryURL(url = supTable3aPath, followlocation = TRUE) 
    
    # read Supplementary Table 3a
    supTable3a <- read_csv(file = supTable3aBin, comment = "#")
    
    # Filter genes in modules that do not have a gene symbol
    module2gene <- supTable3a %>%
        dplyr::select(`Gene Symbol`, Module) %>%
        filter(`Gene Symbol` != "---") %>%
        distinct() %>%
        unstack()
    
    # Map gene names or synonyms to Ensembl IDs and check they are in the
    # dataset
    module2ids <- sapply(module2gene, FUN = function(genes) {
      
      # Generate the gene module expression matrix
      module.exp <- selectGeneExp(exp.mat=ex.mat, gene.set=genes, 
                                  verbose=verbose)
       
      # Extract the ensembl IDs and return
      ensembl.ids <- rownames(module.exp)
      return(value = ensembl.ids)
    })
    
    # Calculate BioAge signature on CMI-PB dataset
    module2score <- t(sapply(module2ids, FUN = function(ids) {
        mat <- ex.mat[ids, ]
        mat <- t(scale(t(mat)))
        return(value = colMeans(mat, na.rm = TRUE))
    }))
    
    # module M1 to M11 are negatively correlated to age, they will have a 
    # negative # sign
    # module M12 to M20 are positively correlated to age, they will have a 
    # positive # sign
    module2sign <- c(rep(-1, times = 11), rep(1, times = 9))
    
    # BioAge = ave(M12:20) - ave(M1:11)
    moduleDF <- as.data.frame(cbind(module2score, module2sign))
    
    if(is.null(predictive.modules)){
       
       # Calculate the average score of positives and negatives modules, separately 
       score <- moduleDF %>%
          pivot_longer(c(-module2sign), names_to = "subject_id", 
                       values_to = "score") %>% 
          group_by(module2sign, subject_id) %>%
          dplyr::summarize(meanScore = mean(score), .groups = 'drop') %>% 
          # substract the average score of modules positively correlated to age 
          # by the average score of module negatively correlated to age
          pivot_wider(names_from = module2sign, values_from = meanScore) %>% 
          group_by(subject_id) %>%
          dplyr::summarise(predictor = `1` - `-1`, .groups = 'drop')
       
       score$subject_id <- as.numeric(score$subject_id)
       
    } else if(!is.null(predictive.modules)){
       
       # Internal function to create BioAge-like score
       calcBioAge <- function(score, sign) {
          bioage <- as.data.frame(cbind(score, sign)) %>%
             pivot_longer(c(-sign), names_to = "subject_id", 
                          values_to = "score") %>% 
             group_by(sign, subject_id) %>%
             dplyr::summarize(meanScore = mean(score), .groups = 'drop') %>%
             ungroup() %>%
             pivot_wider(names_from = sign, values_from = meanScore) %>% 
             group_by(subject_id) %>%
             mutate(`-1` = ifelse(test = "-1" %in% names(.), yes = `-1`, no = 0),
                    `1` = ifelse(test = "1" %in% names(.), yes = `1`, no = 0)) %>%
             dplyr::summarise(predictor = `1` - `-1`, .groups = 'drop')
          return(value = bioage)
       }
       index <- predictive.modules
       score <- calcBioAge(score = module2score[index, , drop = FALSE],
                           sign  = module2sign[index])
       
       score$subject_id <- as.numeric(score$subject_id)
       
    }else {
       stop("Not a valid gene module selected. Choose a number between 1 and 20.")
    }
    
    return(score)
}

###############################################################################

### Functions for evaluating the significance of the results - using AUC scores 
### and Spearman's correlation

###############################################################################

# Function for calculating the area under the ROC curve (AUC) score - 
# N permutations procedure is used for estimating one-tailed P-values
rocEval <- function(data, class, score, seed.value, N.perm,
                    ab.ag, direction){
    
    if(direction %in% c("<", ">")){
        r <- pROC::roc(response=data[[class]], predictor=data[[score]], 
                 direction=direction, 
                 quiet = T)
        
        set.seed(seed.value)
        r.random <- vector("numeric", N.perm)
        for(i in 1:N.perm) {
            r.random[i] <- pROC::roc(response = data[[class]], 
                                     predictor = sample(data[[score]]), 
                                     direction = direction, quiet=T)$auc
        }
        r.p <- (sum(r.random >= r$auc) + 1) / (N.perm + 1)
        r.df <- data.frame(isotype_antigen = ab.ag, 
                           Specificity = r$specificities, 
                           Sensitivity = r$sensitivities) %>% 
            arrange(Sensitivity)
        
    }else if(direction == 'missing'){
        r <- pROC::roc(response=data[[class]], predictor=data[[score]], 
                 direction="auto", 
                 quiet = T)
        
        set.seed(seed.value)
        r.random <- vector("numeric", N.perm)
        for(i in 1:N.perm) {
            r.random[i] <- pROC::roc(response = data[[class]], 
                                        predictor = sample(data[[score]]), 
                                        direction = r$direction, quiet=T)$auc
        }
        r.p <- (sum(r.random >= r$auc) + 1) / (N.perm + 1)
        r.df <- data.frame(isotype_antigen = ab.ag, 
                           Specificity = r$specificities, 
                           Sensitivity = r$sensitivities) %>% 
            arrange(Sensitivity)
        
        # title <- paste("Distribution of AUCs - ",
        #                ab.ag)
        # subtitle <- paste(N.perm, "random resamples of the order of the gene signature scores")
        # 
        # dist.plot <- data.frame(x=r.random) %>%
        #     ggplot() + 
        #     geom_histogram(mapping=aes(x=x, y=..count../sum(..count..)), 
        #                    bins=30) +
        #     scale_y_continuous(labels = scales::percent) +
        #     labs(title=title,
        #          subtitle=subtitle) + 
        #     xlab("AUC score") + 
        #     ylab("Frequency")
        # 
        # plot(dist.plot)
        
    }
    
    return(list('r'=r, 'r.p.value'=r.p, 'r.df'=r.df))
}

# Function for calculating the Spearman's correlation between the gene signature
# score and the day 14 titers
# N permutations procedure is used for estimating one-tailed P-values
corrEval <- function(data, response, score, seed.value, N.perm, 
                     direction, ab.ag, alternative='greater'){
    
        corr <- cor.test(x=data[[score]],
                         y=data[[response]], method = 'spearman',
                         alternative=alternative, exact=FALSE)
        
        corr.random <- vector("numeric", N.perm)
        set.seed(seed.value)
        for(i in 1:N.perm) {
            corr.random[i] <- cor.test(x=sample(data[[score]]),
                                     y=data[[response]], method = 'spearman',
                                     alternative=alternative, exact=FALSE)$estimate  
        }
        
    if(direction %in% c("<", ">")){
        if(direction == "<"){
            corr.p <- (sum(corr.random >= corr$estimate) + 1) / (N.perm + 1)   
        }else if(direction == ">"){
            corr.p <- (sum(corr.random <= corr$estimate) + 1)/ (N.perm + 1) 
        }
    } else if (direction == "missing"){
        corr.p <- (sum(abs(corr.random) >= abs(corr$estimate)) + 1) / (N.perm + 1)
        
        # title <- paste("Distribution of Spearman's correlation coeffient - ",
        #                ab.ag)
        # subtitle <- paste(N.perm, "random resamples of the order of the gene signature scores")
        # 
        # dist.plot <- data.frame(x=abs(corr.random)) %>%
        #     ggplot() + 
        #     geom_histogram(mapping=aes(x=x, y=..count../sum(..count..)), 
        #                    bins=30) +
        #     scale_y_continuous(labels = scales::percent) +
        #     labs(title=title,
        #          subtitle=subtitle) + 
        #     xlab("Spearman's correlation coefficent") + 
        #     ylab("Frequency")
        
        # print(ab.ag)
        # print(abs(corr$estimate))
        # print(corr.p)
        # print(quantile(abs(corr.random), 0.95))
        # print(dist.plot)
        # 
    }
    
    return(list('corr'=corr, 'corr.p.value'=corr.p))
    
}

###############################################################################

### Functions for producing different plots of the evaluation

###############################################################################

# Function for plotting the Spearman's correlations in a heatmap with indication 
# of significance
heatmapPlotter <- function(data, score, save.plot, path, stars,
                           filename, main.title, sub.title=""){
    
    corr.label <- round(data[[score]], 2)
    
    heatmap <- data %>% 
        ggplot(aes_string(x='isotype_antigen', y=1, fill=score, 
                          label=corr.label)) +
        geom_tile() +
        labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation",
             title=main.title, 
             subtitle=sub.title) + 
        scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", 
                             limits=c(-1,1)) +
        theme_classic() +
        geom_text(size=3) +
        geom_text(aes_string(label=stars), color="black", size=5, y=0.55) + 
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) + 
        theme(aspect.ratio=1/10, 
              axis.text.x = element_text(angle = 45, vjust = 0.5, 
                                         hjust=0.4), 
              text=element_text(family="Times New Roman", colour='black'))
    
    if(save.plot){
        ggsave(plot=heatmap, paste0(path, filename, ".png"), w=8,h=2)
    }
}

# Function for plotting of BioAge and day 14 titer values - showing linear 
# regression fit and Spearman's correlation coefficient as well as p-value
corrPlotter <- function(data, x.var, y.var, save.plot, path, filename, 
                        main.title, df.metrics, label.corr, label.corr.pv, 
                        sub.title="", xlabel='', ylabel='Titer value'){
    
    corrPlot <- data %>% 
       ggplot(aes_string(x=x.var, y=y.var)) +
       geom_point(pch = 21, size=1.5, alpha = 0.8, aes(fill=isotype_antigen)) +
       geom_smooth(formula=y~x, method=lm, se=FALSE, fullrange=TRUE, 
                   size=0.7, color='black') +
       theme(strip.background = element_rect(colour="black",
                                             fill="white")) +
       labs(title=main.title,
            subtitle=sub.title) +
       scale_fill_discrete(name='Antibody-antigen pairs') +
       xlab(xlabel) +
       ylab(ylabel) +
       facet_wrap(. ~isotype_antigen) +
       geom_text(aes_string(label=label.corr.pv, x=44, 
                            y = 45 - 0.2*45),
                 data=df.metrics) +
       geom_text(aes_string(label=label.corr, x=44, 
                            y = 45),
                 data=df.metrics) +
       scale_x_continuous(limits = c(0, 50)) +
       scale_y_continuous(limits = c(0, 50))
    
    if(save.plot){
        ggsave(plot=corrPlot, paste0(path, filename, ".png"), w=15, h=11)
    }
}

# Function for plotting boxplots comparing the titer class and the BioAge
boxPlotter <- function(data, response.class, score, save.plot, path, filename, 
                       main.title, sub.title="", xlabel="Response class", 
                       ylabel="Gene signatire score"){
    
    boxplots <- ggplot(data, aes_string(x=response.class, y=score, 
                                        group=response.class)) +
        geom_boxplot(aes_string(fill=response.class), alpha=1, 
                     outlier.colour=NA) +
        geom_dotplot(binaxis = "y", stackdir = "center",
                     aes_string(fill=response.class)) +
        facet_wrap(~isotype_antigen) +
        xlab(xlabel) +
        ylab(ylabel) +
        theme_bw() +
        theme(legend.position="none") +
        theme(panel.border = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size=12),
              panel.spacing = unit(0,"mm"),
              panel.grid.major.x = element_blank(),
              axis.ticks = element_blank()) +
        labs(title=main.title,
             subtitle=sub.title)
    
    if(save.plot){
        ggsave(plot=boxplots, paste0(path, filename, ".png"), w=12, h=11)
    }
}

# Function for plotting ROC curve comparing the titer class and the BioAge
rocPlotter <- function(data, spec, sens, save.plot, main.title, sub.title,
                       df.metrics, label.auc, label.auc.pv, path, filename){
    # Generate AUC plot
    rocplot <- ggplot(data, aes_string(x=spec, y=sens)) +
        geom_line(size=1) +
        geom_abline(intercept = 1, slope = 1, lty=2) +
        scale_x_reverse()  +
        coord_fixed() + 
        theme_bw() + 
        theme(panel.grid = element_blank()) +
        facet_wrap(~isotype_antigen) + 
        geom_text(data=df.metrics, x=-0.20, y=0.25, size=3, 
                  aes_string(label=label.auc), inherit.aes = F) + 
        geom_text(data=df.metrics, x=-0.15, y=0.15, size=3, 
                  aes_string(label=label.auc.pv), inherit.aes = F) +
        labs(title = main.title,
             subtitle = sub.title)
    
    if(save.plot){
        ggsave(plot=rocplot, paste0(path, filename, ".png"), w=12, h=11)
    }
}

###############################################################################

### Main function for evaluting the gene models

###############################################################################


# Function for calculating the gene score for the model 
geneSigCalc <- function(gene.ex.data, post.ab.data=NULL, pre.ab.data=NULL, 
                        gene.sig=NULL, ensembl.sig=NULL, genes.up=NULL, 
                        genes.down=NULL, score='mean', ensembl.flag=FALSE, 
                        filter.ab.ag=NULL, verbose=FALSE){
    
    if(score %in% c('mean', 'geometric', 'missing')){
        
        if(!is.null(gene.sig) & is.null(ensembl.sig)){
           
           # Generate the gene signature expression matrix
           gene.sig.exp <- selectGeneExp(exp.mat=gene.ex.data, gene.set=gene.sig,
                                         verbose=verbose)
           
        } else {
           print(paste("Number of genes in total:", length(unique(ensembl.sig))))
           print(ensembl.sig)
           ensembl.ids <- ensembl.sig
           ensembl.flag <- TRUE
           
           # Subset gene expression data to only include genes from signature
           gene.mask <- toupper(rownames(gene.ex.data)) %in% toupper(ensembl.ids)
           gene.sig.exp <- gene.ex.data[gene.mask, ]
           print(paste("Number of genes in the gene expression matrix:", 
                       nrow(gene.sig.exp)))
           print(rownames(gene.sig.exp))
        }
      
        # Calculate the gene signature score
        if(score=='mean'){
            pred <- get_score(gene.sig.exp)
        }else if(score=='geometric'){
            pred <- arithMean_comb_validate(genes_p=gene.sig.exp, 
                                                  genes_n=NULL,
                                                  scale=1)
        } else if(score=='missing'){
            pred <- calculateGeneSigScore(gene.ex.list=gene.sig.exp, 
                                          gencode=gene.map, 
                                          genes.up=genes.up,
                                          genes.down=genes.down,
                                          ensembl.flag=ensembl.flag)
        } else {
           stop("The selected score is not valid. Please select one of the 
                following: 'mean', ''geometric' or 'missing'")
        }
    }
   
   if(!(is.null(pre.ab.data) & is.null(post.ab.data))){
      
      # Subset antibody titers data
      ab.mask <- post.ab.data$subject_id %in% colnames(gene.ex.data)
      ab.mask.pre <- pre.ab.data$subject_id %in% colnames(gene.ex.data)
      
      # Subset antibody titer data
      ab.data.post <- post.ab.data[ab.mask, ] %>% 
       dplyr::select(subject_id, isotype_antigen, MFI_normalised) %>% 
       dplyr::rename(post.titer = MFI_normalised)
      
      ab.data.pre <- pre.ab.data[ab.mask.pre, ] %>% 
       dplyr::select(subject_id, isotype_antigen, MFI_normalised) %>% 
       dplyr::rename(pre.titer = MFI_normalised)
      
      # Create dataframe with prediction and response classes
      df.prediction <- ab.data.post %>% 
       inner_join(., ab.data.pre, 
                  by=c('subject_id', 'isotype_antigen')) %>% 
       dplyr::group_by(isotype_antigen) %>%
       dplyr::mutate(rank.titer = dense_rank(dplyr::desc(post.titer)),
                     titer.class = ifelse(post.titer >= quantile(post.titer, 
                                                                 0.5), 1, 0),
                     titer.fc = post.titer / pre.titer,
                     rank.titer.fc = dense_rank(dplyr::desc(titer.fc)),
                     titer.fc.class=ifelse(titer.fc >= quantile(titer.fc, 0.5), 
                                           1, 0)) %>%
       full_join(., pred, by='subject_id') %>% 
       dplyr::rename(predictor = predictor) %>% 
       dplyr::mutate(rank.predictor = dense_rank(dplyr::desc(predictor)))
      
      # Isotype-antigens pair of interest can be seelct
      if(!is.null(filter.ab.ag)){
         df.prediction <- df.prediction[(df.prediction$isotype_antigen 
                                         %in% filter.ab.ag), ]
      }
      
   } else {
      
      df.prediction <- pred %>% 
         dplyr::rename(predictor = predictor) %>% 
         dplyr::mutate(rank.predictor = dense_rank(dplyr::desc(predictor)))
   }
    
    return(df.prediction)
}

geneModuleCalc <- function(gene.ex.data, score, post.ab.data=NULL, 
                           pre.ab.data=NULL, module=NULL, filter.ab.ag=NULL, 
                           verbose=FALSE){

    if(score == 'qusage'){
       # Load gene modules
       load('../Study-1-Avey-20117/geneSetDB.rda', envir = .GlobalEnv)

       # Prep GeneSetDB per original code.  Can also use GSEAbase package.
       # Convert from vector of strings to a list
       gene.set <- strsplit(geneSetDB, "\t")      
       module.names <- as.data.frame(x = sapply(gene.set,"[", 1))
       colnames(module.names) <- 'module'
       gene.set <- lapply(gene.set, "[",-1:-2) # remove name and description columns
       gene.set <- lapply(gene.set, function(x){ x[ which( x != "") ] }) # remove empty strings
       names(gene.set) <- module.names$module # adding module names to gene_set
       
       # Selected gene module
       gene.module <- gene.set[[module]]
       
       # Generate the gene module expression matrix
       module.exp <- selectGeneExp(exp.mat=gene.ex.data, gene.set=gene.module, 
                                   verbose=verbose)
       
       # Calculate the average of the gene module, which is termed the gene 
       # module intensity in the original manuscript
       pred <- data.frame(subject_id = colnames(module.exp), 
                          predictor = colMeans(module.exp, na.rm = T))
       pred$subject_id <- as.numeric(pred$subject_id)

   }else if(score == 'bioage'){
      
      # Calculate the BioAge
      pred <- calculateBioAge(ex.mat=gene.ex.data, 
                              gene.map=gene.map, 
                              gene.synonym=gene.synonym,
                              predictive.modules=module,
                              verbose=verbose)
   } else {
      stop("The selected score is not valid. Please select one of the 
                following: 'qusage' or 'bioage'")
   }
   
   if(!(is.null(pre.ab.data) & is.null(post.ab.data))){
      # Subset antibody titers data
      ab.mask <- post.ab.data$subject_id %in% colnames(gene.ex.data)
      ab.mask.pre <- pre.ab.data$subject_id %in% colnames(gene.ex.data)
       
      # Subset antibody titer data
      ab.data.post <- post.ab.data[ab.mask, ] %>% 
         dplyr::select(subject_id, isotype_antigen, MFI_normalised) %>% 
         dplyr::rename(post.titer = MFI_normalised)
       
      ab.data.pre <- pre.ab.data[ab.mask.pre, ] %>% 
         dplyr::select(subject_id, isotype_antigen, MFI_normalised) %>% 
         dplyr::rename(pre.titer = MFI_normalised)
       
      # Create dataframe with prediction and response classes
      df.prediction <- ab.data.post %>% 
         inner_join(., ab.data.pre, 
                    by=c('subject_id', 'isotype_antigen')) %>% 
         dplyr::group_by(isotype_antigen) %>%
         dplyr::mutate(rank.titer = dense_rank(dplyr::desc(post.titer)),
                       titer.class = ifelse(post.titer >= quantile(post.titer, 
                                                                   0.5), 1, 0),
                       titer.fc = post.titer / pre.titer,
                       rank.titer.fc = dense_rank(dplyr::desc(titer.fc)),
                       titer.fc.class=ifelse(titer.fc >= quantile(titer.fc, 0.5), 
                                             1, 0)) %>%
         full_join(., pred, by='subject_id') %>% 
         dplyr::rename(predictor = predictor) %>% 
         dplyr::mutate(rank.predictor = dense_rank(dplyr::desc(predictor))) 
      
      # Isotype-antigens pair of interest can be seelct
      if(!is.null(filter.ab.ag)){
         df.prediction <- df.prediction[(df.prediction$isotype_antigen 
                                         %in% filter.ab.ag), ]
      }
      
   } else {
      
      df.prediction <- pred %>% 
         dplyr::rename(predictor = predictor) %>% 
         dplyr::mutate(rank.predictor = dense_rank(dplyr::desc(predictor)))
      
   }
    
    return(df.prediction)
    
}

mlModelPredict <- function(gene.ex.data=NULL, model, study, post.ab.data=NULL, 
                           pre.ab.data=NULL, cytof.data=NULL, populations=NULL, 
                           filter.ab.ag=NULL, paper.sig=FALSE, verbose=FALSE,
                           meta.data=NULL, year="2020"){
   
   if(study == 'Tsang_2014'){
      
      # Subset cell frequnecy data to the selected populations
      selected.pops <- cytof.data[cytof.data$cell_type_name %in% populations, ] %>% 
         dplyr::select(cell_type_name, subject_id, percent_live_cell) %>% 
         pivot_wider(id_cols = subject_id, names_from = cell_type_name,
                     values_from = percent_live_cell) %>% 
         column_to_rownames("subject_id") %>% 
         as.matrix()
      
      centroids <- model$centroids
      variances <- model$variances
      priors <- model$priors
      
      # DLDA prediction
      dist <- t(apply(selected.pops, 1, function(z) 
         colSums((z-centroids)^2/variances) + priors))
      dist <- dist - rowMeans(dist)
      prob <- safeexp(-0.5*dist)
      prob <- prob/rowSums(prob)
      probli <- prob[,2]
      
      # Change format of class probability predictions
      pred <- data.frame(probli) %>% 
         rownames_to_column('subject_id') %>% 
         dplyr::rename(predictor = probli) %>% 
         drop_na()
      
   pred$predictor <- as.numeric(pred$predictor)
      
      
   } else if(study == 'fourati_2015_NB'){
      
      ### Download and load EM131 ExpressionSet (cf. Fig. 2)
      esetPath <- file.path("https://storage.googleapis.com/em131_20120924",
                            "Fig3ab/EM131.eset.RData")
      
      esetBin <- getBinaryURL(url = esetPath, followlocation = TRUE)
      esetFile <- basename(esetPath)
      esetCon <- file(esetFile, open = "wb")
      writeBin(esetBin, esetCon)
      close(esetCon)
      
      # load RData file
      load(file = esetFile)
      
      # Probe ID to gene symbol 
      # probe.gene.map <- fData(eset)
      probe.gene.map <- read_excel("supplementary_table_4a.xlsx")
      
      
      # The 15-gene signature used as input
      print(names(model$"tables"))
      genes.mat <- probe.gene.map[match(names(model$"tables"),
                                    table = probe.gene.map$ProbeSetID), ]
      
      # Change gene symbol "HP HPR" to "HP"
      genes.mat[genes.mat$GeneSymbol == 'HP HPR', "GeneSymbol"] <- "HP"
      
      # Remove rows without gene symbols
      genes.mat <- genes.mat[genes.mat$GeneSymbol != '---', ]
      gene.sig <- genes.mat$GeneSymbol
      
      # Generate the gene signature expression matrix
      gene.sig.exp <- selectGeneExp(exp.mat=gene.ex.data, gene.set=gene.sig,
                                    verbose=verbose)
      
      # Map Ensembl IDs to gene symbols
      rownames(gene.sig.exp) <- mapvalues(rownames(gene.sig.exp),
                                          from=gene.map$ensembl_gene_id,
                                          to=gene.map$external_gene_name,
                                          warn_missing=FALSE)

      # Map gene symbols to probe IDs
      rownames(gene.sig.exp) <- mapvalues(rownames(gene.sig.exp),
                                         from=genes.mat$GeneSymbol,
                                         to=genes.mat$ProbeSetID,
                                         warn_missing=FALSE)
      
      
      # Scale gene expression as in original study
      gene.sig.exp <- t(scale(t(gene.sig.exp)))
      
      # Predict the class probabilities with the Naive Bayes classifier
      pp <- predict(object = model, newdata = t(gene.sig.exp), type = "raw")
      
      # Use ratio between post-probalities of being responders on non-responders
      # as predictor of response of to vaccine
      pred <- apply(log(pp), MARGIN = 1, FUN = diff)
      pred <- data.frame(cbind(colnames(gene.sig.exp), pred)) %>% 
          dplyr::rename(subject_id = V1,
                        predictor = pred)
      pred$predictor <- as.numeric(pred$predictor)
      
   } else if(study =='fourati_2015_LR'){
      
      # Subset cell frequnecy data to the selected populations
      df.cells <- cytof.data[cytof.data$cell_type_name %in% populations, ] %>% 
         dplyr::select(cell_type_name, subject_id, percent_live_cell) %>% 
         pivot_wider(id_cols = subject_id, names_from = cell_type_name,
                     values_from = percent_live_cell)
      
      if(year=="2020"){
         df.cells <- df.cells %>% 
            dplyr::rename("B cells.PCT IgG+ memory B cells" = Bcells,
                          "B cells.PCT switched IgG+ MK memory B cells" = `ASCs (Plasmablasts)`,
                          "T cells.PCT TEM2 in CD4" = TemCD4,
                          "Innate cells.MdFI CD40 in PDC" = pDC)
         
      }else if(year=="2021"){
         
         df.cells <- df.cells %>% 
            dplyr::rename("B cells.PCT IgG+ memory B cells" = Bcells,
                          "B cells.PCT switched IgG+ MK memory B cells" = `Memory B cells`,
                          "T cells.PCT TEM2 in CD4" = TemCD4,
                          "Innate cells.MdFI CD40 in PDC" = pDC)
         
      }
      
      
      # Predict on the scale of the additive predictors 
      pred <- predict(object=fit, newdata = df.cells, type = "link")
      pred <- data.frame(cbind(df.cells$subject_id, pred)) %>% 
         dplyr::rename(subject_id = V1,
                       predictor = pred) %>% 
         drop_na()
      pred$predictor <- as.numeric(pred$predictor)
      
   } else if(study == 'fourati_2021'){

      # Get the top 500 varying genes from model
      top500 <- colnames(model_list$ptype)
      
      # Identify the top 500 varying genes in the gene expression matrix
      gene.sig.exp <- selectGeneExp(exp.mat=gene.ex.data, gene.set=top500, 
                                    verbose=verbose)
      
      rownames(gene.sig.exp) <- mapvalues(rownames(gene.sig.exp),
                                          from=gene.map$ensembl_gene_id,
                                          to=gene.map$external_gene_name,
                                          warn_missing=FALSE)
      
      # Adjust format gene expression data to fit with Random Forest classifier
      gene.sig.matrix <- t(gene.sig.exp) %>%
         as.matrix()

      pred <- predict(model, newdata=gene.sig.matrix, type="prob") %>% 
           dplyr::rename(predictor = highResponder) %>% 
           rownames_to_column('subject_id')
       
   } else if (study == "Furman_2013") {
      
      # Subset meta data
      meta.day0 <- meta.data %>% 
         filter(planned_day_relative_to_boost == 0)
      
      # Calulcate the age of the subject and round down to nearest integer
      meta.day0$age <- as.numeric(difftime(meta.day0$date_of_boost,
                                           meta.day0$year_of_birth, 
                                           units = "weeks"))/52.25
      # meta.day0$age_rounded <- floor(meta.day0$age)
      
      pred <- meta.day0 %>% 
         dplyr::select(subject_id, age) %>% 
         dplyr::rename(predictor = age)
      
      
   } else if(study == 'Bartholomeus_2018'){
      
      # Transpose dataframe
      gene.ex.data.T <- as.data.frame(t(gene.ex.data)) %>% 
         rownames_to_column("subject_id")
      
       # Remove genes with zero variance, as PCA cannot handle this
       gene.ex.data.filtered <- Filter(var, gene.ex.data.T) %>% 
          column_to_rownames("subject_id")
       
       # Compress features with PCA
       pcares <- prcomp(gene.ex.data.filtered, center=TRUE, scale=TRUE)
       PCAs <- data.frame(pcares$x[,1:5])
       
       pred <- predict(model, newdata=PCAs, usekernel=TRUE, 
                       prior=c(20/34, 14/34))
       pred <- pred$posterior[,2]
       pred <- as.data.frame(pred) %>% 
          dplyr::rename(predictor = pred) %>% 
          rownames_to_column('subject_id')
       
       pred$predictor <- as.numeric(pred$predictor)
   }
   
   pred$subject_id <- as.numeric(pred$subject_id)
   
   if(!(is.null(pre.ab.data) & is.null(post.ab.data))){
   
      if(!is.null(gene.ex.data)){
         # Subset antibody titers data
         ab.mask <- post.ab.data$subject_id %in% colnames(gene.ex.data)
         ab.mask.pre <- pre.ab.data$subject_id %in% colnames(gene.ex.data)
         
      } else if (!is.null(cytof.data)){
         # Subset antibody titers data
         ab.mask <- post.ab.data$subject_id %in% pred$subject_id
         ab.mask.pre <- pre.ab.data$subject_id %in% pred$subject_id
         
      } else {
         # Subset antibody titers data
         ab.mask <- post.ab.data$subject_id %in% pred$subject_id
         ab.mask.pre <- pre.ab.data$subject_id %in% pred$subject_id
         
      }
      
      # Subset antibody titer data
      ab.data.post <- post.ab.data[ab.mask, ] %>% 
         dplyr::select(subject_id, isotype_antigen, MFI_normalised) %>% 
         dplyr::rename(post.titer = MFI_normalised)
      
      ab.data.pre <- pre.ab.data[ab.mask.pre, ] %>% 
         dplyr::select(subject_id, isotype_antigen, MFI_normalised) %>% 
         dplyr::rename(pre.titer = MFI_normalised)
      
      # Create dataframe with prediction and response classes
      df.prediction <- ab.data.post %>% 
         inner_join(., ab.data.pre, 
                    by=c('subject_id', 'isotype_antigen')) %>% 
         dplyr::group_by(isotype_antigen) %>%
         dplyr::mutate(rank.titer = dense_rank(dplyr::desc(post.titer)),
                       titer.class = ifelse(post.titer >= quantile(post.titer, 
                                                                   0.5), 1, 0),
                       titer.fc = post.titer / pre.titer,
                       rank.titer.fc = dense_rank(dplyr::desc(titer.fc)),
                       titer.fc.class=ifelse(titer.fc >= quantile(titer.fc, 0.5), 
                                             1, 0)) %>%
         full_join(., pred, by='subject_id') %>% 
         dplyr::rename(predictor = predictor) %>% 
         dplyr::mutate(rank.predictor = dense_rank(dplyr::desc(predictor)))
      
      # Isotype-antigens pair of interest can be seelct
      if(!is.null(filter.ab.ag)){
         df.prediction <- df.prediction[(df.prediction$isotype_antigen 
                                         %in% filter.ab.ag), ]
      }
   } else {
      
      df.prediction <- pred %>% 
         dplyr::rename(predictor = predictor) %>% 
         dplyr::mutate(rank.predictor = dense_rank(dplyr::desc(predictor)))
      
   }
   
   return(df.prediction)
}

# Function for evaluating the predictive power of a model on the CMI-PB dataset
evaluateModel <- function(df.prediction,
                         add.to.title="", N.perm=1000, 
                         seed.value=123, direction="<", 
                         plot.corr=TRUE,  plot.auc=TRUE, save.plot=TRUE,
                         score.label='Gene signature score',
                         fn1='pred_vs_day14_heatmap',
                         fn2='pred_vs_day14_correlation ',
                         fn3='pred_vs_FC_heatmap',
                         fn4='pred_vs_FC_correlation',
                         fn5='pred_vs_day14_AUC',
                         fn6='pred_vs_day14_boxplot',
                         fn7='pred_vs_FC_AUC',
                         fn8='pred_vs_FC_boxplot',
                         path='Results/'){
    
    # Initialize dataframes for storing the Spearman's correlations and AUC 
    # scores
    df.metrics <- data.frame()
    df.roc <- data.frame()
    df.roc.fc <- data.frame()
    
    # Loop over every antibody-antigen pair and calculate the Spearman's 
    # correlation coefficient
    for(ab.ag in unique(df.prediction$isotype_antigen)){
        
        # Filter the data by antibody-antigen pair predictor and response 
        # variable in dataframe
        df.pred.tmp <- df.prediction %>%
            dplyr::filter(isotype_antigen == ab.ag)
        
        # Calculate correlation between predictor and day 14 titers - evaluate 
        # with N permutations
        corr.results <- corrEval(df.pred.tmp, response='post.titer', 
                                 score='predictor', seed.value=seed.value, 
                                 N.perm=N.perm, direction=direction, 
                                 ab.ag=ab.ag)
        
        corr <- corr.results$corr
        corr.p <- corr.results$corr.p
        
        # Calculate correlation  between predictor and fold-change in antibody 
        # titers - evaluate with N permutations
        corr.fc.results <- corrEval(df.pred.tmp, response='titer.fc', 
                                 score='predictor', seed.value=seed.value, 
                                 N.perm=N.perm, direction=direction, 
                                 ab.ag=ab.ag)
        
        corr.fc <- corr.fc.results$corr
        corr.fc.p <- corr.fc.results$corr.p
        
        # Evaluation the gene signature with ROC curves and AUC scores - 
        # class defined by top and bottom half of the antibody response
        roc.results <- rocEval(df.pred.tmp, class='titer.class', 
                               score='predictor', seed.value=seed.value, 
                               N.perm=N.perm, ab.ag=ab.ag, direction=direction)
        
        r <- roc.results$r
        r.p <- roc.results$r.p.value
        r.df <- roc.results$r.df
        df.roc <- rbind(df.roc, r.df)
        
        # Evaluation the gene signature with ROC curves and AUC scores
        # class defined by standard deviation increase
        roc.fc.results <- rocEval(df.pred.tmp, 
                                  class='titer.fc.class', 
                                  score='predictor', 
                                  seed.value=seed.value, 
                                  N.perm=N.perm, ab.ag=ab.ag, 
                                  direction=direction)
        
        r.fc <- roc.fc.results$r
        r.p.fc <- roc.fc.results$r.p.value
        r.df.fc <- roc.fc.results$r.df
        df.roc.fc <- rbind(df.roc.fc, r.df.fc)

        # Collecting the Spearman's correlation coefficent and AUC as well as p-
        # values
        tmp <- data.frame(isotype_antigen = ab.ag, 
                          corr.value = corr$estimate,
                          label.corr = sprintf("r = %.2g", corr$estimate),
                          corr.p.value = corr.p,
                          label.corr.pv = sprintf("p = %.2g", corr.p),
                          corr.fc.value = corr.fc$estimate,
                          label.corr.fc = sprintf("r = %.2g", corr.fc$estimate),
                          corr.fc.p.value = corr.fc.p,
                          label.corr.fc.pv = sprintf("p = %.2g", corr.fc.p),
                          auc.score = r$auc, 
                          label.auc = sprintf("AUC = %.2f", r$auc),
                          auc.p.value = r.p,
                          label.auc.pv = sprintf("p = %.2g", r.p),
                          auc.score.fc= r.fc$auc, 
                          label.auc.fc = sprintf("AUC = %.2f", 
                                                  r.fc$auc), 
                          auc.p.value.fc = r.p.fc,
                          label.auc.pv.fc = sprintf("p = %.2g", 
                                                     r.p.fc))
        df.metrics = rbind(df.metrics, tmp)
    }
    
    # Save results
    df.results <- df.prediction %>% 
        inner_join(., df.metrics, by='isotype_antigen') %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(titer.class.label = ifelse(titer.class == 1, 
                                                 'Positive', 'Negative'),
                      titer.fc.class.label = ifelse(titer.fc.class == 1, 
                                                     'Positive', 'Negative'),
                      stars.corr.raw = cut(corr.p.value, breaks=c(-Inf, 0.001, 
                                                                  0.01, .05, Inf), 
                                  label=c("***", "**", "*", "")),
                      stars.corr.fc = cut(corr.fc.p.value, breaks=c(-Inf, 0.001, 
                                                                    0.01, 0.05, Inf), 
                                      label=c("***", "**", "*", "")),
                      stars.auc.raw = cut(auc.p.value, breaks=c(-Inf, 0.001, 
                                                                  0.01, .05, Inf), 
                                           label=c("***", "**", "*", "")),
                      stars.auc.fc = cut(auc.p.value.fc, breaks=c(-Inf, 0.001, 
                                                                  0.01, .05, Inf), 
                                           label=c("***", "**", "*", ""))) 
    
    # sig_p = ifelse(corr.p.value < .05, T, F), 
    # p_if_sig = ifelse(corr.p.value <.05, corr.p.value, NA), 
    # r_if_sig = ifelse(corr.p.value <.05, corr.value, NA),
    
    if(plot.corr){
       
       ############## RAW DAY 14 TITERS ############## 
       
       # Title and subtitle the ROC plots
       corr.title <- paste("Correlations between", score.label, 
                           "and day 14 titer")
       
       heat.subtitle <- paste("Stars indicate significant Spearman's correlation 
       coefficients", add.to.title)
       
       # Plot the correlations using a heatmap
       heatmapPlotter(data=df.results, score='corr.value', save.plot=save.plot, 
                      path=path, filename=fn1, main.title=corr.title, 
                      stars="stars.corr.raw",
                      sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                     x=heat.subtitle))
       
       # Subtitle for correlation plots
       corr.subtitle <- "P-values and correlation coefficients are based on 
       Spearman's rank correlation test"
       
       # Plot the correlations using scatter plots
       corrPlotter(data=df.results, x.var='rank.predictor', y.var='rank.titer', 
                   save.plot=save.plot, 
                   path=path, 
                   filename=fn2, 
                   main.title=corr.title, 
                   df.metrics = df.metrics, 
                   label.corr='label.corr', 
                   label.corr.pv='label.corr.pv', 
                   xlabel=paste(score.label, "- ranked"), 
                   ylabel='Titer values - ranked', 
                   sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                  x=corr.subtitle))
       
       ############## FOLD-CHANGE ############## 
       # Title and subtitle the ROC plots
       corr.fc.title <- paste("Correlations between", score.label, 
                           "and fold-change")
       
       # Plot the correlations using a heatmap
       heatmapPlotter(data=df.results, score='corr.fc.value', 
                      save.plot=save.plot, path=path, filename=fn3, 
                      main.title=corr.fc.title, stars="stars.corr.fc",
                      sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                     x=heat.subtitle))
       
       # Plot the correlations using scatter plots
       corrPlotter(data=df.results, x.var='rank.predictor', y.var='rank.titer.fc', 
                   save.plot=save.plot, path=path, filename=fn4, 
                   main.title=corr.fc.title, df.metrics, 
                   label.corr='label.corr.fc', 
                   label.corr.pv='label.corr.fc.pv', 
                   xlabel=paste(score.label, "- ranked"), 
                   ylabel='Fold-changer titers (day 14 / day 0) - ranked', 
                   sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                  x=corr.subtitle))
    }
    
    if(plot.auc){
        
       
       ############## RAW DAY 14 TITERS ############## 
       
       # Title and subtitle the ROC plots
       roc.title <- paste('AUC score and ROC curve for each antibody-antigen 
                          pair', add.to.title)
       roc.subtitle <- "Top half of observed titers on day 14 are considered 
       'positive' and the lower half are considered 'negative'"
       
       # Generate ROC curve when defining the responders from top and bottom 
       # half day 14 antibody titers
       rocPlotter(df.roc, spec='Specificity', sens='Sensitivity',  
                  save.plot=save.plot, path=path, filename=fn5, 
                  df.metrics = df.metrics, label.auc='label.auc', 
                  label.auc.pv='label.auc.pv',
                  main.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                 x=roc.title),
                  sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                 x=roc.subtitle))
       
       # Generate boxplots when defining the responders top and bottom half
       boxplot.title <- paste("Comparing the", score.label, "between 'Negative' 
       and 'Positive' subjects", add.to.title)
       boxPlotter(data=df.results, response.class='titer.class.label', 
                  score='predictor', save.plot=save.plot, path=path, 
                  filename=fn6, xlabel="Response class", ylabel=score.label, 
                  main.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                  x=boxplot.title),
                  sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                  x=roc.subtitle))
       
       
       ############## FOLD-CHANGE ############## 
           
       # Sub-title when comparing the pre- and post-vaccination antibody
       # titers
       subtitle.fc <- "Top half of observed fold-change titer values (day 14 /
       day 0) are considered 'positive' and the lower half are considered 
        'negative'"
       
       # Generate ROC curve when defining the responders based on standard
       # deviation increase in antibody titers
       rocPlotter(df.roc.fc, spec='Specificity', sens='Sensitivity', 
                  path=path, filename=fn7, df.metrics = df.metrics,
                  label.auc='label.auc.fc', save.plot=save.plot,
                  label.auc.pv='label.auc.pv.fc',
                  main.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                  x=roc.title),
                  sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                 x=subtitle.fc))
       
       # Generate boxplots when defining the responders based on standard
       # deviation increase in antibody titers
       boxPlotter(data=df.results, response.class='titer.fc.class.label', 
                  score='predictor', save.plot=save.plot,
                  path=path, filename=fn8,
                  main.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                  x=boxplot.title), 
                  sub.title=gsub(pattern='\\s{2,}', replacement=" ", 
                                 x=subtitle.fc), 
                  xlabel="Response class", ylabel=score.label)
    }
    
    out <- list()
    
    out$results <- df.results
    out$roc <- df.roc
    out$roc.FC <- df.roc.fc
    
    return(out)
}

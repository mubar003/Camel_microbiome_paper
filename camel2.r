# Camel microbiome code V3 ------------------------------------------------
# GAA FAM 231203

setwd("/gdat/fathi/camel/")

# Pre-parsing code (skip this section if doing analysis) ------------------

# Basic QC
qual = read.delim("prokaryotes.summary.tsv",row=1)

par(mfrow=c(1,3))
for (metric in (c("Completeness","Score","Contamination"))) {
  ho=hist(qual[,metric],main=metric,xlab="Percent")
  summary(qual[,metric])
  med=median(ho$counts)
  segments(0,med,9999999,med,col='red',lty=2)
}

setwd("/gdat/fathi/camel/prof/")
xcov = data.frame(); bcov = data.frame()
files = dir(pattern = "*.cov")
for (file in files) for (file in files) {
  t = read.table(file,T,row=1,as.is=T,check=F,comment='')
  xcov[rownames(t),gsub(".cov.*","",file)] = 
    t$XCov_adamantium
  bcov[rownames(t),gsub(".cov.*","",file)] = 
    t$Adamantium_prop
}
xcov[is.na(xcov)]=0
bcov[is.na(bcov)]=0
setwd("..")

# Save the matrix
cat("Bug\t",file="CamCovs.tsv"); write.table(xcov,"CamCovs.tsv",T,F,'\t')
cat("Bug\t",file="CamBCovs.tsv"); write.table(bcov,"CamBCovs.tsv",T,F,'\t')


# Begin stats processing here! --------------------------------------------

# Load required data
setwd("/gdat/fathi/camel/")
xcov = read.delim("CamCovs.tsv",row=1, check=F)
bcov = read.delim("CamBCovs.tsv",row=1, check=F)

# Visualize the data (commented out for efficiency)
#par(mfrow=c(8,8)); for (i in 1:nrow(xcov)) hist(log10(as.numeric(xcov[i,])),main=rownames(xcov)[i])

# Massage the matrix into log10 Relative Abundance
pseudo = min(xcov[xcov > 0]) * .5
ra = xcov; ra[ra==0] = pseudo
ra = log10(sweep(ra,2,colSums(ra),'/'))

# Read in the metadata
camdat = read.delim("Metadata_V2_GAA_FM.tsv",row=1)
ra = ra[,rownames(camdat)]

# Compute alpha diversity (# species approximated by # 95% MAGs with over 0.25X coverage)
adiv = colSums(bcov[,rownames(camdat)] >= .25)
hist(adiv,15)

# Add alpha diversity to metadata
camdat$Alpha_diversity = adiv


# Summary stats and demographics ------------------------------------------

# Function to create summary row for numerical variables
add_summary <- function(data, var_name) {
  return(data.frame(
    Variable = var_name,
    Min = min(data),
    'Q1' = quantile(data, 0.25),
    Median = median(data),
    Mean = round(mean(data),2),
    'Q3' = quantile(data, 0.75),
    Max = max(data)
  ))
}

# Numerical variable summaries
numerical_summaries <- rbind(
  add_summary(camdat$Age_y, "Age"),
  add_summary(camdat$Weight_factor_1_5, "Weight_Factor"),
  add_summary(camdat$Diet_diversity_factor_1_4, "Diet_Diversity"),
  add_summary(camdat$Bristol_factor_1_7, "Bristol_Factor"),
  add_summary(camdat$Alpha_diversity, "Alpha_Diversity")
)

# Write out the table for numerical summaries as TSV
write.table(numerical_summaries, file = "Numerical_Summaries_Table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Function to create a data frame for categorical count summaries
add_categorical_summary <- function(data, var_name) {
  cat_summary <- as.data.frame(table(data))
  names(cat_summary) <- c("Category", "Count")
  cat_summary$Variable <- var_name
  cat_summary <- cat_summary[, c("Variable", "Category", "Count")]
  return(cat_summary)
}

# Categorical variable count summaries
categorical_summaries <- rbind(
  add_categorical_summary(camdat$Gender, "Gender"),
  add_categorical_summary(camdat$Pregnancy, "Pregnancy"),
  add_categorical_summary(camdat$Disease, "Disease_Status"),
  add_categorical_summary(camdat$Living.condition, "Living_Condition"),
  add_categorical_summary(camdat$Diet_grass, "Diet_grass"),
  add_categorical_summary(camdat$Diet_bread, "Diet_bread"),
  add_categorical_summary(camdat$Diet_barley, "Diet_barley"),
  add_categorical_summary(camdat$Diet_milk, "Diet_milk"),
  add_categorical_summary(camdat$Diet_wheat, "Diet_wheat"),
  add_categorical_summary(camdat$Diet_supplement, "Diet_supplement")
  
)

# Adjust the row names
row.names(categorical_summaries) <- NULL

# Write out the table for categorical summaries as TSV
write.table(categorical_summaries, file = "Categorical_Summaries_Table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# PCA visualization -------------------------------------------------------

# Create a principle components plot (beta diversity)
pca = prcomp(ra)
par(mfrow=c(1,1))  # pdf export at 10x7 works well
plot(pca$rotation[,c(1,2)],pch=19,col="#55555588",cex=2)

# Layer on colors
library(RColorBrewer)

plot_pca_with_colors <- function(pca, metadata, color_var, palette_name, transparency = 0.75, scale = 2, 
                                 withheld_samples = NULL, sub="") {
  # Check if withheld_samples is provided and not NULL
  if (!is.null(withheld_samples)) {
    # Ensure withheld_samples is a character vector
    withheld_samples <- as.character(withheld_samples)
    # Exclude the withheld samples from the PCA and metadata
    valid_samples <- !rownames(metadata) %in% withheld_samples
    pca$x <- pca$x[valid_samples, ]
    metadata <- metadata[valid_samples, ]
  }
  
  # Determine if the variable is continuous or discrete based on the palette name
  is_continuous <- palette_name == "gradient"
  
  # Automatically determine the number of colors based on the unique values of color_var
  if (is_continuous) {
    num_colors <- length(unique(metadata[[color_var]]))
    color_palette <- colorRampPalette(c('darkgoldenrod1', 'blue'))(num_colors)
  } else {
    num_colors <- length(unique(metadata[[color_var]]))
    max_colors <- ifelse(palette_name %in% rownames(brewer.pal.info), brewer.pal.info[palette_name, "maxcolors"], num_colors)
    color_palette <- brewer.pal(min(max_colors, num_colors), palette_name)
    if (num_colors > max_colors) {
      warning("Number of unique values exceeds the number of available colors in the palette. Colors will be recycled.")
    }
  }
  
  # Ensure the color variable is a factor and map it to colors
  color_factor <- as.factor(metadata[[color_var]])
  colors_mapped <- color_palette[(as.numeric(color_factor) - 1) %% length(color_palette) + 1]
  
  # Add transparency to colors
  colors_transparent <- adjustcolor(colors_mapped, alpha.f = transparency)
  
  # Plot PCA with the specified colors and scale
  plot(pca$x[, c(1, 2)], pch = 19, col = colors_transparent, cex = scale, 
       xlab = paste("PC1 (", round(pca$sdev[1] / sum(pca$sdev) * 100), "%)", sep=""),
       ylab = paste("PC2 (", round(pca$sdev[2] / sum(pca$sdev) * 100), "%)", sep=""),
       main = paste("Var: ", color_var),sub=sub)
  
  # Add a legend conditionally
  if (is_continuous) {
    if (num_colors > 10) {
      # For continuous variables, show the range
      legend_values <- range(metadata[[color_var]], na.rm = TRUE)
      legend_labels <- c(format(legend_values[1], digits = 2), format(mean(legend_values), digits = 2), format(legend_values[2], digits = 2))
      legend_colors <- color_palette[c(1, num_colors %/% 2, num_colors)]
    } else {
      # For continuous variables with 10 or fewer unique values, show all values
      legend_values <- sort(unique(metadata[[color_var]]))
      legend_labels <- format(legend_values, digits = 2)
      legend_colors <- color_palette[seq_along(legend_values)]
    }
  } else if (num_colors <= 10) {
    # For discrete variables with 10 or fewer levels, show all levels
    legend_values <- levels(color_factor)
    legend_labels <- legend_values
    legend_colors <- color_palette[1:num_colors]
  } else {
    # Do not display a legend if the conditions are not met
    return()
  }
  
  # Display the legend
  legend("bottomright", legend=legend_labels, col=legend_colors, pch=19, 
         pt.cex = 1.5, cex = 0.7, bty = "n", y.intersp = 1.2, x.intersp = 1.2, 
         text.col = legend_colors, inset = c(0.02, 0.02), 
         box.col = NA, bg = adjustcolor("white", alpha.f = 0))
}

# Perform PCA on the t(ra) since features are rows and samples are columns
pca <- prcomp(t(ra), scale. = TRUE)

plot_pca_with_colors(pca, camdat, "Herd_ID", "Set3",withheld_samples = NULL)

# Set up the plot grid
pdf("PCAsO.pdf",11,8)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

# Plotting the PCAs using the function
outliers = "42A"; # NULL# Set to NULL to produce the original. But we want to exclude camel "42A" only.
plot_pca_with_colors(pca, camdat, "Herd_ID", "Set3",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Weight_factor_1_5", "gradient",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Disease", "Set1",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Bristol_factor_1_7", "gradient",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Stool_darkness_factor_1_4", "gradient",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Alpha_diversity", "gradient",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Age_y", "gradient",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Diet_diversity_factor_1_4", "gradient",withheld_samples = outliers)
plot_pca_with_colors(pca, camdat, "Captivity_factor_1_3", "gradient",withheld_samples = outliers)

# Reset par to default
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()

library(vegan)

# Loop through all metadata variables and perform adonis2 tests using the distance matrix (not ra directly, which is log10 transformed)
vars = c("Herd_ID","Age_y","Gender","Pregnancy","Color", "Weight_factor_1_5", "Diet_grass", "Diet_wheat", "Diet_barley", "Diet_bread", "Diet_milk", "Diet_supplement", 
         "Diet_diversity_factor_1_4","Grazing_factor_0_2", "Num_cohoused", "Time_day_factor_1_3", "Temperature_C","Stool_darkness_factor_1_4", "Bristol_factor_1_7", 
         "Disease", "Captivity_factor_1_3", "Alpha_diversity")
adonis_results <- data.frame(row.names = vars)
adonis_results$P_Value = 1
dm = dist(t(ra)) # distance matrix (euclidean because we're using log transformed RA, not counts/simple RA)
set.seed(123)
for (var in vars) {
  ad_test = adonis2(dm ~ camdat[,var], permutations = 9999)
  adonis_results[var,"P_Value"] = ad_test$`Pr(>F)`[1]
  adonis_results[var,"R2"] = ad_test$R2[1]
}
adonis_results = adonis_results[order(adonis_results$P_Value,-adonis_results$R2),]

# Write out the table for adonis results as TSV
write.table(adonis_results, file = "Adonis_Results_Table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# For all significant variables, do a plot_pca_with_colors on just those < p = 0.05 in adonis_results (plus disease if it wasn't included!)
# Filter out variables that are uninteresting (high p value or too few points in class)
sig_vars = rownames(adonis_results)[adonis_results$P_Value < 0.05]
#if (!"Disease" %in% sig_vars) sig_vars = c(sig_vars,"Disease") # in case you want to add disease in; but insig.
sig_vars=sig_vars[sig_vars!="Diet_grass"]

pdf("PCAsO_SigVars.pdf",11,8)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1),mgp=c(2,1,0))
for (var in sig_vars) {
  plot_pca_with_colors(pca, camdat, var, "gradient",withheld_samples = outliers,
                       sub=paste0("P=",adonis_results[var,"P_Value"],", R2=",
                                  signif(adonis_results[var,"R2"],3)))
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()


# Heatmap visualization ---------------------------------------------------

pcor = cor(cbind(Alpha_diversity=adiv,camdat[,c(grep("Age|factor|cohoused|Temp",colnames(camdat)))]),use = "pairwise.complete.obs",method="sp")
pheatmap::pheatmap(pcor)

# Differential abundance --------------------------------------------------

qual = read.delim("prokaryote_summary.csv", row =1)
rownames(qual) = sub(".fa$","",rownames(qual))
qual$genus = sub(";s__.*","",qual$classification)
qual$genus = sub(".*;g__","g__",qual$genus)
qual$genusS = ifelse(qual$genus=="g__" | qual$genus=="Unclassified",paste0("g__",rownames(qual)),qual$genus)

# alternative genus collapse
gensplt = strsplit(read.table("genera_maybe_K2837.reps",F,';','')[,1],'\t',T)
genmap = data.frame()
for (clus in gensplt) for (memb in clus) 
  genmap[memb,"Centroid"]=clus[[1]][1]

qual$genus = ifelse(qual$genus=="g__" | grepl("unclassified",qual$genus,T), 
                    paste0("g__",genmap[rownames(qual),1]),qual$genus)

#View(cbind(qual[,c(1,24,25)],Eq=qual$genus==qual$genusS)) # visualize the new assignments
qual = qual[rownames(xcov),] # Sync by genome name first
xcov.g = rowsum(xcov,qual$genus) # not rowSums()! Here we create a new collapsed table by sum, like aggregate(x,by=list(qual$genus),FUN=sum)

# Get the relative abundance log10 conversion logic from above
pseudo = min(xcov.g[xcov.g > 0]) * .5
ra.g = xcov.g; ra.g[ra.g==0] = pseudo
ra.g = log10(sweep(ra.g,2,colSums(ra.g),'/'))

ra.g = ra.g[,rownames(camdat)] # sync to metadata order too

bcor.g = cor(cbind(Num_species=adiv,camdat[,c(grep("Age|factor|cohoused|Temp",colnames(camdat)))]),t(ra.g),use = "pairwise.complete.obs",method="sp")
#pheatmap::pheatmap(bcor.g) # too messy to display


# Machine learning analysis -----------------------------------------------


# Let's use the genus variables from now on, for machine learning and other things.
library(randomForest)
library(ROCR)


# Let's extend this type of analysis to all variable names in "vars":
# - for each var, if it has just 2 distinct values cast as factor (classification), else regression (store indicator var for this)
# - run random forest basically as above.
# - store the AUC if classification or % var explained if regression
# - determine the top 2 top predictors using variable importance measure in rf: store them, and:
#   * for each of these top 2 predictors, store lm coefficients (if regression) or glm/logistic coefs (if classification), and p-value
# Save final table with columns: variable, regression_or_classification, auc_or_varExplained, topPred1, topPred1_coef, topPred1_pval, topPred2, topPred2_coef, topPred2_pval

# Sort the variable names by how many levels they have
varsP = vars[!vars %in% c("Color","Herd_ID")]
varsP = names(sort(apply(camdat[,varsP],2,function(xx) length(unique(xx)))))

par(mfrow=c(10,6), mar = c(4, 4, 2, 1),mgp=c(2,1,0))
results = data.frame()
set.seed(123)
for (var in varsP) {
  outcome = camdat[,var]
  doClass = length(unique(camdat[,var]))==2
  if (doClass) outcome = factor(outcome)
  rf = randomForest(t(ra.g),outcome,ntree = 1000)
  varImps = rf$importance[order(rf$importance[,1],decreasing = T),1,drop=F]
  if (doClass) {
    pred_ROCR = prediction(rf$votes[,2], outcome)
    roc_ROCR = performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
    auc = performance(pred_ROCR, measure = "auc")@y.values[[1]]
    results[var,"Task"]="Classification"
    results[var,"Perf"]=auc
    
    plot(roc_ROCR, main = var, colorize = T,sub=paste("AUC =",round(auc,2))); abline(a = 0, b = 1)
    
    for (top in 1:2) { # consider the top 2 predictors by RF decNodeImpurity
      bug=rownames(varImps)[top]
      statsy = glm(outcome ~ as.numeric(ra.g[bug,]),family = 'binomial')
      sumstat=summary(statsy); 
      pval = sumstat$coefficients[2,4]; coeff = sumstat$coefficients[2,1]
      results[var,"in"]=levels(outcome)[2]
      results[var,paste0("Pred",top)]=bug
      results[var,paste0("Pred",top,"coef")]=coeff
      results[var,paste0("Pred",top,"pval")]=pval
      
      boxplot(as.numeric(ra.g[bug,]) ~ outcome,main=bug,ylab="",
              xlab=paste0("Imp=",round(varImps[top],2), "; P=",signif(pval,3),", Coef=",signif(coeff,3)) ) 
    }
  } else {
    
    results[var,"Task"]="Regression"
    results[var,"Perf"]=rf$rsq[length(rf$rsq)]
    
    plot(rf$predicted, outcome, col = rgb(red = abs(outcome - rf$predicted)/max(abs(outcome - rf$predicted)), 
         green = 0, blue = 0,alpha = .6), pch=19,main=var,ylab="Actual",xlab="Predicted",
         sub=paste0("VarExp = ",round(rf$rsq[length(rf$rsq)],2) )); abline(a = 0, b = 1,col="blue")
    #abline(lm(outcome~ rf$predicted)) #,col="blue",lwd=2)
    
    for (top in 1:2) { # consider the top 2 predictors by RF decNodeImpurity
      bug=rownames(varImps)[top]
      statsy = lm(as.numeric(ra.g[bug,]) ~ outcome)
      sumstat=summary(statsy); 
      pval = sumstat$coefficients[2,4]; coeff = sumstat$coefficients[2,1]
      #results[var,"in"]=levels(outcome)[2]
      results[var,paste0("Pred",top)]=bug
      results[var,paste0("Pred",top,"coef")]=coeff
      results[var,paste0("Pred",top,"pval")]=pval
      
      plot(as.numeric(ra.g[bug,]) ~ outcome,main=bug,ylab="",
              xlab=paste0("Imp=",round(varImps[top],2), "; P=",signif(pval,3),", Coef=",signif(coeff,3)) ) 
      abline(statsy)
    }
  }
  
  
}
#plot(rf$predicted, outcome, col = rgb(red = abs(outcome - rf$predicted)/max(abs(outcome - rf$predicted)), green = 0, blue = 0))
write.table(results,"Diet_Regression_Classification_Results.tsv",sep="\t",quote=F,row.names=F)

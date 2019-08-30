install.packages("data.table")
install.packages("devtools")
devtools::install_github("moodymudskipper/cutr")
install.packages("corrplot")
install.packages("ggplot2")
library(ggplot2)
library(corrplot)
library(dplyr)
library(data.table)
library(cutr)


data_imputed <- fread("D:/Dossiers Maxime/Etudes/ENSAI/2ème année/Stage Novosibirsk/Project_LDScore/R/imputed/all_data.csv")
data_normRank <- fread("D:/Dossiers Maxime/Etudes/ENSAI/2ème année/Stage Novosibirsk/Project_LDScore/R/normRank/all_data.csv")
data_imputed$chi2 <- qchisq(data_imputed$PV, df=1, low = F)
data_normRank$chi2 <- qchisq(data_normRank$PV, df=1, low = F)

#Filtering
filter <- function(data){
  data <- data[data$MAF > 0.05]
  data <- data[data$iscores > 0.7]
  data <- data[!(data$Chr == 6 & data$pos < 35000000 & data$pos > 25000000)]
  data <- data[data$chi2 < 50]
  data <- data[data$L_aa_new < 1000]
  #data <- data[data$L_ad_new_adjust <= 1 & data$L_ad_new_adjust >= 0.1]
}
data_imputed <- filter(data_imputed)
data_normRank <- filter(data_normRank)

data_normRank$ratio <- data_normRank$Nsnp_used / data_normRank$Nsnp
print(median(data_normRank$ratio))
#---------Correlation matrix and Histograms---------

histo <- function(data){
  data2 <- data[data$L_aa_orig < 600]
  hist(data2$L_aa_orig, main = "histogram of L_aa_orig")
  
  data2 <- data[data$L_ad_orig < 20]
  hist(data2$L_ad_orig, main = "histogram of L_ad_orig")
  
  data2 <- data[data$L_aa < 600]
  hist(data2$L_aa, main = "histogram of L_aa")
  
  data2 <- data[data$L_ad < 20]
  hist(data2$L_ad, main = "histogram of L_ad")
  
  data2 <- data[data$L_aa_new < 600]
  hist(data2$L_aa_new, main = "histogram of L_aa_new")
  
  data2 <- data[data$L_ad_new < 20]
  hist(data2$L_ad_new, main = "histogram of L_ad_new")
  
  data2 <- data[data$L_ad_new_adjust < 20]
  hist(data2$L_ad_new_adjust, main = "histogram of L_ad_new_adjust")

}


correlation <- function(data){
  data_corr <- select(data, MAF, L_aa_orig, L_ad_orig, L_aa, L_ad, L_aa_new, L_ad_new,
                      L_ad_new_adjust, Nsnp, Nsnp_used)
  corrplot(cor(data_corr), p.mat = cor(data_corr), sig.level = -1, insig = "p-value", method = "square")
}

correlation(data_imputed)
histo(data_imputed)

correlation(data_normRank)
histo(data_normRank)




#------------Boxplot--------------
box <- function(data){
  data$chi2resid <- residuals(lm(data$chi2 ~ data$L_aa_new, na.action = na.exclude))
  
  data$L_aa_orig_categ = cut(data$L_aa_orig,
                               breaks = 20,
                               right = FALSE)
  means <- aggregate(data$chi2 ~ data$L_aa_orig_categ, FUN = mean)
  boxplot(data$chi2 ~ data$L_aa_orig_categ, outline = FALSE, main = "Boxplot L_aa_orig")
  points(means[2], col = "red")
  
  data$L_ad_orig_categ = cut(data$L_ad_orig,
                             breaks = 20,
                             right = FALSE)
  means <- aggregate(data$chi2resid ~ data$L_ad_orig_categ, FUN = mean)
  boxplot(data$chi2resid ~ data$L_ad_orig_categ, outline = FALSE, main = "Boxplot L_ad_orig")
  points(means[2], col = "red")
  
  data$L_aa_categ = cut(data$L_aa,
                             breaks = 20,
                             right = FALSE)
  means <- aggregate(data$chi2 ~ data$L_aa_categ, FUN = mean)
  boxplot(data$chi2 ~ data$L_aa_categ, outline = FALSE, main = "Boxplot L_aa")
  points(means[2], col = "red")
  
  data$L_ad_categ = cut(data$L_ad,
                        breaks = 20,
                        right = FALSE)
  means <- aggregate(data$chi2resid ~ data$L_ad_categ, FUN = mean)
  boxplot(data$chi2resid ~ data$L_ad_categ, outline = FALSE, main = "Boxplot L_ad")
  points(means[2], col = "red")
  
  data$L_aa_new_categ = cut(data$L_aa_new,
                            breaks = 20,
                            right = FALSE)
  means <- aggregate(data$chi2 ~ data$L_aa_new_categ, FUN = mean)
  boxplot(data$chi2 ~ data$L_aa_new_categ, outline = FALSE, main = "Boxplot L_aa_new")
  points(means[2], col = "red")
  
  data$L_ad_new_categ = cut(data$L_ad_new,
                               breaks = 20,
                               right = FALSE)
  means <- aggregate(data$chi2resid ~ data$L_ad_new_categ, FUN = mean)
  boxplot(data$chi2resid ~ data$L_ad_new_categ, outline = FALSE, main = "Boxplot L_ad_new")
  points(means[2], col = "red")
  
  data$L_ad_newAdj_categ = cut(data$L_ad_new_adjust,
                        breaks = 20,
                        right = FALSE)
  means <- aggregate(data$chi2resid ~ data$L_ad_newAdj_categ, FUN = mean)
  boxplot(data$chi2resid ~ data$L_ad_newAdj_categ, outline = FALSE, main = "Boxplot L_ad_new_adjust")
  points(means[2], col = "red")
  
}

box(data_normRank)

regression <- function(data){
  print("----------Model 1-------------")
  print(summary(lm(data$chi2 ~ data$L_aa_new, na.action = na.exclude)))
  
  print("----------Model 2-------------")
  print(summary(lm(data$chi2 ~ data$L_aa_new + data$L_ad_orig, na.action = na.exclude)))
  
  print("----------Model 3-------------")
  print(summary(lm(data$chi2 ~ data$L_aa_new + data$L_ad_new, na.action = na.exclude)))
  
  print("----------Model 4-------------")
  print(summary(lm(data$chi2 ~ data$L_aa_new + data$L_ad_new_adjust, na.action = na.exclude)))
  
  print("----------Model 5-------------")
  print(summary(lm(data$chi2 ~ data$L_ad_orig, na.action = na.exclude)))
  
  print("----------Model 6-------------")
  print(summary(lm(data$chi2 ~ data$L_ad_new, na.action = na.exclude)))
  
  print("----------Model 7-------------")
  print(summary(lm(data$chi2 ~ data$L_ad_new_adjust, na.action = na.exclude)))
  
  print("----------Model 8-------------")
  print(summary(lm(data$chi2 ~ data$L_aa_new + data$L_ad_new_adjust + data$L_aa_new*data$L_ad_new_adjust, na.action = na.exclude)))
  
  print("----------Model 9-------------")
  print(summary(lm(data$chi2 ~ data$MAF + data$L_aa_new + data$L_ad_new_adjust + data$L_aa_new*data$L_ad_new_adjust, na.action = na.exclude)))
}

regression(data_normRank)

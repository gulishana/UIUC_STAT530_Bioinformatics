# Problem 1

# (1) compare two different normal distributions
p <- NULL
for (i in 1:1000) {
  sample_A = rnorm(n = 5,mean = 10,sd = 3) 
  sample_B = rnorm(n = 5,mean = 30,sd = 9) 
  t_test_results <- t.test(sample_A,sample_B,paired = FALSE, alternative = "two.sided")
  p[i] <- t_test_results$p.value
}

# first rearrange p-values from the smallest to the largest
p = p[order(p)]

# FDR adjust p-values in p vecoter
m = length(p)
fdr_adjust_p = p.adjust(p,method = "fdr",n=m)
# list the first 30 values of p
head(fdr_adjust_p,n=30)


# then only store the smallest p-values that are less than 0.01 into a new vector x
x = p[which(p<0.01)]
# FDR adjust p-values in x using same total m
fdr_adjust_x = p.adjust(x,method = "fdr",n=m)
# list the first 30 values of FDR-adjusted x
head(fdr_adjust_x,n=30)


# list the last 30 values of FDR-adjusted x
N=length(fdr_adjust_x)
fdr_adjust_x[(N-29):N]
# compare those to the FDR adjusted p-values of the same orginal values from ordered p
fdr_adjust_p[(N-29):N]



# if we only look at the p-value distribution histogram of FDR-adjusted p within the range (0, 0.01)
fdr_adjust_p_below0.01 = fdr_adjust_p[which(fdr_adjust_p<0.01)]
length(fdr_adjust_p_below0.01)
hist(fdr_adjust_p_below0.01)

# if we only look at the p-value distribution histogram of FDR-adjusted x within the range (0, 0.01)
# so as to compare with previous histogram: fdr_adjust_p_below0.01
fdr_adjust_x_below0.01 = fdr_adjust_x[which(fdr_adjust_x<0.01)]
length(fdr_adjust_x_below0.01)
hist(fdr_adjust_x_below0.01)


# (2) compare two same normal distributions
p <- NULL
for (i in 1:1000) {
  sample_A = rnorm(n = 5,mean = 10,sd = 3) 
  sample_B = rnorm(n = 5,mean = 10,sd = 3) 
  t_test_results <- t.test(sample_A,sample_B,paired = FALSE, alternative = "two.sided")
  p[i] <- t_test_results$p.value
}

# first rearrange p-values from the smallest to the largest
p = p[order(p)]

# FDR adjust p-values in p vecoter
m = length(p)
fdr_adjust_p = p.adjust(p,method = "fdr",n=m)
# list the first 30 values of FDR-adjusted p
head(fdr_adjust_p,n=30)


# then only store the smallest p-values that are less than 0.1 into a new vector x
x = p[which(p<0.1)]
# FDR adjust p-values in x using same total m
fdr_adjust_x = p.adjust(x,method = "fdr",n=m)
# list the first 30 values of FDR-adjusted x
head(fdr_adjust_x,n=30)


# list the last 30 values of FDR-adjusted x
N=length(fdr_adjust_x)
fdr_adjust_x[(N-29):N]
# compare those to the FDR adjusted p-values of the same orginal values from ordered p
fdr_adjust_p[(N-29):N]
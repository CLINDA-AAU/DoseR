n <- 19
pval <- vector()
for(i in 1:100){
  a <- rpois(n, 4)
  b <- rpois(n, 2.4)
  
  data <- data.frame(n.bleeds = c(a,b), treatment = c(rep("pre", length(a)),  rep("treat", length(b))))
  ms <- summary(m1 <- glm(n.bleeds ~ treatment, family = "poisson", data=data))
  pval[i] <- ms$coefficients[2,4]
}

tab <- table(pval < 0.05)
tab[2]/sum(tab)*100

## ----eval=FALSE---------------------------------------------------------------
#  sunflowerplot(pressure)

## ----eval=FALSE---------------------------------------------------------------
#  n <- 1000
#  u <- runif(n)#设置随机数
#  x <- 2*(1-u)^(-(1/2)) # U=1-(2/X)^2，X=2*(1-U)^(1/2)
#  hist(x, prob = TRUE, main = expression(f(x)==8/x^3))#直方图显示
#  y <- seq(0, 50, 0.01)#密度函数曲线显示
#  lines(y, 8/(y^3))

## ----eval=FALSE---------------------------------------------------------------
#  n <- 1000 # 随机数个数  先通过beta随机数拟合，将结果与另一种定义方式的结果进行比较
#  y <- rbeta(n,2,2) #另Y=(X+1)/2，Y~Be(2,2),由Y产生随机数
#  x <- 2*y-1 #将Y产生的随机数结果回代
#  hist(x, prob = TRUE, main = expression(f(x)==(3/4)(1-x^2))) #由此得到直方图
#  z <- seq(-1, 1, 0.01) #将结果进行拟合
#  lines(z, (3/4)*(1-z^2))

## ----eval=FALSE---------------------------------------------------------------
#  nu <- 1000  #通过题目给出的公式定义得到相同的结果图像并进行拟合
#  i <- 0
#  for (i in 1:nu) {
#    u1 <- runif(nu,-1,1)
#    u2 <- runif(nu,-1,1)
#    u3 <- runif(nu,-1,1)
#    if(abs(u3[i]) >= abs(u2[i]) && abs(u3[i]) >= abs(u1[i])){
#      u[i] <- u2[i]}
#    else
#     { u[i] <- u3[i]}
#  }
#  hist(u, prob = TRUE, main = expression(f(x)==(3/4)(1-x^2)))
#  z <- seq(-1, 1, 0.01)
#  lines(z, (3/4)*(1-z^2))

## ----eval=FALSE---------------------------------------------------------------
#  n <- 1000
#  u <- runif(n)
#  x <- (2/(1-u)^(1/4))-2 #U=1-(2/(2+X))^4，X=(2/(1-U)^(1/4))-2
#  hist(x, prob = TRUE, main = expression(f(x)==64/(x+2)^5)) # 求导得密度函数
#  y <- seq(0, 20, 0.01) # 得到密度函数曲线
#  lines(y, 64/(y+2)^5)

## ----eval=FALSE---------------------------------------------------------------
#  m <- 1e4
#  x <- runif(m, min=0, max=pi/3)#设置10000个随机数，产生10000个从0到(pi/3)的均匀分布随机数
#  pihat <- mean(sin(x))*(pi/3)#(pi/3)E(sin(X))
#  print(c(pihat,-cos(pi/3) + cos(0)))#分别显示估计和真实的积分值

## ----eval=FALSE---------------------------------------------------------------
#  #先使用the simple Monte Carlo method
#  m <- 1e4
#  x <- runif(m, min=0, max=1)
#  theta.hat <- mean(exp(x))
#  #再使用antithetic variate approach
#  m_<- 1e4/2
#  x1 <- runif(m_, min=0, max=1)
#  x2 <- runif(m_, min=0, max=1)
#  theta_.hat <- (mean(exp(x1))+mean(exp(1-x2)))/2
#  print(c(theta.hat,theta_.hat,exp(1)-1))

## ----eval=FALSE---------------------------------------------------------------
#  m <- 1e4
#  f1 <- function(u)u^(1/4)#g(x)/f_1(x)=1/x^4,g(x)/f_2(x)=exp(-(1/2*x^2))，将其换
#  f2 <- function(u)exp(-(1/2*u^2))# 成（0，1）区间的值更便于计算
#  set.seed(510);  u <- runif(m)
#  T1 <- f1(u); T2 <- f2(u)
#  c(mean(T1), mean(T2))
#  c(sd(T1)/sqrt(m), sd(T2)/sqrt(m))#由结果可知，第二种估计的方差较小，但是不一定是良好的估计，因为你可选择的范围比较多，所以要对结果进行检验

## ----eval=FALSE---------------------------------------------------------------
#  M <- 10000
#  k <- 5 #分成五个区间
#  #由上面的估计结果和接下来的估计结果可以知道，虽然上面的方法可适用的范围比较广，但是也很容易存在较大的偏差。相较起来，本方法的准确性更高，且具有更小的方差，所以本估计方法是一种优良的估计方法
#  r <- M/k
#  N <- 50
#  T2 <- numeric(k)
#  est <- matrix(0, N, 2)
#  g<-function(t)(1-exp(-1))/(1+t^2)*(t>0)*(t<1)#g(x)/f(x)=(1-exp(-1))/(1+x^2)
#  for (i in 1:N) {
#    est[i, 1] <- mean(g(runif(M)))
#    for(j in 1:k)T2[j]<-mean(g(runif(M/k,(j-1)/k,j/k)))
#    est[i, 2] <- mean(T2)
#  }
#  apply(est,2,mean)#分层重要抽样法的适用范围更广
#  apply(est,2,sd)#相对于例5.10的结果来说，分层重要抽样法得到的结果更加精确，并且标准差更小
#  

## ----eval=FALSE---------------------------------------------------------------
#  n <- 20#X服从对数正态分布，Y=ln(x)~N(μ,σ^2)，所以可以直接利用Y进行估计，再代入x即可
#  alpha <- .05#置信度95%的置信区间求得的结果为(mean(y)-sqrt(var(y))* qt(1-alpha/2,df = n-1) / sqrt(n),mean(y)+sqrt(var(y))* qt(1-alpha/2,df = n-1) / sqrt(n))
#  m <- 1000
#  mu <- 0
#  s<-0
#  UCL <- USL <- numeric(m)
#  est <- matrix(0, m, 2)
#  for (i in 1:m) {
#    y <- rnorm(n,mean = mu,sd =2)#y是正态分布随机量
#    est[i, 1] <- mean(y)
#    est[i, 2] <- sqrt(var(y))* qt(1-alpha/2,df = n-1) / sqrt(n)
#    UCL[i] <- est[i, 1] + est[i, 2]
#    USL[i] <- est[i, 1] - est[i, 2]
#    s<-s+(USL[i] < mu && mu < UCL[i] )#统计区间包含真实值的数目
#  }
#  
#  s/m#统计区间包含真实值的比例

## ----eval=FALSE---------------------------------------------------------------
#  n <- 20
#  alpha <- .05
#  m <- 1000
#  mu <- 2
#  s <- 0
#  UCL <- USL <- numeric(m)
#  est <- matrix(0, m, 2)
#  for (i in 1:m) {
#    y <- rchisq(n,df = 2)#令样本改成卡方分布
#    est[i, 1] <- mean(y)
#    est[i, 2] <- sqrt(var(y))* qt(1-alpha/2,df = n-1) / sqrt(n)
#    UCL[i] <- est[i, 1] + est[i, 2]
#    USL[i] <- est[i, 1] - est[i, 2]
#    s<-s+(USL[i] < mu && mu < UCL[i] )#样本抽样方法改变后，得到的数据相对不稳定，置信区间的估计并不良好
#  }
#  s/m

## ----eval=FALSE---------------------------------------------------------------
#  sk <- function(x){
#    #计算样本偏度系数
#    xbar <- mean(x)
#    m3 <- mean((x-xbar)^3)
#    m2 <- mean((x-xbar)^2)
#    return( m3 / m2^1.5)
#  }
#  
#  alpha <- .05
#  n <- 30
#  m <- 2500
#  #epsilon <- c(seq(0.01,.15,.01),seq(.15,1,.05))
#  epsilon1 <- c(seq(0.01,.15,.01),seq(.15,1,.05),seq(1,20,1))#epsilon表示α参数
#  N <- length(epsilon1)
#  pwr <- numeric(N)
#  cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
#  
#  for(j in 1:N){
#    e <- epsilon1[j]
#    sktests <- numeric(m)
#    for (i in 1:m) {
#      x <- rbeta(n,e,e)
#      sktests[i] <- as.integer(abs(sk(x)) >= cv)
#    }
#    pwr[j] <- mean(sktests)
#  }
#  
#  
#  plot(epsilon1,pwr,type = "b",xlab = bquote(eplison1),ylim = c(0,0.1))
#  abline(h = .1,lty = 3)
#  se <- sqrt(pwr * (1-pwr) / m)
#  lines(epsilon1,pwr+se,lty = 3)
#  lines(epsilon1,pwr-se,lty = 3)
#  
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  sk <- function(x){
#    #计算样本偏度系数
#    xbar <- mean(x)
#    m3 <- mean((x-xbar)^3)
#    m2 <- mean((x-xbar)^2)
#    return( m3 / m2^1.5)
#  }
#  
#  alpha <- .05#以t(v)为例
#  n <- 30
#  m <- 2500
#  epsilon2 <- c(seq(1,20,1))
#  N <- length(epsilon2)
#  pwr <- numeric(N)
#  cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
#  
#  for(j in 1:N){
#    e <- epsilon2[j]
#    sktests <- numeric(m)
#    for (i in 1:m) {
#      x <- rt(n,e)
#      sktests[i] <- as.integer(abs(sk(x)) >= cv)
#    }
#    pwr[j] <- mean(sktests)
#  }
#  
#  plot(epsilon2,pwr,type = "b",xlab = bquote(eplison2),ylim = c(0,1))
#  abline(h = 1,lty = 3)
#  se <- sqrt(pwr * (1-pwr) / m)
#  lines(epsilon2,pwr+se,lty = 3)
#  lines(epsilon2,pwr-se,lty = 3)
#  

## ----eval=FALSE---------------------------------------------------------------
#  count5test <- function(x,y){
#    X <- x - mean(x)
#    Y <- y - mean(y)
#    outx <- sum(X > max(Y)) + sum(X < min(Y))
#    outy <- sum(Y > max(X)) + sum(Y < min(X))
#    return(as.integer(max(c(outx,outy)) > 5))
#  }
#  n <- c(20,200,1000)#分别对应小样本、中样本和大样本
#  mu1 <- mu2 <- 0
#  sigma1 <- 1
#  sigma2 <- 1.5
#  m <- 10000
#  power1 <- power2 <- numeric(length(n))
#  set.seed(1234)
#  for(i in 1:length(n)){
#    power1[i] <- mean(replicate(m,expr = {
#    x <- rnorm(n[i],mu1,sigma1)
#    y <- rnorm(n[i],mu2,sigma2)
#    x <- x - mean(x)
#    y <- y - mean(y)
#    count5test(x,y)
#    }))
#     pvalues <- replicate(m,expr={
#      x <- rnorm(n[i],mu1,sigma1)
#      y <- rnorm(n[i],mu2,sigma2)
#      Ftest <- var.test(x, y, ratio = 1,
#                        alternative = c("two.sided", "less", "greater"),
#                        conf.level = 0.945, ...)
#      Ftest$p.value})
#      power2[i] <- mean(pvalues<=0.055)
#  }
#  
#  power1
#  power2

## ----eval=FALSE---------------------------------------------------------------
#  library(MASS)# pressure2, echo=FALSE
#  alpha <- .05
#  n <- c(10,20,30,50,100,500)
#  cv <- qnorm(1-alpha/2,0,sqrt(6/n))
#  
#  
#  skn <- function(x) {#多维下的峰度
#    xbar <- mean(x)
#    n <- nrow(x)
#    skn <- 1/n^2*sum(((x-xbar)%*%solve(cov(x))%*%t(x-xbar))^3)
#    return(skn)
#  }
#  
#  p.reject <- numeric(length(n))
#  m <- 10000
#  s <- matrix(c(1,0,0,1),2,2)
#  for (i in 1:length(n)) {
#    sktests <- numeric(m)
#    for (j in 1:m) {
#      x <- mvrnorm(n[i],c(0,0),s)
#      sktests[j] <- as.integer(abs(skn(x)) >= cv[i])
#    }
#    p.reject[i] <- mean(sktests)
#  }
#  p.reject
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(MASS)
#  n <- c(10,20,30,50,100,500)
#  
#  sk <- function(x){
#    #计算样本偏度系数
#    xbar <- mean(x)
#    m3 <- mean((x-xbar)^3)
#    m2 <- mean((x-xbar)^2)
#    return( m3 / m2^1.5)
#  }
#  
#  
#  skn <- function(x) {
#    xbar <- mean(x)
#    n <- nrow(x)
#    skn <- 1/n^2*sum(((x-xbar)%*%solve(cov(x))%*%t(x-xbar))^3)
#    return(skn)
#  }
#  
#  alpha <- .1
#  n <- 30
#  m <- 2500
#  #s <- matrix(c(1,0,0,1),2,2)
#  epsilon <- c(seq(0,.15,.01),seq(.15,1,.05))
#  N <- length(epsilon)
#  pwr <- numeric(N)
#  cv <- qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
#  
#  for(j in 1:N){
#    e <- epsilon[j]
#    sktests <- numeric(m)
#    for (i in 1:m) {
#     sigma <- sample (c(1,10),replace = TRUE,size = n,prob = c(1-e,e))
#     #s <- matrix(c(1,0,0,1),2,2)
#     x <- rnorm(n,0,sigma)
#     sktests[i] <- as.integer(abs(sk(x)) >= cv)
#    }
#    pwr[j] <- mean(sktests)
#  }
#  
#  
#  plot(epsilon,pwr,type = "b",xlab = bquote(eplison),ylim = c(0,1))
#  abline(h = .1,lty = 3)
#  se <- sqrt(pwr * (1-pwr) / m)
#  lines(epsilon,pwr+se,lty = 3)
#  lines(epsilon,pwr-se,lty = 3)
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(bootstrap)
#  b.cor <- function(x,i) cor(x[i,1],x[i,2])
#  cor.hat <- cor(law$LSAT, law$GPA)
#  n <- nrow(law)
#  bias.jack <- se.jack <- cor.jack <- numeric(n)
#  for (i in 1:n) {
#    cor.jack[i] <- cor(law[(1:n)[-i],1],law[(1:n)[-i],2])
#  }
#  bias.jack <- (n-1)*(mean(cor.jack)-cor.hat)
#  se.jack <- sqrt((n-1)*mean((cor.jack-cor.hat)^2))
#  round(c(original=cor.hat,bias.jack=bias.jack, se.jack=se.jack),3)
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(boot)
#  x <- c(3,5,7,18,43,85,91,98,100,130,230,487)
#  lambdafe<-mean(x)
#  boot.mean <- function(x,i) mean(x[i])
#   # sz <- rexp(n,lambda)
#    de <- boot(data=x,statistic=boot.mean, R = 999)
#    ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
#    boot.ci (de)#

## ----eval=FALSE---------------------------------------------------------------
#  library(bootstrap)
#  x <- cov(scor,scor)
#  b.cor <- function(x,i) cor(x[i,1],x[i,2])
#  lamba <- eigen(x)$values
#  sumlam <- sum(lamba)
#  theta <- lamba[1] / sumlam
#  y <- numeric(length(lamba))
#  B <- 1e4
#  thetastar <- numeric(B)
#  for (i in 1:length(lamba)){
#    y[i] <- lamba[i] / sumlam
#  }
#  
#  for(b in 1:B){
#    xstar <- sample(y,replace=TRUE)
#    thetastar[b] <- mean(xstar)
#    }
#  round(c(Bias.boot=mean(thetastar)-theta, SE.boot=sd(thetastar)),3)

## ----eval=FALSE---------------------------------------------------------------
#  library(DAAG)
#  attach(ironslag)
#  n <- length(magnetic)
#  e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)
#  for (i in 1:n) {
#  
#    for (j in (i+1):n-1){
#      y <- magnetic[-c(i,j)]
#      x <- chemical[-c(i,j)]
#  
#      J1 <- lm(y ~ x)
#      yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i,j)]
#      e1[i,j] <- mean((magnetic[c(i,j)] - yhat1)^2)
#  
#      J2 <- lm(y ~ x + I(x^2))
#      yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(i,j)] + J2$coef[3] *
#        chemical[c(i,j)]^2
#      e2[i,j] <- mean((magnetic[c(i,j)] - yhat2)^2)
#  
#      J3 <- lm(log(y) ~ x)
#      logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
#      yhat3 <- exp(logyhat3)
#      e3[i,j] <-mean(( magnetic[c(i,j)] - yhat3)^2)
#  
#      J4 <- lm(log(y) ~ log(x))
#      logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
#      yhat4 <- exp(logyhat4)
#      e4[i,j] <- mean((magnetic[c(i,j)] - yhat4)^2)
#  
#    }
#  
#  }
#  c(2*sum(e1)/(n*(n-1)),2*sum(e2)/(n*(n-1)),2*sum(e3)/(n*(n-1)),2*sum(e4)/(n*(n-1)))
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(111)#
#  count5test <- function(x, y) {
#  
#  X <- x - mean(x)
#  Y <- y - mean(y)
#  weight <- round(10 * (length(x)/(length(x)+length(y))))
#  outx <- sum(X > max(Y)) + sum(X < min(Y))
#  outy <- sum(Y > max(X)) + sum(Y < min(X))
#  # return 1 (reject) or 0 (do not reject H0)
#  return(as.integer(min(c(outx, outy)) > weight || max(c(outx, outy)) > (10-weight) ))
#  }
#  
#  R <- 10000
#  n1 <- 10
#  n2 <- 20
#  weight <- round(10 * (n1/(n1+n2)))
#  mu1 <- mu2 <- 0
#  sigma1 <- sigma2 <- 1
#  x <- rnorm(n1, mu1, sigma1)
#  y <- rnorm(n2, mu2, sigma2)
#  
#  z <- c(x,y)
#  K <- 1:(n1 + n2)
#  n <- length(x)
#  jg <- numeric(R)
#  for (i in 1:R) {
#    k <- sample(K, size = n, replace = FALSE)
#    x1 <- z[k];y1 <- z[-k]
#    x1 <- x1 - mean(x1) #centered by sample mean
#    y1 <- y1 - mean(y1)
#   jg[i] <- count5test(x1, y1)
#  }
#  
#  q <- mean(jg)
#  round(c(q),3)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12)
#  library(RANN)#Unequal variances and equal expectations
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # what's the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  m <- 1e2; k<-3; p<-2; mu <- 0.3
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=999,
#    sim = "permutation", sizes = N,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  
#  p.values <- matrix(NA,m,3)
#  
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1),ncol=p);
#    y <- cbind(rnorm(n2,0,1.2),rnorm(n2,0,1.2));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#    p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  p.value1 <- mean(p.values[i,1])
#  p.value2 <- mean(p.values[i,2])
#  p.value3 <- mean(p.values[i,3])
#  round(c(p.value1,p.value2,p.value3),3)#nearest NN test and Ball test are generally less powerful than Energy test

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123)
#  library(RANN)#Unequal variances and unequal expectations
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # what's the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  m <- 1e2; k<-3; p<-2; mu <- 0.3
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=999,
#    sim = "permutation", sizes = N,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  
#  p.values <- matrix(NA,m,3)
#  
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1),ncol=p);
#    y <- cbind(rnorm(n2,0.1,1.2),rnorm(n2,0.1,1.2));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#    p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  p.value1 <- mean(p.values[i,1])
#  p.value2 <- mean(p.values[i,2])
#  p.value3 <- mean(p.values[i,3])
#  round(c(p.value1,p.value2,p.value3),3)#Nearest NN test and Ball test are generally less powerful than Energy test
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(125)
#  library(RANN)#Non-normal distributions: t distribution with 1 df (heavy-tailed distribution),
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # what's the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  m <- 1e2; k<-3; p<-2; mu <- 0.3
#  n1 <- n2 <- 10; R<-999; n <- n1+n2; N = c(n1,n2)
#  
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=999,
#    sim = "permutation", sizes = N,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  
#  p.values <- matrix(NA,m,3)
#  
#  for(i in 1:m){
#    x <- matrix(rt(n1*p,n1),ncol=p);
#    y <- cbind(rt(n2,n2),rt(n2,n2));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#    p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  p.value1 <- mean(p.values[i,1])
#  p.value2 <- mean(p.values[i,2])
#  p.value3 <- mean(p.values[i,3])
#  round(c(p.value1,p.value2,p.value3),5)#Energy test and nearest NN test are generally less powerful than  Ball test
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12345)
#  library(RANN)# Non-normal distributions: bimodel distribution (mixture of two normal distributions)
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # what's the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  m <- 1e2; k<-3; p<-2; mu <- 0.1
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=999,
#    sim = "permutation", sizes = N,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1),ncol=p);
#    y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#    p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  p.value1 <- mean(p.values[i,1])
#  p.value2 <- mean(p.values[i,2])
#  p.value3 <- mean(p.values[i,3])
#  round(c(p.value1,p.value2,p.value3),3)#Energy test and Ball test are generally less powerful than nearest NN test

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123456)
#  library(RANN)# Unbalanced samples (say, 1 case versus 10 controls)
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # what's the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  m <- 1e2; k<-3; p<-2; mu <- 0.1
#  n1 <- 5
#  n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=999,
#    sim = "permutation", sizes = N,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1),ncol=p);
#    y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#    p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  p.value1 <- mean(p.values[i,1])
#  p.value2 <- mean(p.values[i,2])
#  p.value3 <- mean(p.values[i,3])
#  round(c(p.value1,p.value2,p.value3),3)# nearest NN test and Ball test are generally less powerful than Energy test

## ----eval=FALSE---------------------------------------------------------------
#  rw.Metropolis <- function(n, sigma, x0, N) {
#  x <- numeric(N)
#  x[1] <- x0
#  mu <- 1
#  b <- 1
#  u <- runif(N,-.5,.5)
#  R <- mu-b*sign(u)*log(1-2*abs(u))
#  k <- 0
#  for (i in 2:N) {
#  y <- rnorm(1, x[i-1], sigma)
#  if (R[i] <= (dt(y, n) / dt(x[i-1], n)))
#  x[i] <- y else {
#  x[i] <- x[i-1]
#  k <- k + 1
#  } }
#  return(list(x=x, k=k))
#  }
#  n <- 4 #degrees of freedom for target Student t dist.
#  N <- 2000
#  sigma <- c(.05, .5, 2, 5)
#  x0 <- 25
#  rw1 <- rw.Metropolis(n, sigma[1], x0, N)
#  rw2 <- rw.Metropolis(n, sigma[2], x0, N)
#  rw3 <- rw.Metropolis(n, sigma[3], x0, N)
#  rw4 <- rw.Metropolis(n, sigma[4], x0, N)
#  #number of candidate points rejected
#   print(c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N))
#  

## ----eval=FALSE---------------------------------------------------------------
#  Gelman.Rubin <- function(psi) {
#    # psi[i,j] is the statistic psi(X[i,1:j])
#    # for chain in i-th row of X
#    psi <- as.matrix(psi)
#    n <- ncol(psi)
#    k <- nrow(psi)
#    psi.means <- rowMeans(psi) #row means
#    B <- n * var(psi.means) #between variance est.
#    psi.w <- apply(psi, 1, "var") #within variances
#    W <- mean(psi.w) #within est.
#    v.hat <- W*(n-1)/n + (B/n) #upper variance est.
#    r.hat <- v.hat / W #G-R statistic
#    return(r.hat)
#  }
#  normal.chain <- function(sigma, N, X1) {
#    #generates a Metropolis chain for Normal(0,1)
#    #with Normal(X[t], sigma) proposal distribution
#    #and starting value X1
#    x <- rep(0, N)
#    x[1] <- X1
#    #u <- runif(N)
#    u0 <- runif(N,-.5,.5)
#    u <- 1-1*sign(u0)*log(1-2*abs(u0))
#    for (i in 2:N) {
#      xt <- x[i-1]
#      y <- rnorm(1, xt, sigma) #candidate point
#      r1 <- dnorm(y, 0, 1) * dnorm(xt, y, sigma)
#      r2 <- dnorm(xt, 0, 1) * dnorm(y, xt, sigma)
#      r <- r1 / r2
#      if (u[i] <= r) x[i] <- y else
#        x[i] <- xt
#    }
#    return(x)
#  }
#  
#  
#  
#  sigma <- .05 #parameter of proposal distribution
#  #sigma <- c(.05, .5, 2, 4)
#  k <- 4 #number of chains to generate
#  n <- 2000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- 25
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#    X[i, ] <- normal.chain(sigma, n, x0)
#  #compute diagnostic statistics
#    psi <- t(apply(X, 1, cumsum))
#      for (i in 1:nrow(psi))
#        psi[i,] <- psi[i,] / (1:ncol(psi))
#    print(Gelman.Rubin(psi))
#  
#  sigma <- .5 #parameter of proposal distribution
#  #sigma <- c(.05, .5, 2, 16)
#  k <- 4 #number of chains to generate
#  n <- 2000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- 25
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#    X[i, ] <- normal.chain(sigma, n, x0)
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#    psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  sigma <- 2 #parameter of proposal distribution
#  #sigma <- c(.05, .5, 2, 16)
#  k <- 4 #number of chains to generate
#  n <- 2000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- 25
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#    X[i, ] <- normal.chain(sigma, n, x0)
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#    psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  
#  sigma <- 5 #parameter of proposal distribution
#  k <- 4 #number of chains to generate
#  n <- 2000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- 25
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#    X[i, ] <- normal.chain(sigma, n, x0)
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#    psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  #Error in if (u[i] <= r) x[i] <- y else x[i] <- xt : Where TRUE/FALSE values are required, missing values cannot be used，In terms of integers, sigma counts up to 5
#  #ˆR < 1.2.

## ----eval=FALSE---------------------------------------------------------------
#  m <-100
#  int <- numeric(m)
#  findpoint <- function(k,a){
#    for (i in 1:m) {
#      tk <- rt(1,k)
#      u <- sqrt(a^2*(k)/(k+1-a^2))
#      int[i] <- as.numeric(tk > u)
#    }
#    return (mean(int))
#  }
#  findpoint2 <- function(k){
#  
#    jg <- numeric(25)
#    f1 <- f2 <- numeric(sqrt(k)*10)
#    j <- 0
#    i <-1
#    while (j < sqrt(k)) {
#      f1[i] <- findpoint(k,j)
#      f2[i] <- findpoint(k-1,j)
#      i <- i + 1
#      j <- j + 0.1
#    }
#    m1 <- lowess(f1, y=NULL, f = 2/3, iter = 3)
#    m2 <- lowess(f2, y=NULL, f = 2/3, iter = 3)
#    for (k in 1:length(f1)) {
#      if (isTRUE(all.equal(f1[k],f2[k])))
#      {
#        jg <- k*0.1
#        break
#        #only find one point(minpoint)
#      }
#  
#    }
#    return (jg)
#  }
#  r <- c(4:25,100,500,1000)
#  jg2 <- numeric(25)
#  for (i in 1:25) {
#    k <- r[i]
#    jg2[i] <- findpoint2(k)
#  }
#  
#  print(c(jg2))
#  #no matter
#  #for (n in 1:25) {
#  
#   # for (j in 1:sqrt(k)) {
#     # f1[j] <- findpoint(k,j)
#    #  f2[j] <- findpoint(k-1,j)
#  #}
#  
#  
#    #f <- findpoint(k-1,a)
#  #plot(f1)
#  #m1 <- lowess(f1, y=NULL, f = 2/3, iter = 3)
#  #lines(lowess(m1$x,m1$y))
#  #plot(m1$x,m1$y)
#  
#  #plot(f2)
#  #m2 <- lowess(f2, y=NULL, f = 2/3, iter = 3)
#  #lines(lowess(m2$x,m2$y))
#  #plot(m2$x,m2$y)
#  #}
#  
#  #jg

## ----eval=FALSE---------------------------------------------------------------
#  na<-444
#  nb<-132
#  noo<-361
#  nab<-63
#  n<-na+nb+noo+nab
#  r0<-sqrt(noo/n)#
#  p0<-0.5
#  q0<-0.5
#  
#  kp <- (nab + 2*na*(1-q0)/(2-p0-2*q0))/(2*n)
#  kq <- (nab + 2*nb*(1-p0)/(2-2*p0-q0))/(2*n)
#  ko <- 1-kp-kq
#  m<-100
#  p<-q<-r<-numeric(m)
#  q[1]<-kq
#  p[1]<-kp
#  r[1]<-r0
#  
#  like_lihood <- numeric(m)
#  for (i in 2:m) {
#    kp<-(nab + 2*na*(1-q[i-1])/(2-p[i-1]-2*q[i-1]))/(2*n)
#    kq<-(nab + 2*nb*(1-p[i-1])/(2-2*p[i-1]-q[i-1]))/(2*n)
#    q[i]<-kq
#    p[i]<-kp
#    r[i]<-r0
#    like_lihood[i] <- (2*p[i]*q[i])^nab+r[i]^(2*noo)+(2*p[i]*r[i])^na+(2*q[i]*r[i])^nb+(p[i]/(2*r[i]))^((na*p[i])/(p[i]+2*r[i]))+(q[i]/(2*r[i]))^((nb*q[i])/(q[i]+2*r[i]))
#  }
#  
#  plot(like_lihood)#they are increasing
#  c(p[m],q[m])

## ----eval=FALSE---------------------------------------------------------------
#  x <- mtcars
#  formu <- list(
#  mtcars$mpg ~ mtcars$disp,
#  mtcars$mpg ~ I(1 / mtcars$disp),
#  mtcars$mpg ~ mtcars$disp + mtcars$wt,
#  mtcars$mpg ~ I(1 / mtcars$disp) + mtcars$wt
#  )
#  l <- length(formu)
#  mt <- numeric(l)
#  for (i in 1:l) {
#    mt[i] <- lapply(x,function(x) lm(formu[[i]]))#y = a + bx
#  }
#  
#  y <- x_true <- matrix(0,4,length(mtcars$disp))#
#  x_true[1,] <- mtcars$disp
#  x_true[2,] <- I(1 / mtcars$disp)
#  x_true[3,] <- mtcars$disp + mtcars$wt
#  x_true[4,] <- I(1 / mtcars$disp) + mtcars$wt
#  for (i in 1:l) {
#    a <- mt[[i]][["coefficients"]][1]
#    b <- mt[[i]][["coefficients"]][2]
#    y[i,] <- a + b * x_true[i,]
#    plot(x_true[i,], y[i,])
#    lines(x_true[i,],y[i,])
#  }
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123)
#  trials <- replicate(100, t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
#  #x <- trials[[i]]$p.value
#  sapply(trials,function(x) list(x$p.value))
#  

## ----eval=FALSE---------------------------------------------------------------
#  #Extra challenge
#  #trials <- replicate(100, t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
#  #x <- trials[[i]]$p.value
#  #for (i in 1:length(trials)) {
#  #i <- 1
#  #  (function(x) print(as.vector(x[[i]]$p.value)))(trials)
#  #}
#   sapply(trials,"[[","p.value")

## ----eval=FALSE---------------------------------------------------------------
#  testlist <- list(iris, mtcars)
#  #testlist <- list(rnorm(1000), runif(1000), rpois(1000,100))
#  #system.time(lmapply(testlist, mean, numeric(1)))
#  #system.time(vapply(testlist, mean, numeric(1)))
#  #system.time(Map(mean, testlist))
#  lmapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
#  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
#  if(simplify == TRUE){return(simplify2array(out))}
#  out <- as.array.default(out)
#  out
#  }
#  lmapply(testlist, mean, numeric(1))
#  m <- lmapply(testlist, mean, numeric(1))
#  class(m)
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(Rcpp)
#  dir_cpp <- '../Rcpp/'
#  # Can create source file in Rstudio
#  sourceCpp(paste0(dir_cpp,"rw_MetropolisC.cpp"))
#  library(microbenchmark)
#  
#  rw.Metropolis <- function(n, sigma, x0, N) {
#  x <- numeric(N)
#  x[1] <- x0
#  mu <- 1
#  b <- 1
#  u <- runif(N,-.5,.5)
#  R <- mu-b*sign(u)*log(1-2*abs(u))
#  k <- 0
#  for (i in 2:N) {
#  y <- rnorm(1, x[i-1], sigma)
#  if (R[i] <= (dt(y, n) / dt(x[i-1], n)))
#  x[i] <- y else {
#  x[i] <- x[i-1]
#  k <- k + 1
#  } }
#  return(list(x=x, k=k))
#  }
#  #
#  #Random number comparison (acceptance rate)
#  #
#  n <- 4 #degrees of freedom for target Student t dist.
#  N <- 2000
#  sigma <- c(.05, .5, 2, 5)
#  x0 <- 25
#  
#  rw1 <- rw.Metropolis(n, sigma[1], x0, N)
#  rw2 <- rw.Metropolis(n, sigma[2], x0, N)
#  rw3 <- rw.Metropolis(n, sigma[3], x0, N)
#  rw4 <- rw.Metropolis(n, sigma[4], x0, N)#number of candidate points rejected
#   print(c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N))
#  rw1C <- rw_MetropolisC(sigma[1], x0, N)
#  rw2C <- rw_MetropolisC(sigma[2], x0, N)
#  rw3C <- rw_MetropolisC(sigma[3], x0, N)
#  rw4C <- rw_MetropolisC(sigma[4], x0, N)
#  print(c(rw1C[[2]]/N, rw2C[[2]]/N, rw3C[[2]]/N, rw4C[[2]]/N))
#  #
#  #In terms of acceptance rate, the latter has a lower acceptance rate
#  #
#  ts <- microbenchmark(rw.MetropolisR=rw.Metropolis(n, sigma[1], x0, N),meanR2=rw_MetropolisC(x0,sigma[1],N))#
#  summary(ts)[,c(1,3,5,6)]#
#  #
#  #In terms of time, the latter takes a lot less time
#  #

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
#  List rw_MetropolisC(double x0, double sigma, int N) {
#    NumericVector x(N);
#    as<DoubleVector>(x)[0] = x0;
#    NumericVector u(N);
#    u = as<DoubleVector>(runif(N));
#    List out(2);
#    int k = 1;
#    for(int i=1;i<(N-1);i++){
#      double y = as<double>(rnorm(1,x[i-1],sigma));
#      if(u[i] <= exp(abs(x[i-1])-abs(y))){
#        x[i] = y;
#        k+=1;
#      }
#      else{
#        x[i] = x[i-1];
#      }
#    }
#    out[0] = x;
#    out[1] = k;
#    return(out);
#  }


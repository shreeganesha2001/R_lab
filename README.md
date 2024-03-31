## prac1
1.
X=c()
T={}
n=length(X)
x1=seq(1/n,1,by=1/n)
y=sort(X)
for(i in 1:n)
 {
   T[i]=(sum(y[1:i])+(n-i)*y[i])/sum(y[1:n])
   }
plot(x1,T,ylim=c(0,1))
lines(x1,x1,type="l")
 
#gamma distribution
 logL1=function(a1)
     {
      b1=a1/mean(X)
       logL1=(n*a1*log(b1))-(n*log(gamma(a1)))-(b1*sum(X))+((a1-1)*sum(log(X)))
       return (logL1)
       }
 opt1=optimize(logL1,c(0,30),maximum=TRUE)
 a1=opt1$maximum ; a1
 b1=a1/mean(X) ; b1
 opt1
 
 #weibull distribution
logL2=function(a2)
   {
     b2=n/sum(X^a2)
     logL2=(n*log(a2))+(n*log(b2))+((a2-1)*sum(log(X)))-n
      return(logL2)
     }
opt2=optimize(logL2,c(0,30),maximum=TRUE)
a2=opt2$maximum ; a2
b2=n/sum(X^a2) ; b2
opt2 

 #exponentiated exponential
logL3=function(b3)
 {
   a3=-n/sum(log(1-exp(-b3*X)))
   logL3=(n*log(a3))+(n*log(b3))+((a3-1)*sum(log(1-exp(-b3*X))))-(b3*sum(X))
   return(logL3)
  }
 opt3=optimize(logL3,c(0,30),maximum=TRUE)
 b3=opt3$maximum ; b3
 a3=-n/sum(log(1-exp(-b3*X))); a3
opt3

hist(y,prob=TRUE,ylim=c(0,1)) 
lines(y,dgamma(y,a1,b1),type="l")
lines(y,dweibull(y,a2,1/b2^(1/a2)),type="l",col=2)
fun3=a3*b3*(1-exp(-b3*y))^(a3-1)*exp(-b3*y)
lines(y,fun3,type="l",col=3)
A=c("Gamma","Weibull","EE")
legend(x="topleft",legend=c(A),col=c(1,2,3),lty=1)


p=1-pweibull(3,a2,1/b2^(1/a2))
p

 #survival function
plot(y,(1-pweibull(y,a2,1/b2^(1/a2))),xlab="failure stress",ylab="survival rate",type="l")
 
#hazard function
plot(y,dweibull(y,a2,1/b2^(1/a2))/(1-pweibull(y,a2,1/b2^(1/a2))),ylab="failure rate",xlab="failure stress",type="l")
  
2.
x = c()
y = sort(x)
n = length(y)
z =c()
t =c()

for(i in 1:n)
{
  
  z[i] = i/n
  t[i] = (sum(y[1:i])+(n-i)*(y[i]))/(sum(y))
  
}

plot(z,t,ylim=c(0,1),main="TTT plot",xlab = "i/n",ylab="T(i/n)")
lines(z,z,type="l")

 
gam = fitdistr(x,"gamma")
gam
     shape         rate   
  0.92159083   0.00989460 
 (0.07768517) (0.00108513)

ga_shape = 0.92159083 
ga_rate =  0.00989460 
gam$loglik

[1] -1178.291


########## Weibull #########

logL1 = function(alpha)
{
  
  beta =  n / sum(x^alpha)
  logL1 =  n*(log(alpha*beta))+(alpha-1)*sum(log(x))-beta*sum(x^alpha)
  return(logL1)
}

opt1 = optimize(logL1,c(0,10),maximum=TRUE)
alpha1=opt1$maximum;alpha1

[1] 0.9245503


beta1 =  n / sum(x^alpha1);beta1

[1] 0.01567414

opt1

$maximum
[1] 0.9245503

$objective
[1] -1177.585


#####  Exponentiated Exponential #######


logL2 = function(beta2)
{
  
  alpha2 = n/(-1*sum(log(1-exp(-beta2*x))))
  logL2 = n*log(alpha2*beta2)+(alpha2-1)*sum(log(1-exp(-beta2*x)))-beta2*sum(x)
  return(logL2)
}

opt2 = optimize(logL2,c(0,10),maximum=TRUE);opt2

$maximum
[1] 0.01020508

$objective
[1] -1178.402

beta2 = opt2$maximum;beta2
[1] 0.01020508
alpha2 = n/(-1*sum(log(1-exp(-beta2*x))));alpha2
[1] 0.9266818

####    PDF    ####

fee = function(x)
{
  fee = alpha2*beta2*(1-exp(-beta2*x))^(alpha2-1)*exp(-beta2*x)
  return(fee)
}

#### CDF ####

Fee = function(x)
{
  Fee = (1-exp(-beta2*x))^alpha2
  return(Fee)
}


########### Rayleigh ###########

logL3 = function(lam)
{
  alpha3 = n/(-1*sum(log(1-exp(-(lam*x)^2))))
  logL3 = n*log(alpha3*(lam^2)*2)+(sum(log(x)))+(alpha3-1)*sum(log(1-exp(-(lam*x)^2)))-sum((lam*x)^2)
  return(logL3)
}

opt3 = optimize(logL3,c(0,10),maximum=TRUE);opt3

$maximum
[1] 0.004508636

$objective
[1] -1193.491

lam = opt3$maximum;lam

[1] 0.004508636

alpha3 = n/(-1*sum(log(1-exp(-(lam*x)^2))));alpha3

[1] 0.3160101


# pdf
fr = function(x)
{
  fr = alpha3*(lam^2)*2*x*exp(-(lam*x)^2)*((1-exp(-(lam*x)^2))^(alpha3-1))
  return(fr)
}

# cdf

Fr = function(x)
{
  Fr = (1-exp(-(lam*x)^2))^(alpha3)
  return(Fr)
}

hist(y,probability = TRUE)
lines(y,dgamma(y,ga_shape,ga_rate),type="l",col=1)
lines(y,dweibull(y,alpha,1/beta^(1/alpha)),type="l",col=2)
lines(y,fee(y),type="l",col=3)
lines(y,fr(y),type="l",col=4)

 

A=c("Gamma","Weibull","Exponentiated Exponentail","Rayleigh")
legend(x="topright",legend = A,col=c(1,2,3,4),lty=1)

1-pweibull(50,alpha,1/beta^(1/alpha))
[1] 0.5579976

######     Survival Function     #######
plot(y,(1-Fr(y)),ylab="survival rate",xlab="failure stress",type="l")
 
########     Hazard Function      #######
plot(y,fr(y)/(1-Fr(y)),ylab="failure rate",xlab="failure stress",type="l")




 
## prac 2
1.
f1=function(abs,b,n){
  return(1-exp(-((n*abs)/b)))
}
n=seq(50,1000,by=50)
cp1=f1(0.01,2,n)
cp2=f1(0.02,2,n)
cp3=f1(0.03,2,n)
cp4=f1(0.04,2,n)
df=data.frame(n,cp1,cp2,cp3,cp4)
df
plot(n,cp1,data=df,ylab="Coverage Probability",type = "l")
lines(n,cp2,data=df,col=2)
lines(n,cp3,data=df,col=3)
lines(n,cp4,data=df,col=4)
A=c("0.01","0.02","0.03","0.04")
legend(x="bottomright",legend = c(A),col=c(1,2,3,4),lty=1,title="abs")
 
x=c(0.01,0.02,0.03,0.04)
mn=floor(-(2*log(0.05)/x))+1
mn

4.
n=50
k=1000
x1=matrix(0,k,n)
T={}
for(i in 1:k){
  for(j in 1:n){
  x1[i,j]=rexp(1,0.5)
     }
  T[i]=sqrt(n)*(mean(x1[i,])-2)/2
  }
hist(T,prob=TRUE)
x=seq(-3,3,by=0.01)
curve(dnorm(x),lwd=2,add=TRUE)
 
5.
n=25
k=1000
x1=matrix(0,k,n)
T={}
for(i in 1:k){
  for(j in 1:n){
    x1[i,j]=rpois(1,1)
  }
  T[i]=sqrt(n)*(mean(x1[i,])-1)
}
hist(T,prob=TRUE)
x=seq(-3,3,by=0.01)
curve(dnorm(x),lwd=2,add=TRUE)


## Preac 8
> df=read.csv(file.choose()) 
> n=length(df$IPG2211A2N);n [1] 397 
> plot.ts(df$IPG2211A2N,ylab="Electric Production",main="Monthly Electricity Production")
 
Seasonal Difference:
 v=var(Xt);v
 yt=diff(Xt,lag=12)
 v1=var(yt);v1
 yt2=diff(yt,lag=12)
 v2=var(yt2);v2

Trend Difference:
v1=var(yt);v1
> y1=diff(yt)
> v2=var(y1);v2
> y2=diff(y1)
> v3=var(y2):v3
> acf(y1,lag.max=40) # 2sqrt(n)=24
> pacf(y1,lag.max=40)
 
> Df=data.frame(
+ p=numeric(),
+ d=numeric(),
+ q=numeric(),
+ P=numeric(),
+ D=numeric(),
+ Q=numeric(),
+ AIC=numeric(),
+ Ljung_pvalue=numeric()
+ )
> p=1
> q=1
> d=1
> for(i in 1:2){
+ for(j in 1:2){
+ for(k in 1:2){
+
model=arima(Xt,order=c(p,d,q),seasonal=list(order=c(i,j,k),period=12,method="CSS",hessian=FALSE))
+ aic1=model$aic
+ residuals=resid(model)
+ result=Box.test(residuals,lag=24,type="Ljung-Box",fitdf=(p+q))
+ Df=rbind(Df,c(p,d,q,i,j,k,aic1,result$p.value))
+ }
+ }
+ }
> colnames(Df)=c("p","d","q","P","D","Q","AIC","Ljung_Pvalue")
> print(Df)
 p d q P D Q AIC Ljung_Pvalue
1 1 1 1 1 1 1 1787.226 2.567233e-03
2 1 1 1 1 1 2 1785.329 1.306633e-02
3 1 1 1 1 2 1 1881.790 1.387547e-10
4 1 1 1 1 2 2 1794.586 1.632148e-02
5 1 1 1 2 1 1 1773.643 9.346528e-02
6 1 1 1 2 1 2 1770.670 1.365475e-01
7 1 1 1 2 2 1 1832.470 3.859641e-03
8 1 1 1 2 2 2 1787.439 1.712291e-01
> d=1
> D=1
> p=1
> P=2
> q=1
> Q=2
> model1=arima(Xt,order=c(p,d,q),seasonal=list(order=c(P,D,Q),period=12));model1
> aic1=model1$aic;aic1
> # Residual Analysis
> residual=model1$residuals
> acf(residual,lag.max=24)
 
> result=Box.test(residual,lag=24,type="Ljung-Box",fitdf=p+q);result

Yuler Walker Equation ( AR(3) Model ):
> acf_values=acf(y1, lag.max = 40, plot = FALSE)
> r=acf_values$acf
> rho=matrix(c(r[2],r[3],r[4]),3,1,byrow=TRUE);rho
 [,1]
[1,] -0.2086338
[2,] -0.2513103
[3,] 0.1219564
> R=matrix(c(1,r[2],r[3],r[2],1,r[2],r[3],r[2],1),3,3,byrow=TRUE);R
 [,1] [,2] [,3]
[1,] 1.0000000 -0.2086338 -0.2513103
[2,] -0.2086338 1.0000000 -0.2086338
[3,] -0.2513103 -0.2086338 1.0000000
> B=solve(R)%*%rho;B
 [,1]
[1,] -0.27684620
[2,] -0.31170925
[3,] -0.01265101
 > ar_model_yw=ar(y1,aic = FALSE, order.max = 3,method = "yule-walker",);ar_model_yw
Call:
ar(x = y1, aic = FALSE, order.max = 3, method = "yule-walker")
Coefficients:
 1 2 3
-0.2768 -0.3117 -0.0127
Order selected 3 sigma^2 estimated as 9.56
> residuals=ar_model_yw$resid
> # Residual Analysis
> result=Box.test(residuals,lag=24,type="Ljung-Box",fitdf=3);result
 Box-Ljung test
data: residuals
X-squared = 91.889, df = 21, p-value = 7.605e-11


2.
> df=read.csv(file.choose());df
> n=length(df$Monthly.beer.production);n
[1] 476
> plot.ts(df$Monthly.beer.production,ylab="Beer Production",main="Monthly Beer Production")
 
> Xt=df$Monthly.beer.production
Seasonal Difference:
> v=var(Xt);v
[1] 1138.302
> yt=diff(Xt,lag=12)
> v1=var(yt);v1
[1] 154.3948
> yt2=diff(yt,lag=12)
> v2=var(yt2);v2
[1] 403.639
Trend Difference:
> v1=var(yt);v1
[1] 154.3948
> y1=diff(yt)
> v2=var(y1);v2
[1] 337.6726
> acf(yt,lag.max=44) # 2sqrt(n)=24
> pacf(yt,lag.max=44)
 
> Df=data.frame(
+ p=numeric(),
+ d=numeric(),
+ q=numeric(),
+ P=numeric(),
+ D=numeric(),
+ Q=numeric(),
+ AIC=numeric(),
+ Ljung_pvalue=numeric()
+ )
> p=2
> q=2
> d=1
> for(i in 1:2){
+ for(j in 1:2){
+ for(k in 1:2){
+ model=arima(Xt,order=c(p,d,q),seasonal=list(order=c(i,j,k),period=12,method="CSS",hessian=FALSE))
+ aic1=model$aic
+ residuals=resid(model)
+ result=Box.test(residuals,lag=44,type="Ljung-Box",fitdf=(p+q))
+ Df=rbind(Df,c(p,d,q,i,j,k,aic1,result$p.value))
+ }
+ }
+ }
> colnames(Df)=c("p","d","q","P","D","Q","AIC","Ljung_Pvalue")
> print(Df)
 p d q P D Q AIC Ljung_Pvalue
1 2 1 2 1 1 1 3445.058 1.110223e-16
2 2 1 2 1 1 2 3446.545 4.440892e-16
3 2 1 2 1 2 1 3533.369 0.000000e+00
4 2 1 2 1 2 2 3448.505 1.110223e-16
5 2 1 2 2 1 1 3445.641 2.553513e-15
6 2 1 2 2 1 2 3447.782 4.440892e-16
7 2 1 2 2 2 1 3506.473 0.000000e+00
8 2 1 2 2 2 2 3449.255 1.665335e-15
> d=1
> D=1
> p=2
> P=1
> q=2
> Q=1

> model1=arima(Xt,order=c(p,d,q),seasonal=list(order=c(P,D,Q),period=12));model1
Call:
arima(x = Xt, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = 12))
Coefficients:
 ar1 ar2 ma1 ma2 sar1 sma1
 -0.6844 -0.3024 -0.3693 -0.424 0.1146 -0.8567
s.e. 0.1122 0.0475 0.1123 0.101 0.0571 0.0324
sigma^2 estimated as 93.54: log likelihood = -1715.53, aic = 3445.06
> aic1=model1$aic;aic1
[1] 3445.058
> residual=model1$residuals
> acf(residual,lag.max=24)
 
> result=Box.test(residual,lag=24,type="Ljung-Box");result
 Box-Ljung test
data: residual
X-squared = 83.926, df = 24, p-value = 1.424e-08
Yuler Walker Equation ( AR(2) Model )
> acf_values=acf(yt, lag.max = 40, plot = FALSE)
> r=acf_values$acf
> rho=matrix(c(r[2],r[3]),2,1,byrow=TRUE);rho
 [,1]
[1,] -0.09167538
[2,] -0.09068540
> R=matrix(c(1,r[2],r[2],1),2,2,byrow=TRUE);R
 [,1] [,2]
[1,] 1.00000000 -0.09167538
[2,] -0.09167538 1.00000000
> B=solve(R)%*%rho;B
 [,1]
[1,] -0.10083647
[2,] -0.09992962
> ar_model_yw=ar(Xt1,aic = FALSE, order.max = 2,method = "yule-walker",);ar_model_yw
Call:
ar(x = Xt1, aic = FALSE, order.max = 2, method = "yule-walker")
Coefficients:
 1 2
-0.1008 -0.0999
Order selected 2 sigma^2 estimated as 152.2
> residuals=ar_model_yw$resid
> # Residual Analysis
> result=Box.test(residuals,lag=40,type="Ljung-Box",fitdf=2);result
 Box-Ljung test
data: residuals
X-squared = 271.82, df = 38, p-value < 2.2e-16

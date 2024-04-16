## prac 1
2.
block = factor(c(rep(1,4),rep(2,3),rep(3,5),rep(4,3)))
treat = factor(c(2,5,3,5,3,4,1,5,3,2,2,4,1,1,4))
y = c(10.3,14.2,9.8,11.4,8.9,12.3,13.3,13.5,10.1,11.3,10.9,12.5,12.8,13,12.7)
dat1 = data.frame(block,treat,y)
b = nlevels(block);b
[1] 4
v = nlevels(treat);v
[1] 5
N = xtabs(~treat+block,data=dat1);N    #Incidence Matrix

     block
treat 1 2 3 4
    1 0 1 0 2
    2 1 0 2 0
    3 1 1 1 0
    4 0 1 1 1
    5 2 0 1 0

T1 = aggregate(dat1$y,by = list(treat = dat1$treat),FUN=sum)
T = matrix(T1$x,v,1)
B1 = aggregate(dat1$y,by = list(block = dat1$block),FUN=sum)
B = matrix(B1$x,b,1)

#Treatment
xt = as.data.frame(table(treat));xt
  treat Freq
1     1    3
2     2    3
3     3    3
4     4    3
5     5    3

 a = xt$Freq;a
[1] 3 3 3 3 3
R = diag(a,v);R
     [,1] [,2] [,3] [,4] [,5]
[1,]    3    0    0    0    0
[2,]    0    3    0    0    0
[3,]    0    0    3    0    0
[4,]    0    0    0    3    0
[5,]    0    0    0    0    3

Rinv =solve(R)

#Block
xb = as.data.frame(table(block));xb
  block Freq
1     1    4
2     2    3
3     3    5
4     4    3

b1 = xb$Freq;b1
[1] 4 3 5 3

k = diag(b1,b);k
     [,1] [,2] [,3] [,4]
[1,]    4    0    0    0
[2,]    0    3    0    0
[3,]    0    0    5    0
[4,]    0    0    0    3

kinv =solve(k)
library(MASS)
c1 = R - N%*%kinv%*%t(N)
RankC = qr(c1)$rank
Rd = RankC+b;Rd
[1] 8

ginvC = ginv(c1)
Q = T-N%*%kinv%*%B

alpha_hat = ginvC%*%Q
beta_hat = kinv%*%(B-(t(N)%*%ginvC%*%Q))

t = v-RankC;t
[1] 1
n = length(y);n
[1] 15

if(RankC == v-1)
{
  cat("Since block design is connected hypotheses can be tested")
  cf = sum(y)^2/n;cf
  sst = t(y)%*%y - cf;sst
  ssb.unadj = t(B)%*%kinv%*%B - cf;ssb.unadj
  sstr.adj = t(alpha_hat)%*%Q;sstr.adj
  
  mstr = sstr.adj/(v-t)
  
  sse = sst - ssb.unadj - sstr.adj;sse
  mse = sse/(n-b-v+t)
  
  sstr.unadj = t(T)%*%Rinv%*%T - cf;sstr.unadj
  ssb.adj = sst - sstr.unadj - sse
  
  msb = ssb.adj/(b-1)
  
  F1 = mstr/mse
  F2 = msb/mse
  
  PT = 1 - pf(F1,v-t,(b-1)*(v-1))
  PB = 1 - pf(F2,b-t,(b-1)*(v-1))
  
  SV1 = c("Treatments(adj)","Blocks(unadj)","Error","Total") 
  SV2 = c("Treatments(unadj)","Blocks(adj)","Error","Total") 
  
  df1 = c(v-t,b-1,n-v-b+t,n-1)
  df2 = c(v-1,b-t,n-v-b+t,n-1)
  
  ss1 =round(c(sstr.adj,ssb.unadj,sse,sst),4)
  ss2 =round(c(sstr.unadj,ssb.adj,sse,sst),4)
  
  ms1 = c(round(mstr,4),"-",round(mse,4),"-")
  ms2 = c("-",round(msb,4),round(mse,4),"-")
  
  f_1 = c(round(F1,4),"-","-","-")
  f_2 = c("-",round(F2,4),"-","-")
  
  p_1 = c(round(PT,4),"-","-","-")
  p_2 = c("-",round(PB,4),"-","-")
  
  print("H_alpha")
  print(data.frame(SV1,df1,ss1,ms1,f_1,p_1))
  
  print("H_beta")
  print(data.frame(SV2,df2,ss2,ms2,f_2,p_2))
}


Since block design is connected hypotheses can be tested[1] "H_alpha"
              SV1 df1     ss1   ms1    f_1    p_1
1 Treatments(adj)   4 24.6679 6.167 8.8856 0.0014
2   Blocks(unadj)   3  4.1338     -      -      -
3           Error   7  4.8583 0.694      -      -
4           Total  14 33.6600     -      -      -
[1] "H_beta"
                SV2 df2     ss2    ms2    f_2    p_2
1 Treatments(unadj)   4 27.9200      -      -      -
2       Blocks(adj)   3  0.8817 0.2939 0.4235 0.7397
3             Error   7  4.8583  0.694      -      -
4             Total  14 33.6600      -      -      -



# ANCOVA

1.
Hypothesis to be tested is,
H0: β = 0  Vs  H1: β ≠ 0
y=c(13,14,12,12,14,12,10,11,12,14,15,14,11,11,10)
x=c(46.5,45.9,49.8,46.1,44.3,48.7,49.0,50.1,48.5,45.2,46.3,47.1,48.9,48.2,50.3)
r=5
t=3
n=r*t;n
[1] 15
ym=matrix(y,ncol=3,byrow = FALSE)
xm=matrix(x,ncol=3,byrow = FALSE)
csx=colSums(xm)
csy=colSums(ym)
rsx=rowSums(xm)
rsy=rowSums(ym)
Correction factors,
cfyy=sum(y)^2/n;cfyy
[1] 2281.667

cfxx=sum(x)^2/n;cfxx
[1] 34072.13
cfxy=sum(x)*sum(y)/n;cfxy
[1] 8817.1
Total sum of squares,
TSSxy=sum(x*y)-cfxy;TSSxy
[1] -33.7
TSSyy=sum(y^2)-cfyy;TSSyy
[1] 35.33333
TSSxx=sum(x^2)-cfxx;TSSxx
[1] 49.296
Treatment sum of squares,
Tryy=(sum(csy^2)/r)-cfyy;Tryy
[1] 3.733333
Trxx=(sum(csx^2)/r)-cfxx;Trxx
[1] 9.796
Trxy=(sum(csx*csy)/r)-cfxy;Trxy
[1] -5.84
Error sum of squares,
 Eyy=TSSyy-Tryy;Eyy
[1] 31.6
Exx=TSSxx-Trxx;Exx
[1] 39.5
 Exy=TSSxy-Trxy;Exy
[1] -27.86
SSE=Eyy-(Exy^2/Exx);SSE
[1] 11.94988
Mean squares error,
MSE=SSE/(n-t-1);MSE
[1] 1.086353
F statistic:
F_cal=(Exy^2/Exx)/MSE;F_cal          
[1] 18.08815
qf(0.05,1,11,lower.tail=FALSE)
[1] 4.844336
Conclusion: Here we observe that Fcal >Ftab, so we reject null hypothesis and conclude that there is effect of covariate on the response variable.

To test H0: α(i) = 0,
Sum of squares adjusted for regression,
 SSE2=TSSyy-(TSSxy^2/TSSxx);SSE2
[1] 12.29516
F0=((SSE2-SSE)/(t-1))/MSE;F0
[1] 0.1589135
qf(0.05,2,11,lower.tail=FALSE)
[1] 3.98




2.
Hypothesis to be tested is,
H0: β = 0  Vs  H1: β ≠ 0
x=c(5,3.5,4.7,4.3,3.8,3,4.2,4.5,4.3,4.1,5.5,4.8,4.8,8.0,7.4,7.1,6.7,5.6,8.1,8.7,8.3,8.4,7.3,8.5,7.8,8.0,8.4,7.7,6.3,8.6,11.8,12.4,12,11.4,10.4,9.2,9,8,7.3,7.2,6.1,6.4,6.4,6.9,5.8,13,13.3,10.7,12.2,11.6,9.7,7,6,7.1,5.3,6.2,9.4,10.6,9,7.6)
y=c(263.7,130.8,382.9,302.5,213.3,132.8,292,315.5,262.4,314.4,310.8,280,331.7,672.5,496,566.7,552.5,397.5,532.3,587.8,520.9,574.3,505,604.6,522.5,555,561.1,791.7,610,710,940.7,990,916.2,935,724.3,611.1,621.7,731,710,604.7,508.8,393,416,400,335.6,983.3,958.8,747.8,866,810.8,950,485.4,395.4,465.4,371.4,402,837.1,901.2,595.7,510)
treat=factor(c(rep('A',15),rep('B',12),rep("c",10),rep("D",8),rep("E",6),rep("F",5),rep("G",4)))
data=data.frame(x,y,treat)
t=7
n=length(x);n
[1] 60
Correction factors,
 cfyy=sum(y)^2/n;cfyy
[1] 18494378
cfxx=sum(x)^2/n;cfxx
[1] 3471.683
cfxy=sum(x)*sum(y)/n;cfxy
[1] 253390.2

Total sum of squares,
TSSxy=sum(x*y)-cfxy;TSSxy
[1] 32007.08
TSSyy=sum(y^2)-cfyy;TSSyy
[1] 3003704
TSSxx=sum(x^2)-cfxx;TSSxx
[1] 388.7173
Treatment sum of squares,
Tryy=sum((aggregate(data$y,by=list(data$treat),FUN=mean)['x'])*(aggregate(data$y,by=list(data$treat),FUN=sum)['x']))-cfyy;Tryy
[1] 2218023
Trxx=sum((aggregate(data$x,by=list(data$treat),FUN=mean)['x'])*(aggregate(data$x,by=list(data$treat),FUN=sum)['x']))-cfxx;Trxx
[1] 297.1311
Trxy=sum((aggregate(data$y,by=list(data$treat),FUN=sum)['x'])*(aggregate(data$x,by=list(data$treat),FUN=mean)['x']))-cfxy;Trxy
[1] 25498.99
Error sum of squares,
 Eyy=TSSyy-Tryy;Eyy
[1] 785680.8
 Exx=TSSxx-Trxx;Exx
[1] 91.58625
Exy=TSSxy-Trxy;Exy
[1] 6508.097
 SSE=Eyy-(Exy^2/Exx);SSE
[1] 323217
Mean squares error,
 MSE=SSE/(n-t-1);MSE
[1] 6215.711
F statistic:
 F_cal=(Exy^2/Exx)/MSE;F_cal
[1] 74.4024
Conclusion: Here we observe that Fcal >Ftab, so we reject null hypothesis and conclude that there is effect of covariate on the response variable.

To test H0: α(i) = 0,
Sum of squares adjusted for regression,
SSE2=TSSyy-(TSSxy^2/TSSxx);SSE2
[1] 368232.7
 F0=((SSE2-SSE)/(t-1))/MSE;F0
[1] 1.207042


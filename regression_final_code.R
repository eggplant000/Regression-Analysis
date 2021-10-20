

#### 회귀분석 프로젝트
#### 
#### 20176735 김수린
#### 2021.06.20.
####

# 데이터 정제 ------------------------------------------------------------------


## load data
raw_data = read.csv("data.csv", header = T)



## 편의를 위해 열이름 변경
names(raw_data) <- c('date', 'day', 'aud', 'avgtem', 'lowtem', 
                     'hitem', 'rain', 'avgWS', 'maxWS', 
                     'avgHum', 'sumSun','sumSR', 'maxSun', 
                     'tem15', 'maxSn', 'avgCl', 'dust')
data = raw_data
head(data)



## 관객수의 크기가 크므로 (천 명) 단위 사용
data['aud'] = data['aud']/1000



## 결측치 확인
apply(is.na(data), 2, sum) # rain 226개, maxSn 359개, dust 27개
                           # sumSun 1개, sumSR 2개, maxSun 1개, tem15 1개



#### 결측치가 적은 변수는 해당 행 삭제
## 240, 241, 284번째 데이터 삭제
data[c(240, 241, 284),]
data = data[-c(240, 241, 284),]
rownames(data) <- NULL  # 인덱스 초기화



#### 결측치가 많은 변수 제외/대체
## maxSn, dust: 제외
## rain: 대체

#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
attach(data)



## 일강수량과 평균상대습도의 상관관계 (r = .55, p < .001)
cor.test(data[!is.na(rain), "rain"], data[!is.na(rain), "avgHum"])
chart.Correlation(data[!is.na(rain), c("rain", 'avgHum')], histogram = FALSE)



## 결측치가 많은 변수 삭제
library(dplyr)
data = data %>% select(-c(rain, maxSn, dust)) 



## 모든 독립변수들 간의 상관관계 확인
chart.Correlation(data[,4:14], histogram=FALSE)



## 각 독립변수와 종속변수의 단순회귀분석 결과를 통해 
## 상관관계가 높은 변수 중 하나만 선택
summary(lm(sqrt(aud)~avgtem, data=data)) #0.308
summary(lm(sqrt(aud)~lowtem, data=data)) #0.198
summary(lm(sqrt(aud)~hitem, data=data)) #0.487

summary(lm(sqrt(aud)~avgWS, data=data)) #0.0319*
summary(lm(sqrt(aud)~maxWS, data=data)) #0.00495**

summary(lm(sqrt(aud)~sumSun, data=data)) #0.18
summary(lm(sqrt(aud)~sumSR, data=data)) #0.22
summary(lm(sqrt(aud)~maxSun, data=data)) #0.104




# 모형 1 (완전 모형) ------------------------------------------------------------



m1 = lm(aud ~ lowtem + maxWS +  avgHum +
                maxSun + tem15 + avgCl, data=data)
summary(m1)


## 회귀분석의 가정 확인
par(mfrow = c(2, 2))
plot(m1) 
# 분산의 동질성 가정이 성립되지 않음



# 모형 2 (sqrt 변환) ----------------------------------------------------------



m2 = lm(sqrt(aud) ~ lowtem + maxWS + avgHum +
                maxSun + tem15 + avgCl, data=data)
summary(m2)

X = model.matrix(m2)
n = nrow(X)
p = ncol(X)


## 회귀분석의 가정 확인
par(mfrow = c(2, 2))
plot(m2) 

## (y_hat, r)
par(mfrow = c(1, 1))
res2 = rstandard(m2)
plot(m2$fitted.values, res2, 
     xlab = 'Fitted values', ylab='Residuals') 

## (x, r)
par(mfrow = c(2, 3))
plot(hitem, res2, xlab="Max temperature", ylab="Residuals")
plot(maxWS, res2, xlab="Max wind speed", ylab="Residuals")
plot(maxSun, res2, xlab="Max solar radiation quantity in a hour", ylab="Residuals")
plot(avgHum, res2, xlab="Avearage relative Humidity", ylab="Residuals")
plot(tem15, res2, xlab="1.5m soil temperature", ylab="Residuals")
plot(avgCl, res2, xlab="Average amount of clouds", ylab="Residuals")
# -> 분산의 동질성 가정 만족


## 영향점 확인
par(mfrow = c(1, 1))

# 1.leverage point
pii = influence(m2)$hat
lev.idx = pii > 2*p/n

color = rep("black", len=n)
color[lev.idx] = 'red'

plot(pii, res2, type='n')
points(pii, res2, pch=21, cex=0.8, bg=color)

# 2.cooks distance
Ci = cooks.distance(m2)
plot(pii, Ci, type='n')
points(pii, Ci, pch=21, cex=0.8, bg=color)
sum(Ci > qf(0.5,p,n-p)) #없음

# 3.dfits
dfits = dffits(m2)
inf.idx = abs(dfits) >= 2*sqrt(p/(n-p))
color[inf.idx] = 'blue'

plot(pii, dfits, type='n', 
     xlab = 'Leverage', ylab = 'DFITS')
points(pii, dfits, pch=21, cex=0.8, bg=color)
text(dfits[which.max(pii)]~pii[which.max(pii)], 
     labels=names(lev.idx[which.max(pii)]),
     cex=0.9, font=2, pos=1)
# index 248 데이터가 다른 데이터들과 확연히 차이나는 영향점임을 확인



# 모형 3 (영향점 제거) -----------------------------------------------------------


data.rem = data[-248,] # 영향점 제거
rownames(data.rem) <- NULL

m3 = lm(sqrt(aud) ~ lowtem + maxWS + maxSun + 
          avgHum + tem15 + avgCl, data=data.rem)
summary(m3)


## 회귀분석의 가정 확인
par(mfrow = c(2, 2))
plot(m3) # 만족


## 영향점 확인
par(mfrow = c(1, 1))

# 1.leverage point
pii = influence(m3)$hat
lev.idx = pii > 2*p/n
color = rep("black", len=n)
color[lev.idx] = 'red'

# 2.cooks distance
Ci = cooks.distance(m3)
sum(Ci > qf(0.5,p,n-p)) #없음

# 3.dfits
dfits = dffits(m3)
inf.idx = abs(dfits) >= 2*sqrt(p/(n-p))
color[inf.idx] = 'blue'

plot(pii, dfits, type='n', 
     xlab = 'Leverage', ylab = 'DFITS')
points(pii, dfits, pch=21, cex=0.8, bg=color)
# 다른 데이터들과 확연히 차이나는 영향점은 없음


## 자기상관(Autocorrelation) 확인
# Durbin-Watson test
library(lmtest)
dwtest(m3, alternative = 'greater')

# runs test
library(randtests)
runs.test(m3$residuals, alternative="left.sided", plot = FALSE)



# 모형 4 (AR(1)) ------------------------------------------------------------



## 표본자기상관계수
acf(m3$residuals, plot=TRUE)
acf(m3$residuals, plot=FALSE)
rho <- acf(m3$residuals, plot=FALSE)[1]$acf[1]

n <- nrow(data.rem)

# Fit a linear model with transformed variables
taud <- sqrt(data.rem$aud[2:n]) - rho*sqrt(data.rem$aud[1:(n-1)])
thitem <- data.rem$hitem[2:n] - rho*data.rem$hitem[1:(n-1)]
tmaxWS <- data.rem$maxWS[2:n] - rho*data.rem$maxWS[1:(n-1)]
tmaxSun <- data.rem$maxSun[2:n] - rho*data.rem$maxSun[1:(n-1)]
tavgHum <- data.rem$avgHum[2:n] - rho*data.rem$avgHum[1:(n-1)]
ttem15 <- data.rem$tem15[2:n] - rho*data.rem$tem15[1:(n-1)]
tavgCl <- data.rem$avgCl[2:n] - rho*data.rem$avgCl[1:(n-1)]
m4 <- lm(taud ~ thitem + tmaxWS + tmaxSun +
           tavgHum + ttem15 + tavgCl)
summary(m4)

acf(m4$residuals, plot=FALSE)

dwtest(m4, alternative = 'greater')
runs.test(m4$residuals, alternative="left.sided", plot = FALSE)



# 모형 5 (요일변수 추가) ----------------------------------------------------------


m5 = lm(sqrt(aud) ~ lowtem + maxWS + maxSun + 
          avgHum + tem15 + avgCl + as.factor(day), data=data.rem)

summary(m5)


X = model.matrix(m5)
n = nrow(X)
p = ncol(X)


## 회귀분석의 가정 확인
par(mfrow = c(2, 2))
plot(m5) # 만족


## 영향점 확인
par(mfrow = c(1, 1))

# 1.leverage point
pii = influence(m5)$hat
lev.idx = pii > 2*p/n
color = rep("black", len=n)
color[lev.idx] = 'red'

# 2.cooks distance
Ci = cooks.distance(m5)
sum(Ci > qf(0.5,p,n-p)) #없음

# 3.dfits
dfits = dffits(m5)
inf.idx = abs(dfits) >= 2*sqrt(p/(n-p))
color[inf.idx] = 'blue'

plot(pii, dfits, type='n', 
     xlab = 'Leverage', ylab = 'DFITS')
points(pii, dfits, pch=21, cex=0.8, bg=color)
# 다른 데이터들과 확연히 차이나는 영향점은 없음



## 모형 선택

library(olsrr)

# 전진선택법
mfwd_aic <- ols_step_forward_aic(m5, details=FALSE)
mfwd_aic$model

# 후진소거법
mbwd_aic <- ols_step_backward_aic(m5, details=TRUE)
mbwd_aic$predictors

# 단계적 선택법
mboth_aic <- ols_step_both_aic(m5, details=TRUE)
mboth_aic$predictors

# 셋 다 같은 결과를 나타냄
# 평균상대습도(avgHum), 1.5m 지중온도(tem15), 평균전운량(avgCl) 삭제

final = mfwd_aic$model
summary(final)



## ANOVA test

final2 = lm(sqrt(aud) ~ maxWS + maxSun + as.factor(day), data=data.rem)

anova(final2, final) # NOT reject
summary(final2)

final3 = lm(sqrt(aud) ~ maxWS + as.factor(day), data=data.rem)
anova(final3,final2) # reject



# 최종 모형 -------------------------------------------------------------------



summary(final2)

# 95% 신뢰구간
confint(final2) 


X = model.matrix(final2)
n = nrow(X)
p = ncol(X)


## 회귀분석의 가정 확인
par(mfrow = c(2, 2))
plot(final2)

## (y_hat, r)
par(mfrow = c(1, 1))
res = rstandard(final2)
plot(final2$fitted.values, res, 
     xlab = 'Fitted values', ylab='Residuals') 

## (x, r)
par(mfrow = c(2, 3))
plot(data.rem$hitem, res, xlab="Max temperature", ylab="Residuals")
plot(data.rem$maxWS, res, xlab="Max wind speed", ylab="Residuals")
plot(data.rem$maxSun, res, xlab="Max solar radiation quantity in a hour", ylab="Residuals")
plot(data.rem$avgHum, res, xlab="Avearage relative Humidity", ylab="Residuals")
plot(data.rem$tem15, res, xlab="1.5m soil temperature", ylab="Residuals")
plot(data.rem$avgCl, res, xlab="Average amount of clouds", ylab="Residuals")



## 영향점 확인
par(mfrow = c(1, 1))

# 1.leverage point
pii = influence(final2)$hat
lev.idx = pii > 2*p/n
color = rep("black", len=n)
color[lev.idx] = 'red'

# 2.cooks distance
Ci = cooks.distance(final2)
sum(Ci > qf(0.5,p,n-p)) #없음

# 3.dfits
dfits = dffits(final2)
inf.idx = abs(dfits) >= 2*sqrt(p/(n-p))
color[inf.idx] = 'blue'

plot(pii, dfits, type='n', 
     xlab = 'Leverage', ylab = 'DFITS')
points(pii, dfits, pch=21, cex=0.8, bg=color)
# 다른 데이터들과 확연히 차이나는 영향점은 없음


## multicollinearity
library(car)
vif(final2)

## autocorrelation
dwtest(final2, alternative = 'greater')
runs.test(final2$residuals, alternative="left.sided", plot = FALSE)

acf(final2$residuals, plot=TRUE)
acf(final2$residuals, plot=FALSE)



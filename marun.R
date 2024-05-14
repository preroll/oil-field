library(data.table)
library(ggplot2)
library(caTools)
library(car)
library(rpart)
library(rpart.plot)


marun <- fread('marun_sample2.csv')

#=============================================================================================#
                                      #DATA CLEANING#
#=============================================================================================#

#find out existence of NA values
summary(marun)
#check categories for Formation
unique(marun$Formation)
unique(marun$`Hole size (in)`)

#non NA set of data
non.na <- marun[-c(3,6)]
#check histograms for FAN300 & FAN600
hist(non.na$FAN300)
hist(non.na$FAN600)
hist(non.na$MIN10GEL)
hist(non.na$MUDLOSSU)
#find medians of FAN300 & FAN600 to input in NA row since data is skewed to left
median(non.na$FAN300) #30
median(non.na$FAN600) #49
median(non.na$MIN10GEL) #5
median(non.na$MUDLOSSU) #25
#input medians since frequency distribution is skewed
marun$FAN300[is.na(marun$FAN300)] <- median(non.na$FAN300)
marun$FAN600[is.na(marun$FAN600)] <- median(non.na$FAN600)
marun$MIN10GEL[is.na(marun$MIN10GEL)] <- median(non.na$MIN10GEL)
marun$MUDLOSSU[is.na(marun$MUDLOSSU)] <- median(non.na$MUDLOSSU)
#check correlations again
cor(marun)

#=============================================================================================#
                                     #DATA EXPLORATION#
#=============================================================================================#

#setting categories for categorical variable Formation
marun$Formation <- factor(marun$Formation) 

#PLOTS USED IN CBA QUESTION 1
#pump pressure ~ pump flow rate (y ~ x) group by hole size
ggplot(marun, aes(y=`Pump flow rate`, x=`Pump pressure`, col=`Hole size (in)`)) + geom_smooth() + geom_point()
#larger hole size, larger pump flow rate, peak ~ 1000 pfr
ggplot(marun, aes(x=factor(marun$`Hole size (in)`), y=`Pump flow rate`, fill=`Hole size (in)`)) + geom_violin() + labs(x='Hole size (in)')
#lost circulation ~ fracture pressure (y ~ x) group by hole size
ggplot(marun, aes(x=`Fracture pressure`, y=`MUDLOSSU`, col=`Hole size (in)`)) + geom_point() + labs(y='Severity of lost circulation')
#viscosity ~ formation (y ~ x)
ggplot(marun, aes(x=Formation, y=MFVIS, fill=Formation)) + geom_violin() + labs(y='Viscosity')

#TRIAL PLOTS WHILE EXPLORING DATA
#northing ~ easting (y ~ x) to visualise the locations of the oil wells
ggplot(marun, aes(x=Easting, y=Northing)) + geom_smooth() + geom_point() 
#FAN600 ~ FAN300 (y ~ x) positive linear
ggplot(marun, aes(x=FAN300, y=FAN600)) + geom_smooth(method=lm)
#viscosity against pump flow rate wrt hole size
ggplot(marun, aes(x=MFVIS, y=`Pump flow rate`, col=`Hole size (in)`)) + geom_point()
#pressure correlations are postive and highly correlated
ggplot(marun, aes(x=`Pore pressure`, y=`Mud pressure (psi)`)) + geom_smooth()
#meterage ~ drill time wrt formation
ggplot(marun, aes(x=DRLTIME, y=METERAGE, col=Formation)) + geom_jitter()
#depth ~ formation
ggplot(marun, aes(x=Formation, y=`Depth (ft)`, fill=Formation)) + geom_violin()

#=============================================================================================#
                                    #LINEAR REGRESSION#
#=============================================================================================#

set.seed(2)
#train test split
train <- sample.split(Y = marun$MUDLOSSU, SplitRatio = 0.7)
trainset <- subset(marun, train == T)
testset <- subset(marun, train == F)

#checking the distribution of MUDLOSSU similar in trainset vs testset
summary(trainset$MUDLOSSU)
summary(testset$MUDLOSSU)

#finding best model
model1 <- lm(MUDLOSSU ~ ., data = trainset)
summary(model1) #0.4292
model2 <- lm(MUDLOSSU ~ .-`Pump pressure`-DRLTIME-RPM-`Depth (ft)`, data = trainset)
summary(model2) #0.4301
model3 <- lm(MUDLOSSU ~ .-`Pump pressure`-DRLTIME-RPM-`Depth (ft)`-`Mud pressure (psi)`-MIN10GEL-Northing-`Pore pressure`-Easting-WOB, data = trainset)
summary(model3) #0.4278
vif(model3) #high vif for FAN300 & FAN600
model4 <- lm(MUDLOSSU ~ .-`Pump pressure`-DRLTIME-RPM-`Depth (ft)`-`Mud pressure (psi)`-MIN10GEL-Northing-`Pore pressure`-Easting-WOB-FAN300, data = trainset)
summary(model4) #0.4231
#check multicollinearity
vif(model4) 

#chose model3 and finding trainset rmse
residuals(model4) 
rmse.model4train <- sqrt(mean(residuals(model4)^2))
summary(abs(residuals(model4)))
#predict on testset
predict.model4test <- predict(model4, newdata = testset)
testset.error <- testset$MUDLOSSU - predict.model4test
#testset rmse
rmse.model4test <- sqrt(mean(testset.error^2))
summary(abs(testset.error))

#Q3 rmse of LIN REG train & test
rmse.model4train #123.7041
rmse.model4test #123.5474

#=============================================================================================#
                            #CLASSIFICATION AND REGRESSION TREE#
#=============================================================================================#

set.seed(2)
#train test split
train <- sample.split(Y = marun$MUDLOSSU, SplitRatio = 0.7)
trainset <- subset(marun, train == T)
testset <- subset(marun, train == F)

#build initial cart model
cart1 <- rpart(MUDLOSSU~., method = 'anova', data=trainset, 
               control=rpart.control(minsplit=2, cp=0))
plotcp(cart1)

dt <- data.table(printcp(cart1))
#number the sequence of the trees
dt[, index := 1:nrow(dt)]
#which index has the minimum error (xerror+xstd)
min_cp_index = min(dt[xerror+xstd == min(xerror+xstd), index])
#find out the error cap 
errorcap <- dt[min_cp_index, xerror+xstd]
errorcap
#find optimal index (which crosses the horizontal line from the left)
optimal_cp_index <- min(dt[xerror<errorcap, index])
optimal_cp_index
#find out optimal CP which is the geometric mean of this index CP and one below it 
cp.opt = sqrt(dt[index == optimal_cp_index, CP] * dt[index == optimal_cp_index-1, CP])
cp.opt
#get the trained model
m.opt <- prune(cart1, cp=cp.opt)
rpart.plot(m.opt)

#trainset rmse
rmse.cart1train <- sqrt(mean(residuals(m.opt)^2))
summary(abs(residuals(m.opt)))
#predict on testset
predict.cart1test<- predict(m.opt, newdata=testset)
testset.error <- testset$MUDLOSSU - predict.cart1test
#testset rmse
rmse.cart1test <- sqrt(mean(testset.error^2))
summary(abs(testset.error))

#Q3 rmse of CART train & test
rmse.cart1train
rmse.cart1test

#summary of optimal CART
summary(m.opt)

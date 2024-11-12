cat_data <- read.csv("~/Desktop/DataAnalytics2022/AdvancedStatisticsUsingR/CA1/cat_hwt.csv")
cat_data
typeof(cat_data)
summary(cat_data)
str(cat_data)
attach(cat_data)
install.packages("rmarkdown")
install.packages("knitr")

library(stats)
library(car)

variables.names <- names(cat_data)
variables.names

#Data was collected on deceased male and female adult cats used for experiments:
#Sex: 1 for male and 2 for female. Nominal Categorical
#Bwt: body weight in kg. Continuous Numerical 
#Hwt: heart weight in g. Continuous Numerical 
#Height: Height in cm. Continuous Numerical 
#Age: Age in years. Continuous Numerical 
#Outdoor: kept outdoors; 1 = Always, 2 = Frequently, 3 = Never. Nominal Categorical

#Q1A 
#Investigate if any of the continuous numerical variables have a linear relationship by producing scatterplots,
#Interpret these plots.

#Scatter plot of Cat's Body Weight v Cat's Heart Weight
plot(Bwt,Hwt, xlab="Cat's Body Weight (kg)", ylab="Cat's Heart Weight (g)", main="Scatter plot of Cat's Body Weight v Cat's Heart Weight", pch=16)
#Plot regression line
abline(lm(Hwt~Bwt,data=cat_data),col='red')
#From this scatter plot we can see the following relationship between Cat's Body Weight v Cat's Heart Weight 
#The data points don't stray to far away from the regression line
#Cat's Body Weight v Cat's Heart Weight appear to have a very strong positive association because in general Cat's Body Weight increases when Cat's Heart Weight increases.

#Scatter plot of Cat's Body Weight v Cat's Height
plot(Bwt,Height, xlab="Cat's Body Weight (kg)", ylab="Cat's Height (cm)", main="Scatter plot of Cat's Body Weight v Cat's Height", pch=16)
#Plot regression line
abline(lm(Height~Bwt,data=cat_data),col='red')
#From this scatter plot we can see the following relationship between Cat's Body Weight v Cat's Height 
#The data points are very scattered and dispersed but with the regression line we can interpret this scatter plot better.
#The data points are greatly dispersed away from the away from the regression line
#Cat's Body Weight v Cat's Height appear to have a weak positive association because in general Cat's Body Weight slightly increases when Cat's Height slightly increases.

#Scatter plot of Cat's Body Weight v Cat's Age
plot(Bwt,Age, xlab="Cat's Body Weight (kg)", ylab="Cat's Age (years)", main="Scatter plot of Cat's Body Weight v Cat's Age", pch=16)
#Plot regression line
abline(lm(Age~Bwt,data=cat_data),col='red')
#From this scatter plot we can see the following relationship between Cat's Body Weight v Cat's Age
#The data points are greatly dispersed away from the away from the regression line
#The data points are very scattered and dispersed but with the regression line we can interpret this scatter plot better.
#Cat's Body Weight v Cat's Age appear to have a weak positive association because in general Cat's Body Weight slightly increases when Cat's Age slightly increases.

#Scatter plot of Cat's Heart Weight v Cat's Height
plot(Hwt,Height, xlab="Cat's Heart Weight (g)", ylab="Cat's Height (cm)", main="Scatter plot of Cat's Heart Weight v Cat's Height", pch=16)
#Plot regression line
abline(lm(Height~Hwt,data=cat_data),col='red')
#From this scatter plot we can see the following relationship between Cat's Heart Weight v Cat's Height
#The data points are moderately dispersed away from the away from the regression line
#Cat's Heart Weight v Cat's Height appear to have a medium positive association because in general Cat's Heart Weight increases when Cat's Height increases.


#Scatter plot of Cat's Heart Weight v Cat's Age
plot(Hwt,Age, xlab="Cat's Heart Weight (g)", ylab="Cat's Age (years)", main="Scatter plot of Cat's Heart Weight v Cat's Age", pch=16)
#Plot regression line
abline(lm(Age~Hwt,data=cat_data),col='red')
#From this scatter plot we can see the following relationship between  Cat's Heart Weight v Cat's Age
#The data points are dispersed away from the away from the regression line
#Cat's Heart Weight v Cat's Age appear to have a very weak association because in general Cat's Heart Weight increases when Cat's Age increases.

#Scatter plot of Cat's Height v Cat's Age
plot(Height,Age, xlab="Cat's Height (cm)", ylab="Cat's Age (years)", main="Scatter plot of Cat's Height v Cat's Age", pch=16)
#Plot regression line
abline(lm(Age~Height,data=cat_data),col='red')
#From this scatter plot we can see the following relationship between Cat's Height v Cat's Age
#The data points are very scattered and dispersed but with the regression line we can interpret this scatter plot better.
#The data points are greatly dispersed away from the away from the regression line
#Cat's Height v Cat's Age appear to have a very weak positive association because in general Cat's Height slightly increases when Cat's Age slightly increases.

#Correlation Matrix
cor(cat_data[,c("Hwt","Bwt", "Height", "Age")])

panel.cor <- function(x, y, digits = 2, prefix = "",
                      cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), 
                digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (abs(r)+0.15))
}

pairs(~Hwt+Bwt+Height+Age,cat_data, upper.panel = panel.smooth, lower.panel= panel.cor)
pairs(~Hwt+Bwt+Height+Age,cat_data, upper.panel = panel.smooth)
#Looking at this plot, there are many linear and non-linear relationships
#Hwt & Bwt has a strong positive linear relationship
#Hwt & Height has a weak positive linear relationship
#Hwt & Age has a very weak positive linear relationship

#Bwt & Height has a weak positive linear relationship
#Bwt & Age has a weak positive linear relationship

#Height & Age has a weak positive curved relationship

#Q1B 
#What is the response variable and why? 
#What is the research question of interest in this dataset?

#Response variable is "Hwt" heart weight in grams because this is something that can't be measured easily.
#Research question: What is the predicted heart weight of a cat? 
#We are interested in exploring the variables to see if they have a significant impact on a cat's heart weight and can we use these variables to predict a cat's heart weight

#Q1C 
#Investigate if any of the categorical /discrete variables seem to have a relationship with the response variable (selected in Q1B)  using boxplots. 
#Interpret these plots. 

#Boxplot of Cat's Heart Weight by Sex
MaleHwt <- Hwt[Sex==1]
MaleHwt
summary(MaleHwt)
FemaleHwt <- Hwt[Sex==2]
FemaleHwt
summary(FemaleHwt)
boxplot(MaleHwt, FemaleHwt, names = c("Male Cat", "Female Cat"), main="Comparing Heart Weight's of Male and Female Cats", ylab="Heart Weight (g)", col=c(4,2))
#From this box plot of Cat's Heart Weight by Cat's Sex, we can see the following information:
#Male cats tend to have a higher heart weight than female cats
#As male cat's heart weight has a higher median and IQR than female cats
#There is a wider spread for male cat's heart weight than female cats
#There s an outlier present in both male and female groups

#Boxplot of Cat's Heart Weight by Cat's Outdoor Categories 
AlwaysHwt <- Hwt[Outdoor==1]
AlwaysHwt
summary(AlwaysHwt)
FreqHwt <- Hwt[Outdoor==2]
FreqHwt
summary(FreqHwt)
NeverHwt <- Hwt[Outdoor==3]
NeverHwt
summary(NeverHwt)
boxplot(AlwaysHwt, FreqHwt, NeverHwt, names = c("Always Outdoors", "Frequently Outdoors", "Never Outdoors"), main="Comparing Heart Weight's by Cat's Outdoor Categories", ylab="Heart Weight (g)", col=c(4,2,3))
#From this box plot of Cat's Heart Weight by Cat's Outdoor Categories, we can see the following information:
#Cat's that are kept always outdoors tend to have a lower heart weight
#As cat's heart weight has a higher median and IQR if the cat was never kept outdoors
#The spreads are similar across all boxplots and median roughly in middle
#There is one outlier present in cat's that are never kept outdoors

#Q1D 
#Using part Q1C, choose one categorical /discrete variable with at least 3 categories to test to see if there is any difference between the means of the response variable (selected in part b) using an one-way ANOVA test. 
#(Note: if categorical variables are not set up as factors, so change to factor prior to running analysis (as.factor())

#Selecting "Outdoor" variable to run a one way anova on "Hwt" heart weight
#As seen from the boxplot above there seems to be a difference in the means
#But also the variability between groups seems similar so it is appropriate to run the one way anova
#We will check assumptions

#Null hypothesis no difference in the true mean of "Hwt" Cat's heart weight for each "Outdoor" group 
#H0 mu_1=mu_2=mu_3, where mu_i is the true  mean of "Hwt" Cat's heart weight for group i 
#Alternative hypothesis at least one Hwt" Cat's heart weight mean is different for each "Outdoor" group 
#H1 not all mu_i are the same, where mu_i is the true  mean of "Hwt" Cat's heart weight for group i 

res1 <- aov(Hwt~as.factor(Outdoor), data = cat_data)
par(mfrow=c(2,1))
plot(res1, which=1:2)
#Interpret Residuals vs Fitted 
#The spread of the three groups looks roughly the same so we are happy with our equal variance assumption. 
#The scatter seems randomly distributed within in the 3 groups, so we are happy with our iid assumption.
#Interpret Normal Q-Q Plot
#The points fall roughly in a straight line indicating that the residuals are roughly normally distributed.

summary(res1)

#p-value = 4.71e-16 which is less than 0.05, so we reject the null hypothesis. 
#We can conclude that at least one  of the Outdoor" groups has a different true mean for "Hwt" Cat's heart weight from the others. 
#i.e. there is a significant relationship between "Outdoor" and Hwt" Cat's heart weight

#Q1E
#Create a boxplot to see if there is any difference between the means of the response variable (selected in Q1B) across the two variables Sex and Outdoor. 
#Create a plot to see if there is an interaction effect for these two variables with the response variable. 
#Interpret these plots. Test to see if any of these factors and/or interactions have a significant relationship with the response variable (selected in Q1B) using an two-way ANOVA test. 
#(Note: categorical variables are not set up as factors, so change to factor prior to running analysis (as.factor())

boxplot(Hwt~Sex+Outdoor,data=cat_data, ylab="Heart Weight (g)")
print(boxplot(Hwt~Sex+Outdoor,data=cat_data, main="Comparing Heart Weight with Sex and Outdoor Status", xlab="Sex (1=Male, 2=Female) : Outdoor (1=Always, 2=Frequently, 3=Never)", ylab="Heart Weight (g)"))
#We can see that male cat's heart weight has a higher median and IQR for all outdoor categories
#We can see that the median heart weight tends to be lower for cats that are kept outdoors
#The spreads are similar across all boxplots and median roughly in middle
#There is an outlier present for a female cat that is never kept outside 
#This outlier appears to match the median heart weight of a male cat that was never kept outside so this possible could be an error.
#This outlier could have possibly been entered into the data as a female when it was actually a male.

#Interaction Plot 
interaction.plot(x.factor = Sex, 
                 trace.factor = Outdoor, 
                 response = Hwt,
                 fun = mean, trace.label ="Outdoor",
                 type = "b", legend = TRUE, 
                 main="Interaction Plot",xlab = "Sex", ylab="Heart Weight (g)",
                 pch=c(1,2),col=c(1,2))

#When we look at the interaction plot, the lines seem roughly parallel for cat's that are frequently and never kept outdoors
#The line for cat's that are always kept outdoors doesn't appear to be parallel with the other outdoor groups.
#Therefore it appears as if the interaction term may be important between these two variables

res2<-aov(Hwt~as.factor(Outdoor)*as.factor(Sex),data=cat_data)

par(mfrow=c(2,1))
plot(res2, which=1:2)

#Interpret Residuals vs Fitted 
#The spread of the groups appear to be roughly the same so we are happy with our equal variance assumption. 
#But there appears to be an outlier in the data represented by data point 144.
#The scatter seems randomly distributed within in the groups, so we are happy with our iid assumption.
#Interpret Normal Q-Q Plot
#The points fall roughly in a straight line indicating that the residuals are roughly normally distributed.
#There appears to be an outlier in the data represented by data point 144.

plot(res2, which=4:5)
#There is no need to worry about outliers as the cook's distance are all below 0.5

#Check if the data is balanced
addmargins(table(Outdoor,Sex))
#The data appears to be unbalanced as there is nearly twice as much males than females 
Anova(res2, type = "2") 

#Since it was unbalanced we used type 2 to adjust for different sizes of groups
#Testing 3 hypothesis here (2 main effects and 1 interaction effect)
#1st H0: The true mean for cat's heart weight is the same across outdoor groups
#1st H1: The true means for cat's heart weight are not the same across outdoor groups
#2nd H0: The true mean for cat's heart weight is the same across males and females
#2nd H1: The true means for cat's heart weight are not the same across males and females
#3rd H0: There is no interaction effect between sex and outdoor groups on cat's heart weight
#3rd H0: There is an interaction effect between sex and outdoor groups on cat's heart weight

#Interpretting the three tests:
#1  p-value < 2.2e-16 which is less than 0.05 so we reject the first null hypothesis
#We can conclude that there is a relationship between cat's heart weight and outdoor groups

#2  p-value = 7.374e-10 which is less than 0.05 so we reject the second null hypothesis
#We can conclude that there is a relationship between cat's heart weight and sex

#3  p-value = 0.0005062 which is less than 0.05 so we reject the third null hypothesis
#There is an interaction effect between sex and outdoor groups on cat's heart weight
#The relationship between cat's heart weight, outdoor groups and sex groups depend on each other

#Q2A
#Fit the “best” simple linear regression model, based on your answers from question 1, justify your choice. 
#Comment on whether the assumptions are satisfied, interpretation of the results and the fit of the model.

hist(Bwt)
#Cat's body weight does NOT appear to be symmetric.
hist(Height)
#Cat's height appears to be roughly symmetric and has a normal distribution.  
hist(Age)
#Cat's age appears to be roughly symmetric and has a normal distribution.

#Pearson
#H0: No Linear relationship between cat's heart weight and height, i.e. r_population correlation = 0
#H1: There is Linear relationship between cat's heart weight and height, i.e. r_population correlation != 0
cor.test(Hwt,Height)
#The p-value = 4.693e-07, There is Linear relationship between cat's heart weight and height.

#H0: No Linear relationship between cat's heart weight and age, i.e. r_population correlation = 0
#H1: There is Linear relationship between cat's heart weight and age, i.e. r_population correlation != 0
cor.test(Hwt,Age)
#The p-value = 6.072e-06, There is Linear relationship between cat's heart weight and age

#Kendall
#Using Kendall method as seen from the histogram of body weight it appears that it is not normally distributed.
#H0: No Linear relationship between cat's heart weight and body weight, i.e. r_population correlation = 0
#H1: There is Linear relationship between cat's heart weight and body weight, i.e. r_population correlation != 0
cor.test(Hwt, Bwt, method = "kendall")
#The p-value < 2.2e-16, There is Linear relationship between cat's heart weight and body weight.

#MODEL1 = cat's heart weight and cat's body weight
#Treating cat's heart weight as response variable
plot(Hwt~Bwt)

#lm(y~x)
slr_Hwt_Bwt = lm(Hwt~Bwt)
slr_Hwt_Bwt
#Check Assumptions
plot(slr_Hwt_Bwt,which=1)
#Random scatter => iid assumption
#Spread is consistent => equal variance assumption
plot(slr_Hwt_Bwt,which=2)
#Straight line therefore residuals are normal. Meets normality assumption.
plot(slr_Hwt_Bwt,which=4)
plot(slr_Hwt_Bwt,which=5)
#No outliers or high leverage points 

summary(slr_Hwt_Bwt)
#Cat's heart weight = -0.3567 + 4.0341 (Cat's body weight)
#For every 1kg increase of Cat's body weight we expect Cat's heart weight to increase on average by 4g
#The intercept can’t be sensibly interpreted as it means if we have 0kg for Cat's body weight then Cat's heart weight is -0.4g.
#H0: True Slope between Cat's heart weight and Cat's body weight = 0 i.e. No Relationship 
#H1: True Slope between Cat's heart weight and Cat's body weight != 0 i.e. Relationship 
#As the p-value < 2.2e-16 which is less than 0.05 we reject the null hypothesis.
#We can conclude that there is a relationship between Cat's heart weight and Cat's body weight
#Almost 65% of variation in Cat's heart weight is being explained by the model.
#Residual standard error = 1.452, therefore there is a small spread around the line.

plot(Hwt~Bwt)
abline(slr_Hwt_Bwt,col=2)

#MODEL2 = cat's heart weight and cat's height 
#Treating cat's heart weight as response variable
plot(Hwt~Height)

#lm(y~x)
slr_Hwt_Height = lm(Hwt~Height)
slr_Hwt_Height
#Check Assumptions
plot(slr_Hwt_Height,which=1)
#Random scatter => iid assumption
#Spread is consistent => equal variance assumption
plot(slr_Hwt_Height,which=2)
#Straight line but there is a slight tail at the end, but we are happy enough with it to say that it meets our normality assumption.
plot(slr_Hwt_Height,which=4)
plot(slr_Hwt_Height,which=5)
#No outliers or high leverage points 

summary(slr_Hwt_Height)
#Cat's heart weight = 5.97745 + 0.19134 (Cat's height)
#For every 1cm increase of Cat's height we expect Cat's heart weight to increase on average by 0.2g
#The intercept can’t be sensibly interpreted as it means if we have 0cm for Cat's height then Cat's heart weight is 6g.
#H0: True Slope between Cat's heart weight and Cat's height = 0 i.e. No Relationship 
#H1: True Slope between Cat's heart weight and Cat's height != 0 i.e. Relationship 
#As the p-value = 4.693e-07 which is less than 0.05 we reject the null hypothesis.
#We can conclude that there is a relationship between Cat's heart weight and Cat's height
#16% of variation in Cat's heart weight is being explained by the model.
#Residual standard error = 2.234, therefore there is a bigger spread around the line than we seen in Model1.

plot(Hwt~Height)
abline(slr_Hwt_Height,col=2)

#MODEL3 = cat's heart weight and cat's age
#Treating cat's heart weight as response variable
plot(Hwt~Age)

#lm(y~x)
slr_Hwt_Age = lm(Hwt~Age)
slr_Hwt_Age
#Check Assumptions
plot(slr_Hwt_Age,which=1)
#Random scatter => iid assumption
#Spread is consistent => equal variance assumption
plot(slr_Hwt_Age,which=2)
#Straight line therefore residuals are normal. Meets normality assumption.
plot(slr_Hwt_Age,which=4)
plot(slr_Hwt_Age,which=5)
#No outliers or high leverage points 

summary(slr_Hwt_Age)
#Cat's heart weight = 7.79764 + 0.23201 (Cat's age)
#For every 1 year increase of Cat's age we expect Cat's heart weight to increase on average by 0.2g
#The intercept can’t be sensibly interpreted as it means if we have 0 years for Cat's age then Cat's heart weight is 7.8g.
#H0: True Slope between Cat's heart weight and Cat's age = 0 i.e. No Relationship 
#H1: True Slope between Cat's heart weight and Cat's age != 0 i.e. Relationship 
#As the p-value = 6.072e-06 which is less than 0.05 we reject the null hypothesis.
#We can conclude that there is a relationship between Cat's heart weight and Cat's age
#13% of variation in Cat's heart weight is being explained by the model.
#Residual standard error = 2.273, therefore there is a bigger spread around the line than we seen in Model1 and it is very similar to Model2.

plot(Hwt~Age)
abline(slr_Hwt_Age,col=2)

#Q2B
#Fit a suitable multiple linear regression model, based on your answers from question 1, justify your choice. 
#Comment on whether the assumptions are satisfied, interpretation of the results and the fit of the model.

#Possible independent variables in the model
#Bwt, Height, Age as none of these variables are collinear with each other and have linear relationships with Hwt

lm1<-lm(Hwt~Bwt+Height+Age,data=cat_data)
#Open 2 visualizations in 1 graph
# 2 Rows and 1 Column
par(mfrow=c(2,1))
plot(lm1,which=1:2)
#Interpret Residuals vs Fitted 
#The spread appears to be roughly the same so we are happy with our equal variance assumption. 
#The scatter seems randomly distributed, so we are happy with our iid assumption.
#Interpret Normal Q-Q Plot
#The points fall in a straight line indicating that the residuals are normally distributed.
plot(lm1,which=4:5)
#No outliers or high leverage points
summary(lm1)
#Hwt = -1.76574 + 3.60366 (Bwt) + 0.06547 (Height) + 0.08102 (Age)
#For every 1kg in cat's body weight, cat's heart weight increases by 3.6g on average keeping cat's height and age constant
#For every 1cm in cat's height, cat's heart weight increases by 0.06547g on average keeping cat's body weight and age constant
#For every 1year in cat's age, cat's heart weight increases by 0.08102g on average keeping cat's body weight and height constant
#Cat's body weight has a more significant relationship with cat's heart weight than cat's height and age
#68% of variation in cat's heart weight is being explained by the model and it is better than our best simple linear regression model.
#Residual standard error = 1.39, therefore there is a small spread around the line and it is better than our best simple linear regression model.

#Q2C
#Are there any variables you would like to add/remove from the model and why? 
#Re-run the multiple linear regression model with these variables, if any, added/removed.  
#(You can do this more than once, in a step wise process if you think appropriate). 
#Compare the fit of this model to the model in Q2B.  
#Perform an F-test to compare the two models, stating your hypothesis and the conclusion of this test. 
#(Hint, use the anova() command). 
#If no variables are removed compare model ran in Q2A with model ran in Q2B.

#Adding Sex and Outdoor categories because they have a relationship with Cat's Heart Weight and may improve our first multiple linear regression model.
lm2<-lm(Hwt~Bwt+Height+Age+as.factor(Sex)+as.factor(Outdoor),data=cat_data)

par(mfrow=c(2,1))
plot(lm2,which=1:2)
#Interpret Residuals vs Fitted 
#The spread appears to be roughly the same so we are happy with our equal variance assumption. 
#The scatter seems randomly distributed, so we are happy with our iid assumption.
#Interpret Normal Q-Q Plot
#The points fall in a straight line indicating that the residuals are normally distributed.
plot(lm2,which=4:5)
#No outliers or high leverage points
#There is no need to run variance inflation factor  as we have proven relationships between all variables and cat's heart weight 
vif(lm2)
#As expected there is no values greater than 10 in the variance inflation factor results.
#Therefore we are not worried about collinearity

summary(lm2)
#Hwt = 0.34088 + 2.80654(Bwt) + 0.05686(Height) + 0.06803(Age) + -0.44678(Sex) + 0.39979(Outdoor) + 1.25937(Outdoor)
#For every 1kg in cat's body weight, cat's heart weight increases by 2.8g on average keeping all other variables constant
#For every 1cm in cat's height, cat's heart weight increases by 0.05686g on average keeping all other variables constant
#For every 1year in cat's age, cat's heart weight increases by 0.06803g on average keeping all other variables constant
#Going from male to female, cat's heart weight decreases by -0.44678g on average keeping all other variables constant
#Going from always outdoor to outdoor frequently, cat's heart weight increases by 0.39979g on average keeping all other variables constant
#Going from always outdoor to never outdoor, cat's heart weight increases by 1.25937g on average keeping all other variables constant


#70% of variation in cat's heart weight is being explained by the model and it is slightly better than our previous multiple linear regression model.
#Residual standard error = 1.352, therefore there is a small spread around the line and it is slightly better than our previous multiple linear regression model.
#We can also see that cat's body weight and going from from always outdoor to never outdoor is very significant
#Cat's height and age are also significant and it seems like going from male to female is less significant 

#F-test
#H0: The first multiple linear regression model without sex and outdoor is the preferred model  
#H1: The last multiple linear regression model with sex and outdoor is the preferred model  
anova(lm1,lm2)
#As the p-value = 0.01416 which is less than 0.05 we reject the null hypothesis.
#We can conclude that adding at least sex or outdoor to the model improves the fit of the model.
AIC(lm1,lm2)
BIC(lm1,lm2)
#AIC suggests that the last multiple linear regression model with sex and outdoor is the better model
#BIC suggests that the first multiple linear regression model without sex and outdoor is the better model

#Likelihood Ratio Test
install.packages("lmtest")
library(lmtest)
lrtest(lm2, lm1)
#H0: The second model and the first model fit the data equally well. i.e. first model is a better fit
#H1: The second model and the first model do not fit the data equally well. i.e. second model is a better fit
#As the p-value = 0.01416 which is less than 0.05 we reject the null hypothesis.
#This indicates that the second model and the first model do not fit the data equally well. As a result, we should employ the second model.

lrtest(lm2, lm1)

#Q2D
#Conclude your overall results. 
#Use your preferred model to predict the fitted values for your response variable
#Calculate the residual term if you are given the following data:
#Sex = 2, Bwt = 2.6, Hwt = 11.2, Height = 23.5, Age = 10, Outdoor = 2
#Sex = 1, Bwt = 3.1, Hwt = 12.5, Height = 25.4, Age = 13.2, Outdoor = 2

#Hwt = 0.34088 + 2.80654(Bwt) + 0.05686(Height) + 0.06803(Age) + -0.44678(Sex) + 0.39979(Outdoor) + 1.25937(Outdoor)
AnsHwt1 = 0.34088 + (2.80654 * 2.6) + (0.05686 * 23.5) + (0.06803 * 10) + (-0.44678 * 2) + (0.39979 * 2) + (1.25937 * 2)
AnsHwt1
#ANS: 12.07915
AHwt1 <- 11.2
#Difference between actual and predicted cat's heart weight
Diff1 <- AHwt1-AnsHwt1
Diff1
#ANS = -0.879154

AnsHwt2 = 0.34088 + (2.80654 * 3.1) + (0.05686 * 25.4) + (0.06803 * 13.2) + (-0.44678 * 1) + (0.39979 * 2) + (1.25937 * 2)
AnsHwt2
#ANS: 14.25493
AHwt2 <- 12.5
Diff2 <- AHwt2-AnsHwt2
Diff2
#ANS = -1.754934




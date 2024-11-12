# Cat Data<br>
Load in the data “Cat_Hwt.csv”.<br>
Data was collected on deceased male and female adult cats used for experiments: Sex: 1 for male and 2 for female.<br>
Bwt: body weight in kg.<br>
Hwt: heart weight in g.<br>
Height: Height in cm<br>
Age: Age in years.<br>
Outdoor: kept outdoors; 1 = Always, 2 = Frequently, 3 = Never<br><br>
## Relationships<br>
a) Investigate if any of the continuous numerical variables have a linear relationship by producing scatterplots, interpret these plots.<br>
b) What is the response variable and why? What is the research question ofinterest in this dataset?<br>
c) Investigate if any of the categorical /discrete variables seem to have a relationship with the response variable (selected in part b) using boxplots. Interpret these plots.<br>
d) Using partc),choose one categorical/discrete variable with at least 3 categories to test to see if there is any difference between the means of the response variable (selected in part b) using an one-way ANOVA test. (Note: if categorical variables are not set up as factors, so change to factor prior to running analysis (as.factor())<br>
e) Create a boxplot to see if there is any difference between the means of the response variable (selected in part b) across the two variables Sex and Outdoor. Also, create a plot to see if there is an interaction effect for these two variables with the response variable. Interpret these plots. Test to see if any of these factors and/or interactions have a significant relationship with the response variable (selected in part b) using an two-way ANOVA test. (Note: categorical variables are not set up as factors, so change to factor prior to running analysis (as.factor())<br>
## Multiple Linear Regression<br>
a) Fit the “best” simple linear regression model, based on your answers from question 1, justify your choice. Comment on whether the assumptions are satisfied, interpretation of the results and the fit of the model.
     1<br>
b) Fit a suitable multiple linear regression model,based on your answers from question 1, justify your choice. Comment on whether the assumptions are satisfied, interpretation of the results and the fit of the model.<br>
c) Are there any variables you would like to add/remove from the model and why? Re-run the multiple linear regression model with these variables, if any, added/removed. (You can do this more than once, in a stepwise process if you think appropriate). Compare the fit of this model to the model in part b. Perform an F-test to compare the two models, stating your hypothesis and the conclusion of this test. (Hint, use the anova() command). If no variables are removed compare model ran in part a with model ran in part b.<br>
d) Conclude your overall results. Use your preferred model to predict the fitted values for your response variable and calculate the residual term if you are given the following data:<br>
Sex Bwt Hwt Height Age Outdoor <br>2 2.6 11.2 23.5 10 2<br>
1 3.1 12.5 25.4 13.2 2

# Ionosphere Data:<br>
Load the mlbench package and load in the Ionosphere dataset from this package. Use help and str to understand the data that was collected on radar data at hospital in Labrador.
set.seed(153)<br><br>
(a) Have a look at the data, are the any variables to remove before running any models?<br>
(b) Createatrainingandtestdataset(70:30) <br><br>
## Decision Trees:<br>
(c) Create a decision tree for the train data with Class as the response and all of the other variables bar the variables discussed in (a) as predictors. Is this a classification tree or a regression tree?<br>
(d) Create a diagram of the decision tree created in (c). Interpret the tree diagram. Is this a useful visualisation for the data?<br>
(e) Using print function or otherwise, answer the following about the decision tree created in part (c):<br>
How many terminal nodes are there?<br>
What is the minimum number of observations in these terminal nodes?<br>
(f) Use the predict function to find the predicted values for the test dataset. Create a confusion matrix. What is the accuracy of this model?<br>
(g) Create a ROC plot to show the sensitivity vs specificity of the model. Find the Area Under the Curve (AUC). Interpret this.<br><br>
## Ensemble techniques Trees:<br>
(h) Use the bagging function on the training data to predict Class. What are the important variables used in this technique.<br>
(i) Use the predict function to find the predicted values for the test dataset using the model in (h). Create a confusion matrix. What is the accuracy of this model?<br>
(j) Use the random forest technique on the training data to predict Class.Whatarethe important variables used in this technique.<br>
(k) Use the plot function on your output of the randomForest function in (j). What does it tell you?<br>
(l) Predict the response for the test set and create the confusion matrix and calculate the accuracy. How does it compare with the model obtained for bagging in part (h)?<br>
(m) Perform boosting on the training data to predict Class. What are the important variables used in this technique.<br>
(n) Predict the response for the test set and create a confusion matrix. What is the accuracy of this model? How does it compare with the bagging and the random forest models?<br>
(o) Create a ROC plot to show the sensitivity vs specificity of all the models from the previous section. Find the Area Under the Curve (AUC) for each. Interpret this.
Support Vector Machines:<br>
(p) Have a look at the data briefly, do you think a linear or radial kernel is more appropriate given the visualizations.<br>
(q) Use tune.svm to select the best hyperparameter values for the svm model with the kernel selected in part q. Run this model on training dataset<br>
(r) Predict the response for the testset and create a confusion matrix. What is the accuracy of this model?<br><br>
## Overall:<br>
(s) Based on all the models performed here and the different measures of performance, which of these techniques would you recommend, giving reasons.
                        

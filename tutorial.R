Input =("
 Speaker  Likert
        Pooh      3
        Pooh      5
        Pooh      4
        Pooh      4
        Pooh      4
        Pooh      4
        Pooh      4
        Pooh      4
        Pooh      5
        Pooh      5
        Piglet    2
        Piglet    4
        Piglet    2
        Piglet    2
        Piglet    1
        Piglet    2
        Piglet    3
        Piglet    2
        Piglet    2
        Piglet    3
        Tigger    4
        Tigger    4
        Tigger    4
        Tigger    4
        Tigger    5
        Tigger    3
        Tigger    5
        Tigger    4
        Tigger    4
        Tigger    3
        ")

Data = read.table(textConnection(Input),header=TRUE)


### Order levels of the factor; otherwise R will alphabetize them

Data$Speaker = factor(Data$Speaker,
                      levels=unique(Data$Speaker))


### Create a new variable which is the likert scores as an ordered factor

Data$Likert.f = factor(Data$Likert,
                       ordered = TRUE)

library(lsmeans)
library(ordinal)
model = clm(Likert.f ~ Speaker,
            data = Data)

lsmeans(model,
        ~ Speaker, 
        adjust="tukey")

lsmeans(model,
        ~ Speaker, 
        adjust="tukey")       ### Tukey-adjusted comparisons


###  Check the data frame

library(psych)

headTail(Data)

str(Data)

summary(Data)


### Remove unnecessary objects

rm(Input)
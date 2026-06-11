# IrrTechInvest
R code for the submitted manuscript "Optimizing irrigation investment for water sustainability and food security under climate change"

This file contains multiple functions for panel data model analysis. It includes:

(1) regressionModel: this is to conduce panal data model analyisis

(2) estYld: this is to calculate yield changes with inputs of future climate and CO2

(3) incorpTechBenifit: this is to incoporate technological adjustment factors into model and estimate impacts

(4) simTechEffect: this is to simulate both irrigation technology upgrades and irrigatgion technology expansion 

(5) optComb: This is a function to optimize using a vector programming algorithm to speed-up calculation

    optIrrInvest: This is a function to call the optComb function with inputs of technology combinations, regional change in yield, irrigation and technology for certain crop 

In each function, I have added the comment for each input variable in the R file. 



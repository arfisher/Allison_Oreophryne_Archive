library(HH)
library(dplyr)

# create data
data = data.frame("Not_at_all_angry"=c(0.11,0.08,0.09,0.08,0.09,0.12,0.05,0.08),"Not_very_angry"=c(0.75,0.75,0.74,0.70,0.78,0.68,0.86,0.71),"Fairly_angry"=c(0.13,0.14,0.16,0.17,0.11,0.18,0.06,0.19),"Very_angry"=c(0.02,0.02,0.02,0.05,0.02,0.02,0.03,0.01),"Region"=c("North","Midlands","East","London","South","Wales","Scotland","Northern_Ireland"),"England" = c("England", "England", "England", "England", "England", "Not England", "Not England", "Not England"))

# make stacked bar chart
likert(Region ~.|England, layout=c(1,2), data, positive.order = TRUE, 
       scales=list(y=list(relation="free")),
       strip.left=strip.custom(bg="gray97"),
       strip=FALSE,
       as.percent = "noRightAxis", ReferenceZero = 2.5,
       main = 'Angry levels in different regions', 
       ylab = "Region", xlab = "Percentage",
       sub= list("Angry Level Rating",x=unit(.6, "npc")))
       
       
       
       


likert(Gender ~ ., data_2, ReferenceZero = 3, as.percent = "noRightAxis", main = "Satisfaction of Trump")       



my_data <- as_tibble(iris)

my_data %>% 
	summarize( 
		count= n(), 
		mean_sep = mean(Sepal.Length, na.rm=T), 
		mean_pet = mean(Petal.Length, na.rm=T)
		)

my_data %>% 
	group_by(Species) %>%
	summarize( 
		count= n(), 
		mean_sep = mean(Sepal.Length, na.rm=T), 
		mean_pet = mean(Petal.Length, na.rm=T)
		)

library(tidyverse)


NiVdata <- read_csv("data/niv_case_reports.csv")

#Plot cases by year
NiVdata %>%
  ggplot()+
  geom_col(aes(x=`start_year`, y = `Total Cases`))+
  xlab("Year outbreak started")

#Plot case fatality rate per year (perhaps better to do this per outbreak?)

NiVdata %>%
  ggplot()+
  geom_col(aes(x=`start_year`, y = `deceased`/`Total Cases`))+
  xlab("Year outbreak started")+
  ylab("Case fatality rate per year")

NiVdata$outbreakID <- paste("Outbreak", seq(1,nrow(NiVdata),1), NiVdata$start_year)


NiVdata %>%
  ggplot()+
  geom_histogram(aes(x = `deceased`/`Total Cases`))+
  xlab("Case fatality rate")+
  ylab("Number of outbreaks")

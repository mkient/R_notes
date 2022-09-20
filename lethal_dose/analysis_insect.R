
##### Data analysis paper #### 
#### Loading required library ####
library(tidyverse)
library(car)
library(emmeans)
library(ecotox)


### Loading data ####
wd = "E:/Papiers/Papier_Iro/data_R/"
setwd(wd)
df_base <- read.table("df_base.csv",header = T,sep = ";",dec = ".",
                       na.strings = NA)

names(df_base)

# mean(df_base$mort)
# min(df_base$mort)
# max(df_base$mort)
# var(df_base$mort)
# sd(df_base$mort)

## Normality test####
shapiro.test(df_base$mort)
qqnorm(df_base$mort,ylim=c(0,27), main = "Inactive")
qqline(df_base$mort)


#### Model of mortality rate #####
# data handling
df_base$populations <- as.factor(df_base$populations)
df_base$Traitements <- as.factor(df_base$Traitements)
df_base$mort_rate <- 100*(df_base$mort/df_base$total)
## without control data
df_base1 <- subset(df_base,
                   Traitements=="Bactivec"|
                     Traitements=="Griselesf"|
                     Traitements=="Bactivec+Griselesf")

df_base1$Traitements <- as.factor(df_base1$Traitements)
df_base1$populations <- as.factor(df_base1$populations)
df_base1$Replique <- as.factor(df_base1$Replique)
## Normality test
shapiro.test(df_base1$mort_rate)
qqnorm(df_base1$mort_rate,ylim=c(0,100), main = "Inactive")
qqline(df_base1$mort_rate)
#model
glm1 <- glm(cbind(mort, vivant)~Traitements*populations,
            data=df_base1, family = binomial("logit"))
summary(glm1)
car::Anova(glm1, type=3)
anova(glm1)
drop1(glm1,.~., test = "Chisq")
coef(glm1)

emmeans(glm1, pairwise~populations, type="response")
emmeans(glm1, pairwise~Traitements, type="response")
emmeans(glm1, pairwise~Traitements|populations, type="response")
emmeans(glm1, pairwise~populations|Traitements, type="response")

## Graph ~ geom_bar *** for mortality rate####
##data handling for graphic
tbl1 <- df_base%>%
  dplyr::group_by(populations,Traitements,dose)%>%
  dplyr::summarise(mort = mean(mort_rate),sd=sd(mort_rate),
                   se=sd(mort_rate)/sqrt(4))
tbl1[,4:6][is.na(tbl1[,4:6])] <- 0
## An_gambiae
gamb1 <- subset(tbl1, populations=="An_gambiae")
ag <- gamb1[1:18,]
min(ag$mort)
ag$dose <- as.factor(ag$dose)
ggp1 <- ggplot(data = ag,aes(x=dose,y=mort))+
  geom_bar(aes(fill = dose),stat="identity")+
  geom_errorbar(aes(ymin = mort - se, ymax = mort + se),
                width=0.2,)+
  facet_wrap(~Traitements)

ggp1 + 
  geom_hline(aes(yintercept = 90,linetype = "seuil de sensibilit? = 90%"),
             col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Reds")+
  labs(title="Taux de mortalit? d'Anopheles gambiae sl", 
       x="Dose (ppm)", y = "Taux de mortalit? (%)")

ggp1 + 
  geom_hline(aes(yintercept = 90, linetype = "Seuil de \nsensibilit?"),col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank(), legend.title = element_blank())+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Reds")+
  labs(x="Dose (ppm)", y = "Taux de mortalit? (%)")


## Culex
culex <- subset(tbl1, populations=="Culex_pipens")
cx1<- culex[1:18,]
cx1$dose <- as.factor(cx1$dose)
ggp2 <- ggplot(data = cx1,aes(x=dose,y=mort))+
  geom_bar(aes(fill = dose),stat="identity")+
  geom_errorbar(aes(ymin = mort - se, ymax = mort + se),
                width=0.2,)+
  facet_wrap(~Traitements)

ggp2 + 
  geom_hline(aes(yintercept = 90,linetype = "seuil de sensibilit? = 90%"),
             col="darkred")+
  scale_fill_brewer(palette = "Reds")+
  labs(title="Taux de mortalit? de Culex pipiens", 
       x="Dose (ppm)", y = "Taux de mortalit? (%)")

ggp2 + 
  geom_hline(aes(yintercept = 90,linetype = "seuil de sensibilit? = 90%"),
             col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Reds")+
  labs(title="Taux de mortalit? de Culex pipiens", 
       x="Dose (ppm)", y = "Taux de mortalit? (%)")

ggp2+
  geom_hline(aes(yintercept = 90,linetype = "seuil de \nsensibilit?"),
             col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank(), legend.title = element_blank())+
  scale_fill_brewer(palette = "Reds")+
  guides(fill = FALSE)+
  labs(x="Dose (ppm)", y = "Taux de mortalit? (%)")



### Lethal dose analyses #########
### LC_probit ~~ An gambiae #########
## Bactivec ~ An gambiae 
dbact1  <- subset(df_base1,populations=="An_gambiae"&
                          Traitements=="Bactivec")
bact1 <- LC_probit((mort/total)~log10(dose),weights = total,
                  p=c(50,95),data = dbact1  ,conf_level = 0.95)

##Griselesf ~ An gambiae 
dgris1  <- subset(df_base1,populations=="An_gambiae"&
                          Traitements=="Griselesf")
grisel <- LC_probit((mort/total)~log10(dose),weights = total,
                  p=c(50,95),data = dgris1 ,conf_level = 0.95)

## Bactivec+Griselesf ~ An gambiae 
dbactb_gr1 <- subset(df_base1,populations=="An_gambiae"&
                           Traitements=="Bactivec+Griselesf")
bact_gr1 <- LC_probit((mort/total)~log10(dose),weights = total,
                    p=c(50,95),data = dbactb_gr1 ,conf_level = 0.95)

## LC plot ~ regression curve*** An gambiae
gamb2 <- subset(df_base1,populations=="An_gambiae")

p0 <- ggplot(data = gamb2,
             aes(x = log10(dose), y =(mort/total))) +
  geom_point()+ 
  geom_smooth(method = "glm",
              method.args = list(family = binomial(link = "probit")),
              aes(weight = total), colour = "#0000FF", se = TRUE)+
  facet_wrap(~Traitements)

######## comments
# p0 + scale_fill_brewer(palette = "Reds")+
#   labs(title="Courbe de regression des mortalites d'An gambiae en fonction de la dose",
#        x="Log10(dose) ", y = "Taux de mortalit? ")

# p1 <- ggplot(data = bg,
#              aes(x = log10(dose), y =(mort/total))) +
#   geom_point()+
#   geom_smooth(method = "glm",
#               method.args = list(family = binomial(link = "probit")),
#               aes(weight = total), colour = "#0000FF", se = TRUE)
#
# p1+ geom_segment(aes(x = log10(bact$dose[1]), y = 0,xend = log10(bact$dose[1]),
#                      yend = 0.50, colour = "segment"))+
#   geom_segment(aes(x = -0.5, y = 0.50, xend = log10(bact$dose[1]),
#                    yend = 0.50, colour = "segment"))+
#   geom_segment(aes(x = log10(bact$dose[2]), y = 0,xend = log10(bact$dose[2]),
#                    yend = 0.950, colour = "segment"))+
#   geom_segment(aes(x = -0.5, y = 0.950, xend = log10(bact$dose[2]),
#                    yend = 0.950, colour = "segment"))

### LC_probit ~~~ Culex pipiens ####
##Bactivec ~ culex 
dbact2  <- subset(df_base1,populations=="Culex_pipens"&
                Traitements=="Bactivec")
bact2 <- LC_probit((mort/total)~log10(dose),weights = total,
                   p=c(50,95),data = dbact2 ,conf_level = 0.95)

##Griselesf ~ Culex_pipens 
dgris2  <- subset(df_base1,populations=="Culex_pipens"&
                Traitements=="Griselesf")
grise2 <- LC_probit((mort/total)~log10(dose),weights = total,
                    p=c(50,95),data = dgris2 ,conf_level = 0.95)

## Bactivec+Griselesf ~ Culex_pipens 
dbactb_gr2 <- subset(df_base1,populations=="Culex_pipens"&
               Traitements=="Bactivec+Griselesf")
bact_gr2 <- LC_probit((mort/total)~log10(dose),weights = total,
                     p=c(50,95),data = dbactb_gr2 ,conf_level = 0.95)

## LC plot ~ regression curve *** Culex pipiens
culex_p<- subset(df_base1,populations=="Culex_pipens")

p1 <- ggplot(data = culex_p,
             aes(x = log10(dose), y =(mort/total))) +
  geom_point()+ 
  geom_smooth(method = "glm",
              method.args = list(family = binomial(link = "probit")),
              aes(weight = total), colour = "#0000FF", se = TRUE)+
  facet_wrap(~Traitements)

### Data table LC_probit#####
LCg <- rbind(bact1,grisel,bact_gr1)
LCc <- rbind(bact2,grise2,bact_gr2)
Larvicide <- c("Bactivec","Bactivec","Griselesf","Griselesf",
               "Bactivec+Griselesf", "Bactivec+Griselesf",
               "Bactivec","Bactivec","Griselesf","Griselesf",
               "Bactivec+Griselesf", "Bactivec+Griselesf")
pop <- rep(c("An_gambiae","Culex_pipens"),c(6,6))
LC <- rbind(LCg,LCc)
LC <- LC%>%add_column(data.frame(pop,Larvicide),
                      .before = "p")
# Save as excel table ## 
write.table(LC, file = "Lethal_dose.csv",qmethod = "double",
            row.names = FALSE,sep=";")

###### end of the analyses ######


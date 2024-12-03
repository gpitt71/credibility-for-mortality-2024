library(data.table)
library(xtable)
dt1 <- fread("C:/Users/gpitt/Documents/Postdoc/Torino/Mortality/results/gender_validation_testing_selection_l2.csv")
dt1[["group"]] <- "G2"
dt2 <- fread("C:/Users/gpitt/Documents/Postdoc/Torino/Mortality/results/scandinavian_validation_testing_selection_l2.csv")
dt2[["group"]] <- "G1"
dt3 <- fread("C:/Users/gpitt/Documents/Postdoc/Torino/Mortality/results/ita_validation_testing_selection_l2.csv")
dt3[["group"]] <- "G3"


dt <- rbind(dt2,
            dt1,
            dt3)

tmp <- dt[order(forecasting_horizon),
          -c('Metric.x',
             'Metric.y')][,c("forecasting_horizon",
                             "group",
                             "model",
                             "family",
                             "mse_val",
                             "Value")]

tmp[,model:=toupper(model)][,mse_val := mse_val *10e03][,Value := Value *10e03]


tmp[, family := fcase(
  family == 0, "S1",
  family == 1, "S2",
  family == 2, "S3"
)]

colnames(tmp) <- c("Forecasting Horizon",
                   "Group",
                   "Selected Model",
                   "Family",
                   "MAE (val)",
                   "MAE (test)")



print(xtable(tmp, digits=4),include.rownames = FALSE)



# Models rankings ----

library(tidyr)

dt1 <- fread("C:/Users/gpitt/Documents/Postdoc/Torino/Mortality/results/gender_validation_l2.csv")
dt1[["group"]] <- "G2"
dt2 <- fread("C:/Users/gpitt/Documents/Postdoc/Torino/Mortality/results/scandinavian_validation_l2.csv")
dt2[["group"]] <- "G1"
dt3 <- fread("C:/Users/gpitt/Documents/Postdoc/Torino/Mortality/results/ita_validation_l2.csv")
dt3[["group"]] <- "G3"


dt <- rbind(dt2,
            dt1,
            dt3)


tmp <- dt %>%
  pivot_longer(cols=c("mse_0",
                      "mse_1",
                      "mse_2"),
               names_to = c("Strategy"),
               values_to = c("Value")) %>%
  mutate(Strategy = as.numeric(gsub("[^0-9]", "", Strategy))) %>%
  as.data.table()


tmp[, Strategy := fcase(
  Strategy == 0, "S1",
  Strategy == 1, "S2",
  Strategy == 2, "S3"
)][, model := toupper(model)][,tracer := paste0(model,"-",Strategy)]

tmp[,Rank := rank(Value), by = .(forecasting_horizon, group)]


library(ggplot2)


tmp[,yax := paste0(group, "-", as.character(forecasting_horizon))]

tmp <- tmp %>% mutate(Rank=as.factor(Rank),
                      yax=factor(yax, levels = c("G1-1",
                                                 "G2-1",
                                                 "G3-1",
                                                 "G1-5",
                                                 "G2-5",
                                                 "G3-5",
                                                 "G1-12",
                                                 "G2-12",
                                                 "G3-12")),
                      tracer=factor(tracer, levels = c("APC-S1",
                                                       "LC-S1",
                                                       "RH-S1",
                                                       "APC-S2",
                                                       "LC-S2",
                                                       "RH-S2",
                                                       "APC-S3",
                                                       "LC-S3",
                                                       "RH-S3")))

ggplot(tmp %>% arrange(forecasting_horizon, as.numeric(gsub("[^0-9]", "", group))),
       aes(x = tracer, y = yax, fill = Rank)) +
  geom_tile(aes(fill = Rank), colour = "black") +
  ggtitle(" ") +
  theme_classic()+
  geom_text(aes(label = Rank))+
  scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","#4169E1"))(9), na.value="#EEEEEE", name="Rank") +
  xlab(" ") +
  ylab(" ")+
  theme(axis.text.y = element_text(size=15),
        axis.text.x  = element_text(size=15,
                                    angle = 90, hjust = 1),
        legend.text=element_text(size=20))

ggsave(paste0("C:\\Users\\gpitt\\Pictures\\mortality_project\\rankings_mortality.pdf"),
       width = 10,
       height= 5)






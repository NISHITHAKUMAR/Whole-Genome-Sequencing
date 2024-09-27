

library(dplyr)
library(ggplot2)

path = '/Users/nishitha/Dropbox/My Mac (Nishitha MacBook Pro)/Desktop/NEC/Metrics_VIOLIN"

setwd(dirname(path))

dir.create("Violin_Plots")
setwd("Violin_Plots")

# AUC
png("AUC_Violin.png", height = 10, width = 4, units = "in", res = 300)
ggplot(data = data, aes(x=Antibiotic, y = AUC_Mean, fill= Antibiotic))+
  theme(axis.text.y = element_text(face="bold", size = 14),
        axis.title.y = element_text(face="bold", size = 14),
        axis.title.x = element_text(face="bold", size = 14))+
  labs(y="AUC")+
  geom_violin()+
  geom_boxplot(width=0.1)
  coord_flip()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()


# Accuracy
png("Accuracy_Violin.png", height = 10, width = 4, units = "in", res = 400)
ggplot(data = data, aes(x=Antibiotic, y = Acc_Mean, fill= Antibiotic))+
  theme(axis.text.y = element_text(face="bold", size = 14),
        axis.title.y = element_text(face="bold", size = 14),
        axis.title.x = element_text(face="bold", size = 14))+
  labs(y="Accuracy")+
  geom_violin()+
  geom_boxplot(width=0.1)+
  coord_flip()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

dev.off()

# Sensitivity
png("Sensitivity_Violin.png", height = 10, width = 4, units = "in", res = 400)
ggplot(data = data, aes(x=Antibiotic, y = Sens_Mean, fill= Antibiotic))+
  theme(axis.text.y = element_text(face="bold", size = 14),
        axis.title.y = element_text(face="bold", size = 14),
        axis.title.x = element_text(face="bold", size = 14))+
  labs(y="Sensitivity")+
  geom_violin()+
  geom_boxplot(width=0.1)+
  coord_flip()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

# Specificity
png("Specificity_Violin.png", height = 10, width = 4, units = "in", res = 400)
ggplot(data = data, aes(x=Antibiotic, y = Spec_Mean, fill= Antibiotic))+
  theme(axis.text.y = element_text(face="bold", size = 14),
        axis.title.y = element_text(face="bold", size = 14),
        axis.title.x = element_text(face="bold", size = 14))+
  labs(y="Specificity")+
  geom_violin()+
  geom_boxplot(width=0.1)+
  coord_flip()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

# Precision
png("Precision_Violin.png", height = 10, width = 4, units = "in", res = 400)
ggplot(data = data, aes(x=Antibiotic, y = Prec_Mean, fill= Antibiotic))+
  theme(axis.text.y = element_text(face="bold", size = 14),
        axis.title.y = element_text(face="bold", size = 14),
        axis.title.x = element_text(face="bold", size = 14))+
  labs(y="Precision")+
  geom_violin()+
  geom_boxplot(width=0.1)+
  coord_flip()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

# Kappa
png("Cohens_Kappa_Violin.png", height = 10, width = 4, units = "in", res = 400)
ggplot(data = data, aes(x=Antibiotic, y = Kappa_Mean, fill= Antibiotic))+
  theme(axis.text.y = element_text(face="bold", size = 14),
        axis.title.y = element_text(face="bold", size = 14),
        axis.title.x = element_text(face="bold", size = 14))+
  labs(y="Cohens Kappa")+
  geom_violin()+
  geom_boxplot(width=0.1)+
  coord_flip()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()
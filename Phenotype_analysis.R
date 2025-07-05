library(ggplot2)
library(tidyverse)
library(lme4)
library(car)
library(plotrix)

setwd("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/phenotype_analysis/")

phen_32 <- read_csv("phenotypes_rna_seq.txt") |> 
  mutate(pop = str_replace(pop, "P4", "P18")) |>
  mutate(pop = str_replace(pop, "P2", "P4")) |>
  mutate(pop = str_replace(pop, "P3", "P26")) |>
  mutate(pop = str_replace(pop, "P5", "P14")) |>
  mutate(popfam = paste(pop,fam, sep = "_")) |>
  mutate_at(c('Temp', 'pop', 'popfam'), as.factor)

phen_42 <- read_csv("phenotypes_metamorphosis.txt") |> 
  mutate(pop = str_replace(pop, "P4", "P18")) |>
  mutate(pop = str_replace(pop, "P2", "P4")) |>
  mutate(pop = str_replace(pop, "P3", "P26")) |>
  mutate(pop = str_replace(pop, "P5", "P14")) |> 
  mutate(popfam = paste(pop,fam, sep = "_")) |>
  mutate_at(c('Temp', 'pop', 'popfam'), as.factor)


##### RNAseq individuals gosner stage 32

## Liver
phen_32_temp_15_liver <- phen_32 |> 
  group_by(fam, pop, Temp) |> 
  filter(Temp == 15 & Tissue == "liver") |> 
  ungroup() 
  
phen_32_temp_20_liver <- phen_32 |> 
  group_by(fam, pop, Temp) |> 
  filter(Temp == 20 & Tissue == "liver") |> 
  mutate(LP_mean_fam = mean(LP), MAM_mean_fam = mean(mass)) |> 
  ungroup() 

phen_32_15_joined_liver <- left_join(phen_32_temp_15_liver,phen_32_temp_20_liver, by="popfam") |> 
  na.omit() 

phen_32_15_joined_liver_plast <- phen_32_15_joined_liver |> 
  mutate(LP_plasticity = LP.x - LP_mean_fam, MAM_plasticity =  mass.x - MAM_mean_fam) |> 
  write_delim("phen_32_liver_plasticity.txt", delim="\t")

phen_32_15_joined_liver_final <- phen_32_15_joined_liver_plast |> 
  distinct(ID.x, .keep_all = TRUE)
  
print(phen_32_15_joined_liver_final %>% 
        group_by(ID.x) %>%
        summarize(n()), n=50)

## Muscle
phen_32_temp_15_muscle <- phen_32 |> 
  group_by(fam, pop, Temp) |> 
  filter(Temp == 15 & Tissue == "muscle") |> 
  mutate(popfam = paste(pop,fam, sep = "_")) |> 
  ungroup() 

phen_32_temp_20_muscle <- phen_32 |> 
  group_by(fam, pop, Temp) |> 
  filter(Temp == 20 & Tissue == "muscle") |> 
  mutate(LP_mean_fam = mean(LP), MAM_mean_fam = mean(mass)) |> 
  mutate(popfam = paste(pop,fam, sep = "_")) |> 
  ungroup() 

phen_32_15_joined_muscle <- left_join(phen_32_temp_15_muscle,phen_32_temp_20_muscle, by="popfam") |> 
  na.omit() 

phen_32_15_joined_muscle_plast <- phen_32_15_joined_muscle |> 
  mutate(LP_plasticity = LP.x - LP_mean_fam, MAM_plasticity =  mass.x - MAM_mean_fam) |> 
  write_delim("phen_32_muscle_plasticity.txt", delim="\t")

phen_32_15_joined_muscle_plast_final <- phen_32_15_joined_muscle_plast |> 
  distinct(ID.x, .keep_all = TRUE) |> 
  mutate_at(c('Temp.x', 'pop.x', 'popfam'), as.factor)


### Testing significance of pop and temp using linear mixed effect models with family nested in population as a random effect.
LP_32 <-lmer(LP ~ pop*Temp + (1|popfam), phen_32)
summary(LP_32)
Anova(LP_32)

LP_32_plasticity <-lmer(LP_plasticity ~ pop.x + (1|popfam), phen_32_15_joined_muscle_plast_final)
summary(LP_32_plasticity)
Anova(LP_32_plasticity)

MAM_32 <-lmer(mass ~ pop*Temp + (1|popfam), phen_32)
summary(MAM_32)
Anova(MAM_32)

MAM_32_plasticity <-lmer(MAM_plasticity ~ pop.x + (1|popfam), phen_32_15_joined_muscle_plast_final)
summary(MAM_32_plasticity)
Anova(MAM_32_plasticity)


### Plotting phenotypes for stage 32

phen_32_to_plot <- phen_32 |> 
  group_by(pop, Temp) |> 
  summarize(LP_mean=mean(LP), LP_err = std.error(LP), LP_sd = sd(LP), MAM_mean = mean(mass), MAM_err = std.error(mass), MAM_sd = sd(mass))

plast_32_to_plot  <-  phen_32_15_joined_muscle_plast_final |> 
  group_by(pop.x, Temp.x) |>
  summarize(LP_plasticity_mean=mean(LP_plasticity), LP_plasticity_err = std.error(LP_plasticity), LP_plasticity_sd = sd(LP_plasticity), MAM_plasticity_mean = mean(MAM_plasticity), MAM_plasticity_err = std.error(MAM_plasticity), MAM_plasticity_sd = sd(MAM_plasticity))


ggplot(phen_32_to_plot, aes(x = pop, y = LP_mean, colour = Temp)) + 
  geom_errorbar(aes(ymin = LP_mean - LP_err, ymax = LP_mean + LP_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Larval period (d)") +
  theme_bw()
ggsave("LP_32.pdf", dpi=600)

ggplot(phen_32_to_plot, aes(x = pop, y = MAM_mean, colour = Temp)) + 
  geom_errorbar(aes(ymin = MAM_mean - MAM_err, ymax = MAM_mean + MAM_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Mass at stage 32 (g)") +
  theme_bw()
ggsave("MAM_32.pdf", dpi=600)

ggplot(plast_32_to_plot, aes(x = pop.x, y = LP_plasticity_mean, colour = Temp.x)) + 
  geom_errorbar(aes(ymin = LP_plasticity_mean - LP_plasticity_err, ymax = LP_plasticity_mean + LP_plasticity_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Larval period plasticity (d)") +
  theme_bw() + theme(legend.position = "none")
ggsave("LP_plasticity_32.pdf", dpi=600)

ggplot(plast_32_to_plot, aes(x = pop.x, y = MAM_plasticity_mean, colour = Temp.x)) + 
  geom_errorbar(aes(ymin = MAM_plasticity_mean - MAM_plasticity_err, ymax = MAM_plasticity_mean + MAM_plasticity_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Mass at stage 32 plasticity (g)") +
  theme_bw() + theme(legend.position = "none")
ggsave("MAM_plasticity_32.pdf", dpi=600)


##### Common garden experiment until metamorphosis
phen_42_temp_15 <- phen_42 |> 
  group_by(fam, pop, Temp) |> 
  filter(Temp == 15) |> 
  ungroup() 

phen_42_temp_20 <- phen_42 |> 
  group_by(fam, pop, Temp) |> 
  filter(Temp == 20) |> 
  mutate(LP_mean_fam = mean(LP), MAM_mean_fam = mean(mass)) |> 
  ungroup() 

phen_42_15_joined_ <- left_join(phen_42_temp_15,phen_42_temp_20, by="popfam") |> 
  na.omit() 

phen_42_15_joined_plast <- phen_42_15_joined_ |> 
  mutate(ID.x = paste(pop.x,fam.x,ind.x, sep = "_")) |> 
  mutate(LP_plasticity = LP.x - LP_mean_fam, MAM_plasticity =  mass.x - MAM_mean_fam)

phen_42_15_joined_final <- phen_42_15_joined_plast |> 
  distinct(ID.x, .keep_all = TRUE)

print(phen_42_15_joined_final %>% 
        group_by(ID.x) %>%
        summarize(n()), n=200)


### Testing significance of pop and temp using linear mixed effect models with family nested in population as a random effect.
LP_42 <-lmer(LP ~ pop*Temp + (1|popfam), phen_42)
summary(LP_42)
Anova(LP_42)

LP_42_plasticity <-lmer(LP_plasticity ~ pop.x + (1|popfam), phen_42_15_joined_final)
summary(LP_42_plasticity)
Anova(LP_42_plasticity)

MAM_42 <-lmer(mass ~ pop*Temp + (1|popfam), phen_42)
summary(MAM_42)
Anova(MAM_42)

MAM_42_plasticity <-lmer(MAM_plasticity ~ pop.x + (1|popfam), phen_42_15_joined_final)
summary(MAM_42_plasticity)
Anova(MAM_42_plasticity)


### Plotting phenotypes for stage 42
phen_42_to_plot <- phen_42 |> 
  group_by(pop, Temp) |> 
  summarize(LP_mean=mean(LP), LP_err = std.error(LP), LP_sd = sd(LP), MAM_mean = mean(mass), MAM_err = std.error(mass), MAM_sd = sd(mass))

plast_42_to_plot  <-  phen_42_15_joined_final |> 
  group_by(pop.x, Temp.x) |>
  summarize(LP_plasticity_mean=mean(LP_plasticity), LP_plasticity_err = std.error(LP_plasticity), LP_plasticity_sd = sd(LP_plasticity), MAM_plasticity_mean = mean(MAM_plasticity), MAM_plasticity_err = std.error(MAM_plasticity), MAM_plasticity_sd = sd(MAM_plasticity))


ggplot(phen_42_to_plot, aes(x = pop, y = LP_mean, colour = Temp)) + 
  geom_errorbar(aes(ymin = LP_mean - LP_err, ymax = LP_mean + LP_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Larval period (d)") +
  theme_bw()
ggsave("LP_42.pdf", dpi=600)

ggplot(phen_42_to_plot, aes(x = pop, y = MAM_mean, colour = Temp)) + 
  geom_errorbar(aes(ymin = MAM_mean - MAM_err, ymax = MAM_mean + MAM_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Mass at stage 42 (g)") +
  theme_bw()
ggsave("MAM_42.pdf", dpi=600)

ggplot(plast_42_to_plot, aes(x = pop.x, y = LP_plasticity_mean, colour = Temp.x)) + 
  geom_errorbar(aes(ymin = LP_plasticity_mean - LP_plasticity_err, ymax = LP_plasticity_mean + LP_plasticity_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Larval period plasticity (d)") +
  theme_bw() + theme(legend.position = "none")
ggsave("LP_plasticity_42.pdf", dpi=600)

ggplot(plast_42_to_plot, aes(x = pop.x, y = MAM_plasticity_mean, colour = Temp.x)) + 
  geom_errorbar(aes(ymin = MAM_plasticity_mean - MAM_plasticity_err, ymax = MAM_plasticity_mean + MAM_plasticity_err), width = 0.15, linewidth  = 0.75) +
  geom_point(shape = 16, size  = 2) +
  xlab("Population") + ylab("Mass at stage 42 plasticity (g)") +
  theme_bw() + theme(legend.position = "none")
ggsave("MAM_plasticity_42.pdf", dpi=600)




LP_42 <- left_join((phen_42 |> 
  group_by(pop) |> 
  filter(Temp == 15) |> 
  summarize(LP_15 = mean(LP))), 
  (phen_42 |> 
     group_by(pop) |> 
     filter(Temp == 20) |> 
     summarize(LP_20 = mean(LP))), by = "pop") |> 
  mutate(perc_LP = 1 - LP_20/LP_15, days_change_LP = LP_15 - LP_20)


MAM_42 <- left_join((phen_42 |> 
             group_by(pop) |> 
             filter(Temp == 15) |> 
             summarize(MAM_15 = mean(mass))), 
          (phen_42 |> 
             group_by(pop) |> 
             filter(Temp == 20) |> 
             summarize(MAM_20 = mean(mass))), by = "pop") |> 
  mutate(perc_mass = (1 - MAM_20/MAM_15)*100, weight_change_MAM = MAM_15 - MAM_20)


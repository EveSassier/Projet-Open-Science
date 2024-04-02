## Code Alice GALLOT et Eve SASSIER, Projet Open Science

rm(list=ls())

install.packages(c("ape", "phytools", "RPANDA", "egg", "car", "maps"))

library("devtools")
devtools::install_github("jinyizju/V.PhyloMaker2")


library(ape)
library(phytools)
library(tidyverse)
library(RPANDA)
library(V.PhyloMaker2)
library(ggplot2)
library(egg)
library(car)
library(maps)


#### Etape 1 : Tri des données  ####

##on extrait les donnees du fichier csv sous forme de dataframe
data <- read.csv("rsbl20210478_si_002.csv")
data_phylo_plant <- read.csv("data_phylo_plant.csv", sep=";")
tree_phylo_plant <- read.tree(file="phylo_plant.tre", keep.multi = FALSE)
tree_cytb <- read.nexus("cytb.nex")
localisation <- read.csv("localisation.csv", sep= ";")
auteurs <- read.csv(file="auteurs.csv", sep=";")

## on extrait deux vecteurs avec respectivement toutes les especes de chauves souris
## et tous les genre de plantes 
plant_genus <- unique(data$Plan_genus)
plant_genus <- as.factor(plant_genus)
bat_species <- unique(data$Bat_specie)
bat_species <- as.factor(bat_species)

## on cree une matrice nulle dont les colonnes correspondent aux genres de plantes et
## les lignes correspondent aux esp?ces de chauve-souris afin de repr?senter les interactions
## entre chauve-souris et plantes (relation binaire)
matrice_interactions <- matrix(0, nrow=length(bat_species), ncol=length(plant_genus))
colnames(matrice_interactions) <- plant_genus
rownames(matrice_interactions) <- bat_species


##on cree la phylogenie des genres de plantes avec V.PhyloMaker2
#plot(tree_phylo_plant)

#On supprime la 2e branche des 2 Piper, Passiflora, et Manilkara (car 2 espèces pour 1 seul genre)
#On supprime aussi la branche "X", qui n'apporte aucune information
tree_phylo_plant <- drop.tip(tree_phylo_plant, c("Piper_umbellatum", "Passiflora_acreana", "Manilkara_zapota", "X"), 
                             trim.internal = TRUE, subtree = FALSE, root.edge = 1, 
                             rooted = is.rooted(tree_phylo_plant), collapse.singles = FALSE, interactive = FALSE)


# On supprime le 2nd membre des noms d'espece, pour n'avoir que le nom du genre
tree_phylo_plant$tip.label<- sapply(strsplit(split="_", tree_phylo_plant$tip.label), "[[", 1)

#on vérifie que la taille des donnees de l'arbre phylogenetique et la matrice d'interactions correspond
#setdiff(tree_phylo_plant$tip.label, plant_genus)

##on renomme les genres des especes pour qu'ils correspondent aux genres de la matrice d'interactions
#après avoir fait l'abre
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Gesneria"] <- "Gesneriaceae"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Trianaea"] <- "Trianea"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Chrysophyllum"] <- "Crysophyllum"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Brosimum"] <- "Brosinum"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Cucurbita"] <- "Cucurbitaceae"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Tovomita"] <- "Tomovita"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Melastoma"] <- "Melastomataceae"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Ottonia"] <- "Otonia"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Manilkara"] <- "Achras"
tree_phylo_plant$tip.label[tree_phylo_plant$tip.label=="Piper"] <- "Pothomorphe"
#setdiff(tree_phylo_plant$tip.label, plant_genus)
#setdiff(plant_genus, tree_phylo_plant$tip.label)

### On crée l'arbre sur la phylogénie des chauves-souris (on pourra aussi par la suite se contenter d'importer tree_cytb_modif)

#tree_cytb[["tip.label"]] # affiche toutes les espèces de l'arbre tree_cytb (la commande "tree_cytb$tip.label" marche aussi)
#rownames(matrice_interactions) #affiche les noms de toutes les espèces de chauves-souris qui sont dans la matrice d'interactions
#rownames(matrice_interactions) %in% tree_cytb$tip.label #pour vérifier que toutes nos espèces de chauves-souris sont bien dans l'arbre choisi
# -> ce n'est pas le cas ! 

# On teste donc avec tous les autres arbres dispos, pour voir lequel contient le plus d'espèces que nous avons dans notre matrice d'interactions
# sum((rownames(matrice_interactions) %in% tree_cytb$tip.label)==FALSE) # on fait ça pour chaque arbre
# c'est bien l'arbre cytb qui contient le plus d'espèces qu'on souhaite, mais il y en a qd même 13 qui ne sont pas répertoriées dans cet arbre.
# Les espèces absentes de l'arbre mais présentes dans la matrice d'interactions sont :
# Hsunycteris_thomasi, Glossophaga_sp1, Phyllostomus_elongatus, Brachyphylla_nana_pumila, Sturnira_sp1, Dermanura_sp1, Vampyriscus_brocki,
# Hylonycteris_underwoodi, Artibeus_sp, Platyrrhinus_sp, Artibeus_glaucus, Phyllostomus_discolor, Anoura_cultrata

## Création d'un nouvel arbre, ne contenant que les espèces qui nous intéressent
## (càd présentes dans la matrice d'interactions):
tree_cytb_modif <- drop.tip(tree_cytb, setdiff(tree_cytb$tip.label, rownames(matrice_interactions)), 
                            trim.internal = TRUE, subtree = FALSE, root.edge = 1,
                            rooted = is.rooted(tree_cytb), collapse.singles = FALSE,
                            interactive = FALSE)
#plot(tree_cytb_modif)


#### Etape 2 : A l'echelle regionale ####

## on utilise un subdata_wide en fonction values to Latitude 
subdata <- select(data, Bat_specie, Plan_genus, Latitude)
subdata_wide <- pivot_wider(subdata, names_from = Plan_genus, values_from = Latitude)

## Toute entree NULL dans le subdata_wide signifie qu'il n'y a pas d'interactions entre
## le genre de plante (colonne) et l'espece de chauve-souris (wide)
## on affecte dans la matrice d'interaction 0 si il n'y a pas d'interaction 1 sinon
for (i in 2:ncol(subdata_wide)){
  for (j in 1:nrow(subdata_wide)){
    if (subdata_wide[j, i]=="NULL"){matrice_interactions[j, i-1]<-0}
    else{matrice_interactions[j, i-1]<- 1}
  }
}


#On enlève la colonne X de la matrice d'interactions
matrice_interactions <- matrice_interactions[,-8]

setdiff(colnames(matrice_interactions), tree_phylo_plant$tip.label)
matrice_interactions <- as.data.frame(matrice_interactions)

#On somme les colonnes Achras et Manilkara et on place le résultat dans la colonne Achras
#(pas de 2 dans le résultat donc pas besoin de les transformer en 1)
matrice_interactions <- mutate(matrice_interactions, Achras = Achras + Manilkara)

#On enlève la colonne "Manilkara" de la matrice d'interactions 
matrice_interactions <- matrice_interactions[,-49]

#On somme les colonnes Piper et Pothomorphe et on place le résultat dans la colonne Pothomorphe
matrice_interactions <- mutate(matrice_interactions, Pothomorphe = Piper + Pothomorphe)

#On transforme les 2 en 1 dans la colonne sommée "Pothomorphe"
o <- nrow(matrice_interactions)
for (i in 1:o){
  if (matrice_interactions$Pothomorphe[i]==2){
    matrice_interactions$Pothomorphe[i]=1
  }
  else if (matrice_interactions$Pothomorphe[i]==1){
    matrice_interactions$Pothomorphe[i]=1
  }
  else {matrice_interactions$Pothomorphe[i]=0}
}

#On enlève la colonne "Piper" de la matrice d'interactions
matrice_interactions <- matrice_interactions[,-3]

# on supprime les colonnes de la matrice qui correspondent pour lesquelles on n'a pas l'information
# phylogénétique, càd les noms de genre qui correspondent à des noms de famille
genus_interet <- intersect(colnames(matrice_interactions), tree_phylo_plant$tip.label)
matrice_interactions <- matrice_interactions %>% select(all_of(genus_interet))


## Plantes
#On realise le test de Mantel plantes a echelle regionale
matrice_inter_phylo_plant = matrice_interactions
test_Mantel_plant <- phylosignal_network(network = matrice_inter_phylo_plant, tree_A= tree_phylo_plant, method = "Jaccard_binary", correlation = "Pearson", nperm = 10000)
test_Mantel_plant <- as.data.frame(test_Mantel_plant)
test_Mantel_plant <- t(test_Mantel_plant)

## Chauve-souris
# on supprime les lignes de la matrice qui correspondent aux espèces pour lesquelles on n'a pas l'information
# phylogénétique, c'est-à-dire les espèces de chauve-souris qui ne sont pas dans l'abre :
batspecies_interet <- intersect(rownames(matrice_interactions), tree_cytb_modif$tip.label)
matrice_inter_phylo_bat <- matrice_interactions[c(batspecies_interet),]

# test de Mantel chauves-souris a echelle regionale :
test_Mantel_bat <- phylosignal_network(network = t(matrice_inter_phylo_bat), tree_A=tree_cytb_modif, method = "Jaccard_binary", correlation = "Pearson", nperm = 10000)
test_Mantel_bat <- as.data.frame(test_Mantel_bat)
test_Mantel_bat <- t(test_Mantel_bat)

# test de Mantel chauves-souris et plantes a echelle regionale :
matrice_inter_both <- matrice_inter_phylo_bat[rownames(matrice_inter_phylo_bat), colnames(matrice_inter_phylo_plant)]
#View(matrice_inter_both)

test_Mantel_both <- phylosignal_network(network = matrice_inter_both, tree_A=tree_phylo_plant, tree_B = tree_cytb_modif, method = "Jaccard_binary", correlation = "Pearson", nperm = 10000)
test_Mantel_both <- as.data.frame(test_Mantel_both)
test_Mantel_both <- t(test_Mantel_both)
test_Mantel_both

#### Etape 3 : A l'echelle locale  ####
## Chercher le signal pour chaque localisation

#on cree un vecteur avec toutes les latitudes differentes et on cree une boucle for pour explorer le signal phylogenetique a chaque latitude
latitude <- unique(data$Latitude)
longitude <- unique(data$Longitude)
L_lat <-length(latitude) #L_lat definit le nombre de latitudes differentes donc le nombre d'iterations de la boucle

#Création d'un dataframe vide pour stocker tous les résultats des tests de Mantel plantes et chauve-souris
stock_tests_both_lat <- as.data.frame(setNames(replicate(8,numeric(0), simplify = F),
                                               c("nb_A","nb_B","mantel_cor_A","pvalue_upper_A","pvalue_lower_A","mantel_cor_B","pvalue_upper_B","pvalue_lower_B")))
for (i in 1:L_lat){
  #on fait une boucle for pour toutes les localisations
  lat <- latitude[i]
  
  #pour chaque localisation on crée une nouvelle matrice d'interactions (ces étapes sont valables pour tous les tests phylogénétiques)
  submatrice <- filter(data, data$Latitude==lat)
  submatrice <- as.data.frame(submatrice) #faire transformation en dataframe sinon la suite fonctionne pas
  #On transforme les Piper en Pothomorphe, et Manilkara en Achras dans la colonne nommée "Plan_genus"
  submatrice <- subset(submatrice, submatrice$Plan_genus != "X")
  p <- nrow(submatrice)
  for (j in 1:p){
    if (submatrice$Plan_genus[j]=="Piper"){
      submatrice$Plan_genus[j]="Pothomorphe"
    }
    else if (submatrice$Plan_genus[j]=="Manilkara"){
      submatrice$Plan_genus[j]="Achras"
    }
  }
  submatrice$count = 1  #on rajoute une colonne pour pouvoir faire le pivot wider
  submatrice <- pivot_wider(submatrice, id_cols = Bat_specie, names_from = Plan_genus, values_from = count)
  submatrice <- as.data.frame(submatrice)
  submatrice[is.na(submatrice)] = 0
  #On remplace (s'il y en a) les NULL et les 0 par des 0, et on remplace les c(1,1) et les 1 par des 1 :
  for (k in 2:ncol(submatrice)){
    for (l in 1:nrow(submatrice)){
      if (submatrice[l, k]=="NULL"){submatrice[l, k]<-0}
      else if (is.list(submatrice[l, k])){submatrice[l, k]<-1}
      else if (submatrice[l, k]==0){submatrice[l, k]<-0}
      else{submatrice[l, k]<- 1}
    }
  }
  #pour être sûrs que la submatrice contient bien des nombres, et pas des listes (sinon erreurs):
  #On conserve d'abord les noms des espèces car ils vont être supprimés
  names_bats <- submatrice$Bat_specie
  #Puis on transforme la submatrice pour que ses colonnes soient des "numeric" et pas des "list"
  submatrice <- as.data.frame(lapply(submatrice, as.numeric))
  
  ## on renomme les lignes
  rownames(submatrice) = names_bats
  submatrice = submatrice[,-1]

  
  ## Matrice et arbre plantes :
  
  #on supprime dans la submatrice les genres qui ne sont pas presents dans la phylogenie (la fonction intersect ne conserve que les especes communes aux deux colonnes comparees)
  #en réalité pour les plantes tous les genres sont presents dans la phylogenie
  #donc matrice_inter_phylo_plant_lat = submatrice
  genus_interet_lat <- intersect(colnames(submatrice), tree_phylo_plant$tip.label)
  matrice_inter_phylo_plant_lat <- submatrice %>% select(all_of(genus_interet_lat))
  #View(matrice_inter_phylo_plant_lat)
  
  #On crée un nouvel arbre de plantes en enlevant les genres pas présents dans la matrice :
  tree_plants_lat <- drop.tip(tree_phylo_plant, setdiff(tree_phylo_plant$tip.label, colnames(matrice_inter_phylo_plant_lat)),
                              trim.internal = TRUE, subtree = FALSE, root.edge = 1, rooted = is.rooted(tree_phylo_plant),
                              collapse.singles = FALSE, interactive = FALSE)

  ## Matrice et arbre chauve-souris :
  
  #on supprime dans la submatrice les espèces qui ne sont pas presentes dans la phylogenie (la fonction intersect ne conserve que les especes communes aux deux colonnes comparees)
  batspecies_interet_lat <- intersect(rownames(submatrice), tree_cytb_modif$tip.label)
  matrice_inter_phylo_bat_lat <- submatrice[c(batspecies_interet_lat),]
  #View(matrice_inter_phylo_bat_lat)
  
  #On crée un nouvel arbre de chauve-souris en enlevant les espèces pas présentes dans la matrice :
  tree_bats_lat <- drop.tip(tree_cytb_modif, setdiff(tree_cytb_modif$tip.label, rownames(matrice_inter_phylo_bat_lat)), 
                            trim.internal = TRUE, subtree = FALSE, root.edge = 1, 
                            rooted = is.rooted(tree_cytb_modif), collapse.singles = FALSE, interactive = FALSE)

  ## Test du signal phylogenetique chauve-souris et plantes :
  
  matrice_inter_both_lat <- matrice_inter_phylo_bat_lat[rownames(matrice_inter_phylo_bat_lat), colnames(matrice_inter_phylo_plant_lat)]
  #View(matrice_inter_both_lat)
  
  test_Mantel_both_lat <- phylosignal_network(network = matrice_inter_both_lat, tree_A = tree_plants_lat, tree_B = tree_bats_lat, method = "Jaccard_binary", correlation = "Pearson", nperm = 10000)
  test_Mantel_both_lat <- as.data.frame(test_Mantel_both_lat)
  test_Mantel_both_lat <- t(test_Mantel_both_lat)
  stock_tests_both_lat[i,] <- test_Mantel_both_lat
} 

## Ajout des latitudes au dataframe des résultats
stock_tests_both_lat$latitude <- latitude
auteurs <- as.data.frame(auteurs)
stock_tests_both_lat <- bind_cols(stock_tests_both_lat, auteurs)
stock_tests_both_lat$Premier.auteur <- as.factor(stock_tests_both_lat$Premier.auteur)

## Résultats de tous les tests de Mantel à l'échelle locale
View(stock_tests_both_lat)

#### Etape 4 : Representation graphique ####

# on crée un vecteur significatif qui est TRUE si la pvalue < 0.05 et FALSE sinon pour les pvalue A et B
# et un autre pour la pvalue < 0.1

for (i in 1:L_lat){
  stock_tests_both_lat$significatif_A_0.05[i] <- stock_tests_both_lat$pvalue_upper_A[i] <= 0.05 ;
  stock_tests_both_lat$significatif_B_0.05[i] <- stock_tests_both_lat$pvalue_upper_B[i] <= 0.05 ;
  stock_tests_both_lat$significatif_A_0.1[i] <- stock_tests_both_lat$pvalue_upper_A[i] <= 0.1 ;
  stock_tests_both_lat$significatif_B_0.1[i] <- stock_tests_both_lat$pvalue_upper_B[i] <= 0.1 ;
  stock_tests_both_lat$nb_A_sup_5[i] <- stock_tests_both_lat$nb_A[i] >=5 ; 
  stock_tests_both_lat$nb_B_sup_5[i] <- stock_tests_both_lat$nb_B[i] >=5 ; 
  stock_tests_both_lat$nb_A_sup_10[i] <- stock_tests_both_lat$nb_A[i] >=10 ; 
  stock_tests_both_lat$nb_B_sup_10[i] <- stock_tests_both_lat$nb_B[i] >=10 ; 
}

# on crée de nouvelles colonnes pour rassembler l'information de significativité et de nombre minimum (on le fait avec 10 et 5)
for (i in 1:L_lat){
  stock_tests_both_lat$sigA0.1_and_sup10[i] <- stock_tests_both_lat$significatif_A_0.1[i] & stock_tests_both_lat$nb_A_sup_10[i]
  stock_tests_both_lat$sigB0.1_and_sup10[i] <- stock_tests_both_lat$significatif_B_0.1[i] & stock_tests_both_lat$nb_B_sup_10[i]
  stock_tests_both_lat$sigA0.1_and_sup5[i] <- stock_tests_both_lat$significatif_A_0.1[i] & stock_tests_both_lat$nb_A_sup_5[i]
  stock_tests_both_lat$sigB0.1_and_sup5[i] <- stock_tests_both_lat$significatif_B_0.1[i] & stock_tests_both_lat$nb_B_sup_5[i]
}

##on represente la pvalue en fonction de la latitude
pvalue_A_graphe <- ggplot(stock_tests_both_lat, aes(latitude, pvalue_upper_A))+
  geom_point(aes(color=sigA0.1_and_sup10))+
  theme_bw()
#pvalue_A_graphe

pvalue_B_graphe <- ggplot(stock_tests_both_lat, aes(latitude, pvalue_upper_B))+
  geom_point(aes(color=sigB0.1_and_sup10))+
  theme_bw()
#pvalue_B_graphe

## on represente le coef correlation en fonction de la latitude
coef_cor_A_graphe <-ggplot(stock_tests_both_lat, aes(latitude, mantel_cor_A))+
  geom_point(aes(color=sigA0.1_and_sup10))+
  theme_bw()
#coef_cor_A_graphe

coef_cor_B_graphe <-ggplot(stock_tests_both_lat, aes(latitude, mantel_cor_B))+
  geom_point(aes(color=sigB0.1_and_sup10))+
  theme_bw()
#coef_cor_B_graphe

## on represente le coef de correlation en fonction du nombre de taxons dans l'echantillon  
coef_cor_A_graphe_nb <- ggplot(stock_tests_both_lat, aes(nb_A, mantel_cor_A))+
  geom_point(aes(color=sigA0.1_and_sup10))+
  theme_bw()
#coef_cor_A_graphe_nb

coef_cor_B_graphe_nb <- ggplot(stock_tests_both_lat, aes(nb_B, mantel_cor_B))+
  geom_point(aes(color=sigB0.1_and_sup10))+
  theme_bw()
#coef_cor_B_graphe_nb

significatif_local <- ggarrange(pvalue_A_graphe, pvalue_B_graphe,
                                coef_cor_A_graphe, coef_cor_B_graphe, 
                                coef_cor_A_graphe_nb, coef_cor_B_graphe_nb, ncol=2, nrow=3)
ggsave(plot = significatif_local, file = "representation_pval_mantel_cor.pdf", width=25, height = 30, units = "cm")


#### Etape 5 : Representation climatique locale ####

for (i in 1:(ncol(localisation))){
  localisation[, i] <-gsub(",", ".", localisation[, i])
  localisation[, i] <- as.numeric(localisation[, i])
}

stock_tests_both_lat <- left_join(stock_tests_both_lat, localisation)
#View(stock_tests_both_lat)
write.csv2(stock_tests_both_lat, file = "resultat_test_mantel_local.csv", row.names = FALSE)

corr_A_precipitation <- ggplot(stock_tests_both_lat, aes(precipitation, mantel_cor_A))+
  geom_point(aes(color=sigA0.1_and_sup10))+
  theme_bw()
#corr_A_precipitation

corr_B_precipitation <- ggplot(stock_tests_both_lat, aes(precipitation, mantel_cor_B))+
  geom_point(aes(color=sigB0.1_and_sup10))+
  theme_bw()
#corr_B_precipitation

corr_A_temperature <- ggplot(stock_tests_both_lat, aes(temperature, mantel_cor_A))+
  geom_point(aes(color=sigA0.1_and_sup10))+
  theme_bw()
#corr_A_temperature

corr_B_temperature <- ggplot(stock_tests_both_lat, aes(temperature, mantel_cor_B))+
  geom_point(aes(color=sigB0.1_and_sup10))+
  theme_bw()
#corr_B_temperature

corr_A_NPP <- ggplot(stock_tests_both_lat, aes(npp, mantel_cor_A))+
  geom_point(aes(color=sigA0.1_and_sup10))+
  labs(x="Net Primary Production", y="coefficient de correlation A")+
  theme_bw()
#corr_A_NPP

corr_B_NPP <- ggplot(stock_tests_both_lat, aes(npp, mantel_cor_B))+
  geom_point(aes(color=sigB0.1_and_sup10))+
  labs(x="Net Primary Production", y="coefficient de correlation B")+
  theme_bw()
#corr_B_NPP

world <- map_data("world")

coordonnees_A <- ggplot()+
  geom_polygon(data=world, aes(x=long, y=lat, group=group))+
  geom_point(data = stock_tests_both_lat, aes(longitude, latitude, color=sigA0.1_and_sup10))+
  coord_cartesian(xlim = c(-100, -35), ylim = c(-30, 20))+
  theme_bw()
coordonnees_A

coordonnees_B <- ggplot()+
  geom_polygon(data=world, aes(x=long, y=lat, group=group))+
  geom_point(data = stock_tests_both_lat, aes(longitude, latitude, color=sigB0.1_and_sup10))+
  coord_cartesian(xlim = c(-100, -35), ylim = c(-30, 20))+
  theme_bw()
coordonnees_B


climatic_local <- ggarrange(coordonnees_A, coordonnees_B, corr_A_precipitation, corr_B_precipitation, 
                            corr_A_temperature, corr_B_temperature, corr_A_NPP, corr_B_NPP,
                            ncol=2, nrow=4)
ggsave(plot = climatic_local, file = "climatic_local_sig0.1_and_sup10.pdf", width=36, height = 49, units = "cm")


#### Etape 6 : Tests statistiques ####
## Pour sigA0.1_and_sup5 et sigA0.1_and_sup10, pas besoin de modèle linéaire mixte, car les latitudes concernées ne sont pas
## incluses. On fait donc des régressions linéaires simples pour sigA0.1. 

# On ne s'intéresse qu'aux latitudes significatives (p < 0.1) et avec nb_A ou nb_B > 10 (ou > 5)
# -> creation de nouveaux dataframes pour ces latitudes :
stock_tests_both_lat_sigA0.1_and_sup10 <- filter(stock_tests_both_lat, sigA0.1_and_sup10==TRUE)
stock_tests_both_lat_sigB0.1_and_sup10 <- filter(stock_tests_both_lat, sigB0.1_and_sup10==TRUE)
stock_tests_both_lat_sigA0.1_and_sup5 <- filter(stock_tests_both_lat, sigA0.1_and_sup5==TRUE)
stock_tests_both_lat_sigB0.1_and_sup5 <- filter(stock_tests_both_lat, sigB0.1_and_sup5==TRUE)

## Régression linéaire precipitations :

reg_lin_A_sup10_precipitation <- lm(mantel_cor_A~precipitation, data=stock_tests_both_lat_sigA0.1_and_sup10)
reg_lin_A_sup5_precipitation <- lm(mantel_cor_A~precipitation, data=stock_tests_both_lat_sigA0.1_and_sup5)

# Visualisation :
corr_A_precipitation_sigandmin10 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup10, aes(precipitation, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_precipitation_sigandmin10

corr_A_precipitation_sigandmin5 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup5, aes(precipitation, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_precipitation_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
plot(reg_lin_A_sup10_precipitation,2)
shapiro.test(residuals(reg_lin_A_sup10_precipitation))
# l'hypothèse de normalité est vérifiée

# on teste l'hypothèse de linéarité (avec chaque test):
plot(reg_lin_A_sup10_precipitation,1) 

# on teste l'hypothèse d'homogénéité des variances (avec chaque test):
plot(reg_lin_A_sup10_precipitation, 3) 
ncvTest(reg_lin_A_sup10_precipitation)
# l'hypothèse d'homogénéité est vérifiée

# Résultats des régressions linéaires précipitations :
summary(reg_lin_A_sup10_precipitation) # pas de lien significatif pour reg_lin_A_sup10_precipitation
summary(reg_lin_A_sup5_precipitation) # pas de lien significatif pour reg_lin_A_sup5_precipitation

## Régression linéaire température :

reg_lin_A_sup10_temperature <- lm(mantel_cor_A~temperature, data=stock_tests_both_lat_sigA0.1_and_sup10)
reg_lin_A_sup5_temperature <- lm(mantel_cor_A~temperature, data=stock_tests_both_lat_sigA0.1_and_sup5)

# Visualisation :
corr_A_temperature_sigandmin10 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup10, aes(temperature, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_temperature_sigandmin10

corr_A_temperature_sigandmin5 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup5, aes(temperature, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_temperature_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
plot(reg_lin_A_sup5_temperature,2)
shapiro.test(residuals(reg_lin_A_sup5_temperature))
# l'hypothèse de normalité est vérifiée

# on teste l'hypothèse de linéarité (avec chaque test):
plot(reg_lin_A_sup5_temperature,1) 

# on teste l'hypothèse d'homogénéité des variances (avec chaque test):
plot(reg_lin_A_sup5_temperature, 3) 
ncvTest(reg_lin_A_sup5_temperature)
# l'hypothèse d'homogénéité est vérifiée

# Résultats des régressions linéaires température :
summary(reg_lin_A_sup10_temperature) # pas de lien significatif pour reg_lin_A_sup10_temperature
summary(reg_lin_A_sup5_temperature) # pas de lien significatif pour reg_lin_A_sup5_temperature

## Régression linéaire NPP :

reg_lin_A_sup10_npp <- lm(mantel_cor_A~npp, data=stock_tests_both_lat_sigA0.1_and_sup10)
reg_lin_A_sup5_npp <- lm(mantel_cor_A~npp, data=stock_tests_both_lat_sigA0.1_and_sup5)

# Visualisation :
corr_A_npp_sigandmin10 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup10, aes(npp, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_npp_sigandmin10

corr_A_npp_sigandmin5 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup5, aes(npp, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_npp_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
plot(reg_lin_A_sup5_npp,2)
shapiro.test(residuals(reg_lin_A_sup5_npp))
# l'hypothèse de normalité est vérifiée

# on teste l'hypothèse de linéarité (avec chaque test):
plot(reg_lin_A_sup10_npp,1) 

# on teste l'hypothèse d'homogénéité des variances (avec chaque test):
plot(reg_lin_A_sup5_npp, 3) 
ncvTest(reg_lin_A_sup5_npp)
# l'hypothèse d'homogénéité est vérifiée

# Résultats des régressions linéaires NPP :
summary(reg_lin_A_sup10_npp) # pas de lien significatif pour reg_lin_A_sup10_npp
summary(reg_lin_A_sup5_npp) # pas de lien significatif pour reg_lin_A_sup5_npp

## Régression linéaire latitude :

reg_lin_A_sup10_latitude <- lm(mantel_cor_A~latitude, data=stock_tests_both_lat_sigA0.1_and_sup10)
reg_lin_A_sup5_latitude <- lm(mantel_cor_A~latitude, data=stock_tests_both_lat_sigA0.1_and_sup5)

# Visualisation :
corr_A_latitude_sigandmin10 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup10, aes(latitude, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_latitude_sigandmin10

corr_A_latitude_sigandmin5 <- ggplot(stock_tests_both_lat_sigA0.1_and_sup5, aes(latitude, mantel_cor_A)) +
  geom_point(color="#00BFC4") +
  geom_smooth(method = "lm") +
  theme_bw()
corr_A_latitude_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
plot(reg_lin_A_sup5_latitude,2)
shapiro.test(residuals(reg_lin_A_sup5_latitude))
# l'hypothèse de normalité est vérifiée

# on teste l'hypothèse de linéarité (avec chaque test):
plot(reg_lin_A_sup5_latitude,1) 

# on teste l'hypothèse d'homogénéité des variances (avec chaque test):
plot(reg_lin_A_sup5_latitude, 3) 
ncvTest(reg_lin_A_sup5_latitude)
# l'hypothèse d'homogénéité est vérifiée

# Résultats des régressions linéaires latitude :
summary(reg_lin_A_sup10_latitude) # pas de lien significatif pour reg_lin_A_sup10_latitude
summary(reg_lin_A_sup5_latitude) # pas de lien significatif pour reg_lin_A_sup5_latitude

## Pour sigB0.1_and_sup5 et sigB0.1_and_sup10, les latitudes 8,8666 et 4,8798 correspondent aux références
## 11 et 12, concernées par le problème d'indépendance (premier auteur identique). 
## On fait donc un modèle linéaire mixte pour sigB0.1 en prenant en compte le premier auteur comme "random variable".

## Modèle linéaire mixte precipitations :

lmm_B_sup10_precipitation <- lme(mantel_cor_B~precipitation, random = ~ 1 |Premier.auteur,
                                 data=stock_tests_both_lat_sigB0.1_and_sup10, method = "REML")

lmm_B_sup5_precipitation <- lme(mantel_cor_B~precipitation, random = ~ 1 |Premier.auteur,
                                data=stock_tests_both_lat_sigB0.1_and_sup5, method = "REML")

# Visualisation :
corr_B_precipitation_sigandmin10 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup10, aes(precipitation, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=precipitation), as.data.frame(t(fixef(lmm_B_sup10_precipitation))))
  theme_bw()
corr_B_precipitation_sigandmin10

corr_B_precipitation_sigandmin5 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup5, aes(precipitation, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=precipitation), as.data.frame(t(fixef(lmm_B_sup5_precipitation))))+
  theme_bw()
corr_B_precipitation_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
qqnorm(residuals(lmm_B_sup10_precipitation))
qqline(residuals(lmm_B_sup10_precipitation))
shapiro.test(residuals(lmm_B_sup10_precipitation))
# l'hypothèse de normalité est vérifiée

## Résultats du modèle linéaire mixte précipitations :
summary(lmm_B_sup10_precipitation) 
summary(lmm_B_sup5_precipitation) # problème : on obtient NaN en p-value.
# Cause probable : trop petit nombre de données par rapport au nombre de groupes (8 groupes pour 9 données), résultant en
# une trop faible variation intra-groupes par rapport à la variation inter-groupes, ce qui nous donne 0 degrés de liberté
# pour la valeur de la pente, et donc pose problème pour le calcul de la p-value, d'où le NaN.
# Solution: se référer à la p-value obtenue avec un modèle linéaire simple lm().

# essai avec une régression linéaire normale :
reg_lin_B_sup10_precipitation <- lm(mantel_cor_B~precipitation, data=stock_tests_both_lat_sigB0.1_and_sup10)
reg_lin_B_sup5_precipitation <- lm(mantel_cor_B~precipitation, data=stock_tests_both_lat_sigB0.1_and_sup5)
summary(reg_lin_B_sup10_precipitation)
summary(reg_lin_B_sup5_precipitation)
# on obtient bien une p-value avec ce test. Elle n'est pas significative (valeurs de pente et d'intercept semblables à
# celles obtenues avec le lmm, donc on peut faire confiance à cette p-value). Le R² associé est faible.

## Modèle linéaire mixte température :

lmm_B_sup10_temperature <- lme(mantel_cor_B~temperature, random = ~1 | Premier.auteur, 
                               data=stock_tests_both_lat_sigB0.1_and_sup10)
lmm_B_sup5_temperature <- lme(mantel_cor_B~temperature, random =~1 | Premier.auteur, 
                              data=stock_tests_both_lat_sigB0.1_and_sup5)

# Visualisation :
corr_B_temperature_sigandmin10 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup10, aes(temperature, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=temperature), as.data.frame(t(fixef(lmm_B_sup10_temperature))))+
  theme_bw()
corr_B_temperature_sigandmin10

corr_B_temperature_sigandmin5 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup5, aes(temperature, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=temperature), as.data.frame(t(fixef(lmm_B_sup5_temperature))))+
  theme_bw()
corr_B_temperature_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
qqnorm(residuals(lmm_B_sup10_temperature))
qqline(residuals(lmm_B_sup10_temperature))
shapiro.test(residuals(lmm_B_sup10_temperature))
# l'hypothèse de normalité est vérifiée

## Résultats du modèle linéaire mixte température :
summary(lmm_B_sup10_temperature) 
summary(lmm_B_sup5_temperature) # problème: on obtient NaN en p-value. Cause et solution : voir dans "precipitations".

# essai avec une régression linéaire normale :
reg_lin_B_sup10_temperature <- lm(mantel_cor_B~temperature, data=stock_tests_both_lat_sigB0.1_and_sup10)
reg_lin_B_sup5_temperature <- lm(mantel_cor_B~temperature, data=stock_tests_both_lat_sigB0.1_and_sup5)
summary(reg_lin_B_sup10_temperature)
summary(reg_lin_B_sup5_temperature)
# on obtient bien une p-value avec ce test. Elle n'est pas significative, et le R² associé est faible.

## Modèle linéaire mixte NPP :

lmm_B_sup10_npp <- lme(mantel_cor_B~npp, random =~1 |Premier.auteur, data=stock_tests_both_lat_sigB0.1_and_sup10)
lmm_B_sup5_npp <- lme(mantel_cor_B~npp, random = ~1 | Premier.auteur, data=stock_tests_both_lat_sigB0.1_and_sup5)

# Visualisation :
corr_B_npp_sigandmin10 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup10, aes(npp, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=npp), as.data.frame(t(fixef(lmm_B_sup10_npp))))
theme_bw()
corr_B_npp_sigandmin10

corr_B_npp_sigandmin5 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup5, aes(npp, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=npp), as.data.frame(t(fixef(lmm_B_sup5_npp))))
theme_bw()
corr_B_npp_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
qqnorm(residuals(lmm_B_sup10_npp))
qqline(residuals(lmm_B_sup10_npp))
shapiro.test(residuals(lmm_B_sup10_npp))
# l'hypothèse de normalité est vérifiée

# Résultats du modèle linéaire mixte NPP :
summary(lmm_B_sup10_npp) 
summary(lmm_B_sup5_npp) # problème : on obtient NaN en p-value. Cause et solution : voir dans "precipitations".

# essai avec une régression linéaire normale :
reg_lin_B_sup10_npp <- lm(mantel_cor_B~npp, data=stock_tests_both_lat_sigB0.1_and_sup10)
reg_lin_B_sup5_npp <- lm(mantel_cor_B~npp, data=stock_tests_both_lat_sigB0.1_and_sup5)
summary(reg_lin_B_sup10_npp)
summary(reg_lin_B_sup5_npp)
# on obtient bien une p-value avec ce test. Elle n'est pas significative, et le R² associé est faible.

## Modèle linéaire mixte latitude :

lmm_B_sup10_latitude <- lme(mantel_cor_B~latitude, random = ~1 | Premier.auteur, data=stock_tests_both_lat_sigB0.1_and_sup10)
lmm_B_sup5_latitude <- lme(mantel_cor_B~latitude, random =~1 | Premier.auteur, data=stock_tests_both_lat_sigB0.1_and_sup5 )

# Visualisation :
corr_B_latitude_sigandmin10 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup10, aes(latitude, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=latitude), as.data.frame(t(fixef(lmm_B_sup10_latitude))))
theme_bw()
corr_B_latitude_sigandmin10

corr_B_latitude_sigandmin5 <- ggplot(stock_tests_both_lat_sigB0.1_and_sup5, aes(latitude, mantel_cor_B)) +
  geom_point(color="#00BFC4") +
  geom_abline(aes(intercept=`(Intercept)`, slope=latitude), as.data.frame(t(fixef(lmm_B_sup5_latitude))))
theme_bw()
corr_B_latitude_sigandmin5

# on teste l'hypothèse de normalité (avec chaque test):
qqnorm(residuals(lmm_B_sup10_latitude))
qqline(residuals(lmm_B_sup10_latitude))
shapiro.test(residuals(lmm_B_sup10_latitude))
# l'hypothèse de normalité est vérifiée

## Résultats du modèle linéaire mixte latitude :
summary(lmm_B_sup10_latitude) 
summary(lmm_B_sup5_latitude) # problème : on obtient NaN en p-value. Cause et solution : voir dans "precipitations".

# essai avec une régression linéaire normale :
reg_lin_B_sup10_latitude <- lm(mantel_cor_B~latitude, data=stock_tests_both_lat_sigB0.1_and_sup10)
reg_lin_B_sup5_latitude <- lm(mantel_cor_B~latitude, data=stock_tests_both_lat_sigB0.1_and_sup5)
summary(reg_lin_B_sup10_latitude)
summary(reg_lin_B_sup5_latitude)
# on obtient bien une p-value avec ce test. Elle n'est pas significative, et le R² associé est faible.

## Conclusion de cette partie : les tests statistiques n'ont pas permis de mettre en évidence une corrélation à l'échelle locale
## entre signal phylogénétique et conditions climatiques/latitude, ni pour le signal phylogénétique des plantes, 
## ni pour celui des chauve-souris.

regressions_lineaires <- ggarrange(corr_A_precipitation_sigandmin10, corr_B_precipitation_sigandmin10, 
                                   corr_A_temperature_sigandmin10, corr_B_temperature_sigandmin10,
                                   corr_A_npp_sigandmin10, corr_B_npp_sigandmin10,
                                   corr_A_latitude_sigandmin10, corr_B_latitude_sigandmin10,
                                   ncol=2, nrow=4)
ggsave(plot = regressions_lineaires, file = "regressions_lineaires.pdf", width=36, height = 49, units = "cm")


#### Etape 7 : Representation correlation A par rapport correlation B ####

## on represente le coef de correlation de B en fonction du coef de correlation de A : 
# visuellement on change la couleur et la forme en fonction du nombre de taxons et de la significativité du test 
A_B_graphe_pvalue_sup10 <- ggplot(stock_tests_both_lat, aes(mantel_cor_A, mantel_cor_B))+
  geom_point(aes(color=sigB0.1_and_sup10, shape=sigA0.1_and_sup10))+
  labs(y="Coeff de corrélation du test de Mantel chauves-souris", x="Coeff de corrélation du test de Mantel plantes")+
  theme_bw()
A_B_graphe_pvalue_sup10

A_B_graphe_pvalue_sup5 <- ggplot(stock_tests_both_lat, aes(mantel_cor_A, mantel_cor_B))+
  geom_point(aes(color=sigB0.1_and_sup5, shape=sigA0.1_and_sup5))+
  labs(y="Coeff de corrélation du test de Mantel chauves-souris", x="Coeff de corrélation du test de Mantel plantes")+
  theme_bw()
A_B_graphe_pvalue_sup5
# Ce qu'on voit sur ce graphe c'est que le signal phylogénétique n'est pas symétrique : il y a beaucoup de latitudes où on a
# du signal phylogénétique pour les chauve-souris, mais pas pour les plantes (donc à ces latitudes, les chauve-souris proches
# entre elles ont tendance à interagir avec des plantes proches, mais aux mêmes latitudes, les plantes proches n'interagissent
# pas avec des chauve-souris proches). Il n'y a que 2 latitudes pour lesquelles c'est l'inverse, et 2 latitudes pour lesquelles
# c'est symétrique.
# Cela signifie que les causes du signal phylogénétique ne sont pas les mêmes chez les chauve-souris et chez les plantes.

corrA_corrB_graphe <- ggarrange(A_B_graphe_pvalue_sup10, ncol=1, nrow=1)
ggsave(plot = corrA_corrB_graphe, file = "corrA_corrB_graphe.pdf", width=15, height = 15, units = "cm")

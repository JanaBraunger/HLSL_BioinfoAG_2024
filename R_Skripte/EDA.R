# Proteom-Daten geben darüber Auskunft, welche Proteine wie stark in der gemessenen Probe exprimiert, d.h. vorhanden, sind.
  
# Falls ihr euch eigene Datensätze suchen wollt: PRIDE ist eine öffentliche Datenbank für Proteom-Datensätze.
# https://www.ebi.ac.uk/pride/
# 
# 
# In der "exploratory data analysis (EDA)" Phase versuchen wir erstmal uns einen generellen Überblick über die Daten zu verschaffen.
# Wie viele Gene werden exprimiert?
#   Wie viele Proben wurden gemessen?
#   Welche Vergleiche könnten interessant sein?
#   Welche Proteine werden besonders stark/schwach exprimiert?
#   Welche Proteine unterscheiden sich stark zwischen verschiedenen Proben?
#   
#   Verschiedene Arten von Graphen sind dabei hilfreich, z.B.:
# - Balkendiagramme (bar plots)
# - Histogramme
# - box plots
# - Heatmaps
# - volcano plot
# 
#   Außerdem benötigen wir einfache statistische Tests:
#   - T-Test
#   




# Datensätze einlesen:
# verschiedene Lungenkrebs-Zelllinien wurden gemessen (RNAseq und Protein)
# Wir fokussieren uns vorerst nur auf die Protein-Daten.
protein_daten = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_proteome.csv", row.names = 1)


head(protein_daten)
#column names = cell line names
#row names = gene names
#values = Intensität der Proteine

dim(protein_daten)

#Beschreibung der Zell-Linien
metadata = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_description.csv", row.names = 1)
head(metadata)
dim(metadata)
# Welche Informationen haben wir zu den Zellen?


#Stimmen die Namen überein?
#Die %in% und table() Funktionen
colnames(protein_daten)
metadata$model_name


colnames(protein_daten) %in% metadata$model_name
table(colnames(protein_daten) %in% metadata$model_name)


#Welche Organe sind abgedeckt? Die unique() Funktion
unique(metadata$tissue) #nur Lunge

# Für welche weiteren Spalten könnte die unique() Funktion interessant sein?

# Die table() Funktion
table(metadata$tissue_status)
#81x metastasis, 79x Tumor
table(metadata$cancer_type)
#84 x NSCLC
#57 x SCLC
# Es gibt verschiedene Lungenkrebs-Arten

table(metadata$gender)
table(metadata$ethnicity)
table(metadata$smoking_status)

# Histogramme:
# Histogramme zeigen die Verteilung von Daten, d.h. welche Werte wie oft vorkommen
# Nur für numerische Daten geeignet!
hist(metadata$age_at_sampling)

# verschiedene Daten-Typen (Klassen)
# character (Text), numeric (Zahlen), logical (True/False)
class(metadata$age_at_sampling)
class(protein_daten$A427)

hist(as.numeric(metadata$age_at_sampling))


# 
#Wie viele Proteine wurden in jeder Probe gemessen?
# Funktionen: !is.na(), sum(), apply()
# barplot() und boxplot()




barplot(apply(rna_seq, 2, function(x) sum(!is.na(x))))
boxplot(apply(rna_seq, 2, function(x) sum(!is.na(x))))




# Welche Fragestellungen ergeben sich aus der EDA?
#Ein interessanter Vergleich könnte zum Beispiel sein, ob sich beim Lungenkrebs die primären Tumore von den Metastasen unterscheiden?
#Oder NSCLC von SCLC? Oder Raucher von Nicht-Rauchern? Oder Frauen vs. Männer?


# Heatmaps können dabei helfen, verschiedene Gruppen zu identifizieren:
# Die "base" Funktion heatmap() erfüllt diesen Zweck auch, aber ich finde das package pheatmap etwas schöner und einfacher
# 
#install.packages("pheatmap")
library(pheatmap)
pheatmap(na.omit(protein_daten))

#Wir können die columns mit hilfreichen Informationen beschriften:
#Dazu nutzen wir die metadaten
rownames(metadata) = metadata$model_name
metadata$age_at_sampling = as.numeric(metadata$age_at_sampling)

pheatmap(na.omit(protein_daten), scale = "row", show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100),
         annotation_col = metadata[,c("cancer_type"), drop = F])

# Falls die Farben nicht ideal sind:
unique(metadata$cancer_type)
bessere_farben = list(cancer_type = c())

pheatmap(na.omit(protein_daten), scale = "row", show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100),
         annotation_col = metadata[,c("cancer_type"), drop = F], annotation_colors = bessere_farben)

# Wie könnten wir mehrere Informationen gleichzeitig beschriften?


# Statistik:Welche Gene sind unterschiedlich stark exprimiert?
# z.B. SCLC vs. NSCLC
# 
# Wir wollen den "fold-change" und den "p-value" für jedes Gen herausfinden
# --> Wie viel mal mehr/weniger ist jedes Gen exprimiert und ist dies signifikant?

rownames(metadata) # Das sind alle Namen der Zelllinien
head(protein_daten[,c("ABC1","DMS114")])
head(protein_daten[,rownames(metadata)[1:4]])

#wir wollen nur Zelllinien behalten, die entweder SCLC oder NSCLC sind:
rownames(metadata)[metadata$cancer_type == "Small Cell Lung Carcinoma"]

protein_sclc_nsclc = protein_daten[,c(rownames(metadata)[metadata$cancer_type == "Small Cell Lung Carcinoma"],rownames(metadata)[metadata$cancer_type == "Non-Small Cell Lung Carcinoma"])]
head(protein_sclc_nsclc)

pheatmap(na.omit(protein_sclc_nsclc), scale = "row", show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100),
         annotation_col = metadata[,c("cancer_type"), drop = F], annotation_colors = bessere_farben)



# t.test() Funktion:
#

# fold-change = mean(group1)/mean(group2)
# 1. Berechne für jedes Gen den Mittelwert in Gruppe 1 (SCLC)
# mean(), apply()
# 2. Berechne für jedes Gen den Mittelwert in Gruppe 2 (NSCLC)
# 3. Teile die beiden Mittelwerte für jedes Gen durcheinander
# mean1-mean2
# Wie sieht die Verteilung dieser "fold-changes" aus?

# Volcano Plot

# Gene Set Enrichment


## TRANSKRIPTOM:
# Sehen wir ähnliche Unterschiede im Proteom?
# Wie sehr korrelieren Transkriptom und Proteom?
# Gibt es Gene, die zwar stark exprimiert, aber anscheinend nicht translatiert werden? Oder andersherum?




#RNA-seq Daten: Lungenkrebs-Zelllinien
rna_seq = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_RNA_seq.csv", row.names = 1)




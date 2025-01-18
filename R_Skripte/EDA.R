# Proteom-Daten geben darüber Auskunft, welche Proteine wie stark in der gemessenen Probe exprimiert, d.h. vorhanden, sind.
  
# Falls ihr euch eigene Datensätze suchen wollt: PRIDE ist eine öffentliche Datenbank für Proteom-Datensätze.
# https://www.ebi.ac.uk/pride/
# 
# 
# In der "exploratory data analysis (EDA)" Phase versuchen wir erstmal uns einen generellen Überblick über die Daten zu verschaffen.
#   Wie viele Proteine wurden gemessen?
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
#   Außerdem benötigen wir einfache statistische Mittel:
#   - Mittelwert
#   - Standardabweichung
#   - T-Test
#   




# Datensätze einlesen:
# verschiedene Lungenkrebs-Zelllinien wurden gemessen (RNAseq und Protein)
# Wir fokussieren uns vorerst nur auf die Protein-Daten.
# 
#AUFGABE: Lade den Datensatz von github herunter, ändere den Pfad und lies den Datensatz ein.
protein_daten = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_proteome.csv", row.names = 1)


head(protein_daten)
#column names = cell line names
#row names = gene names
#values = Intensität der Proteine

dim(protein_daten)

#Beschreibung der Zell-Linien
##AUFGABE: Lade die Metadaten von github herunter, ändere den Pfad und lies den Datensatz ein.

metadata = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_description.csv", row.names = 1)
head(metadata)
dim(metadata)
# Welche Informationen haben wir zu den Zellen?




#Welche Organe sind abgedeckt? Die unique() Funktion
unique(metadata$tissue) #nur Lunge

# Für welche weiteren Spalten könnte die unique() Funktion interessant sein?

# Die table() Funktion
table(metadata$tissue_status)
#81x metastasis, 79x Tumor

# AUFGABE: Wende die unique() und table() Funktionen auf weitere Spalten an, die dich interessieren.



# Histogramme:
# Histogramme zeigen die Verteilung von Daten, d.h. welche Werte wie oft vorkommen
# Nur für numerische Daten geeignet!
hist(metadata$age_at_sampling)


# verschiedene Daten-Typen (Klassen)
# character (Text), numeric (Zahlen), logical (True/False)
class(metadata$age_at_sampling)


hist(as.numeric(metadata$age_at_sampling))

# AUFGABE: Wie könntest du Verteilung von nicht-numerischen Daten darstellen? z.B. von metadata$ethnicity
# Tipp: table() + barplot()


# Boxplots sind auch hilfreich, um die Verteilung von Daten anzuschauen und zu vergleichen:
boxplot(as.numeric(metadata$age_at_sampling))

# AUFGABE: Wie könnten wir die Altersverteilung nicht von allen Lungenkrebs-Fällen, sondern nur NSCLC, darstellen?
# TEIL2: Wie könnten wir die ALtersverteilung von NSCLC und SCLC nebeneinander vergleichen?


# Welche Fragestellungen ergeben sich aus der EDA?
#Ein interessanter Vergleich könnte zum Beispiel sein, ob sich beim Lungenkrebs die primären Tumore von den Metastasen unterscheiden?
#Oder NSCLC von SCLC? Oder Raucher von Nicht-Rauchern? Oder Frauen vs. Männer?


# Heatmaps können dabei helfen, verschiedene Gruppen zu identifizieren:
# Die "base" Funktion heatmap() erfüllt diesen Zweck auch, aber ich finde das package pheatmap etwas schöner und einfacher
# 
#install.packages("pheatmap")
library(pheatmap)
pheatmap(na.omit(protein_daten), show_rownames = F, scale = "row")

# Bessere Farben:

#Wir können die columns mit hilfreichen Informationen beschriften:
#Dazu nutzen wir die metadaten
rownames(metadata) = metadata$model_name
metadata$age_at_sampling = as.numeric(metadata$age_at_sampling)


# Auch von diesen Beschriftungen können wir die Farben anpassen:
unique(metadata$cancer_type)
bessere_farben = list(cancer_type = c())



# Statistik:Welche Gene sind unterschiedlich stark exprimiert?
# z.B. SCLC vs. NSCLC
# 
# Wir wollen den "fold-change" und den "p-value" für jedes Gen herausfinden
# --> Wie viel mal mehr/weniger ist jedes Gen exprimiert und ist dieser Unterschied signifikant?

rownames(metadata) # Das sind alle Namen der Zelllinien

# Wiederholung: Auf bestimmte Teile einer Tabelle zugreifen.



#wir wollen nur Zelllinien behalten, die entweder SCLC oder NSCLC sind:
rownames(metadata)[metadata$cancer_type == "Small Cell Lung Carcinoma"]

protein_sclc_nsclc = protein_daten[,c(rownames(metadata)[metadata$cancer_type == "Small Cell Lung Carcinoma"],rownames(metadata)[metadata$cancer_type == "Non-Small Cell Lung Carcinoma"])]
head(protein_sclc_nsclc)

pheatmap(na.omit(protein_sclc_nsclc), scale = "row", show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100),
         annotation_col = metadata[,c("cancer_type"), drop = F], annotation_colors = bessere_farben)



# t.test() Funktion:
# Einfaches Beispiel:

# Erklärung apply():

# Jetzt wollen wir die t.test() Funktion auf alle Reihen unseres Datensatzes anwenden
# Dabei soll SCLC gegen NSCLC getestet werden
# In welchen Spalten sind die SCLC/NSCLC Zellen?
length(rownames(metadata)[metadata$cancer_type == "Small Cell Lung Carcinoma"])
# In den ersten 57 Spalten sind SCLC
SCLC_spalten = c(1:57)
length(rownames(metadata)[metadata$cancer_type == "Non-Small Cell Lung Carcinoma"])
# In den nächsten 84 Spalten sind NSCLC
NSCLC_spalten = c(58:(57+84))

protein_sclc_nsclc$p_Wert = apply(protein_sclc_nsclc, 1, function(x) t.test(x[SCLC_spalten], x[NSCLC_spalten])$p.value)
# Der T-Test kann nicht mit zu vielen fehlenden Werten umgehen.

protein_sclc_nsclc = na.omit(protein_sclc_nsclc)
#Jetzt müsste es Funktionieren

# AUFGABE: Stelle die Verteilung der p-Werte dar!


# fold-change = mean(SCLC)-mean(NSCLC)
# AUFGABE:
# 1. Berechne für jedes Protein den Mittelwert in Gruppe 1 (SCLC_spalten)
# Tipps: mean(), apply()
# 2. Berechne für jedes Gen den Mittelwert in Gruppe 2 (NSCLC_spalten)
# 3. Ziehe jetzt die beiden Mittelwerte für jedes Protein voneinander ab
# 
protein_sclc_nsclc$FC = protein_sclc_nsclc$SCLC_mean - protein_sclc_nsclc$NSCLC_mean

# Wie sieht die Verteilung dieser "fold-changes" aus?


# Volcano Plot
plot(protein_sclc_nsclc$FC, -log10(protein_sclc_nsclc$p_Wert))

# interessante Proteine beschriften:

top10up = head(rownames(protein_sclc_nsclc)[order(protein_sclc_nsclc$FC, decreasing = T)],10)
text(protein_sclc_nsclc[top10up, "FC"], -log10(protein_sclc_nsclc[top10up, "p_Wert"]), labels = top10up)


# AUFGABE: Suche dir 2 Proteine aus und beschrifte diese im Plot.

# Gene set enrichment:
#install.packages("BiocManager")
#BiocManager::install("ReactomePA")
#install.packages("org.Hs.eg.db")
library(ReactomePA)
library(org.Hs.eg.db)

gene_list = protein_sclc_nsclc$FC
names(gene_list) = mapIds(org.Hs.eg.db, rownames(protein_sclc_nsclc), 'ENTREZID', 'SYMBOL')
gene_list_sorted = gene_list[order(gene_list, decreasing = T)]



GSEA_Ergebnis = gsePathway(gene_list_sorted, 
                           pvalueCutoff = 0.05)

dotplot(GSEA_Ergebnis, showCategory=20, x = "NES", color = "p.adjust", title = "Unterschiede zwischen SCLC und NSCLC")

## TRANSKRIPTOM:
# Sehen wir ähnliche Unterschiede in den Transkriptom-Daten?
# Wie sehr korrelieren Transkriptom und Proteom?

#RNA-seq Daten: Lungenkrebs-Zelllinien
rna_seq = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_RNA_seq.csv", row.names = 1)

# 



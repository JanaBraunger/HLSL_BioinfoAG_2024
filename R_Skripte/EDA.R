# Transkriptom-Daten geben uns Auskunft darüber, welche Gene wie stark exprimiert werden. Nach der Prozessierung der Roh-Sequenzierdaten, erhalten wir als Ergebnis eine Tabelle, in der jedes Gen für jede gemessene Probe einen Expressionswert hat. Dieser ist z.B. TPM (transcripts per kilobase million) oder FPKM (fragments per kilobase million).
# 
# Falls ihr euch eigene Datensätze suchen wollt:
# GEO ist eine öffentliche Datenbank für RNAseq-Datensätze.
# https://www.ncbi.nlm.nih.gov/geo/
  
# PRIDE ist eine öffentliche Datenbank für Proteom-Datensätze.
# https://www.ebi.ac.uk/pride/
# 
# 
# In der "exploratory data analysis (EDA)" Phase versuchen wir erstmal uns einen generellen Überblick über die Daten zu verschaffen.
# Wie viele Gene werden exprimiert?
#   Wie viele Proben wurden gemessen?
#   Welche Vergleiche könnten interessant sein?
#   Welche Gene werden besonders stark/schwach exprimiert?
#   Welche Gene unterscheiden sich stark zwischen verschiedenen Proben?
#   
#   Verschiedene Arten von Graphen sind dabei hilfreich, z.B.:
#   - Balkendiagramme (bar plots)
# - Histogramme
# - box plots
# - Heatmaps
# - volcano plot
# 
# 
# Wir fangen an, indem wir die benötigten packages installieren und laden

#install.packages("pheatmap")
library(pheatmap)


# Datensätze einlesen:
#RNA-seq Daten: Lungenkrebs-Zelllinien
rna_seq = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_RNA_seq.csv", row.names = 1)

head(rna_seq)
#column names = cell line names
#row names = gene names
#values = TPM of gene expression

dim(rna_seq)

#Proteom Daten:
#proteome = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_proteome.csv", row.names = 1)
#head(proteome)
#column names = cell line names
#row names = gene names
#values = protein intensities
#dim(proteome)

#Beschreibung der Zell-Linien
metadata = read.csv("K:/Ergebnisse/LS_testing/1_PhD/HLSL/transcriptome_proteome_dataset/lung_cell_lines_description.csv", row.names = 1)
head(metadata)
dim(metadata)

#do we have information for all cell lines?
#Die %in% und table() Funktionen
table(colnames(rna_seq) %in% metadata$model_name)


#table(colnames(proteome) %in% metadata$model_name)



#welche Informationen haben wir über die Zell-Linien?
colnames(metadata)

#Welche Organe sind abgedeckt? Die unique() Funktion
unique(metadata$tissue) #nur Lunge

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

# Histogram
hist(as.numeric(metadata$age_at_sampling))
# 
#Wie viele mRNAs wurden in jeder Probe gemessen?
# Funktionen: !is.na(), sum(), apply()
# barplot() und boxplot()
barplot(apply(rna_seq, 2, function(x) sum(!is.na(x))))
boxplot(apply(rna_seq, 2, function(x) sum(!is.na(x))))


# Wie viele Proteine wurden in jeder Probe gefunden?
#(apply(proteome, 2, function(x) sum(!is.na(x))))
#boxplot(apply(proteome, 2, function(x) sum(!is.na(x))))

# Welche Fragestellungen ergeben sich aus der EDA?
#Ein interessanter Vergleich könnte zum Beispiel sein, ob sich beim Lungenkrebs die primären Tumore von den Metastasen unterscheiden?
#Oder NSCLC von SCLC? Oder Raucher von Nicht-Rauchern? Oder Frauen vs. Männer?


# Heatmaps können dabei helfen, verschiedene Gruppen zu identifizieren:
pheatmap(na.omit(rna_seq), scale = "row", show_rownames = F, show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100))

#Wir können die columns mit hilfreichen Informationen beschriften:
#Dazu nutzen wir die metadaten
rownames(metadata) = metadata$model_name
metadata$age_at_sampling = as.numeric(metadata$age_at_sampling)

pheatmap(na.omit(rna_seq), scale = "row", show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100),
         annotation_col = metadata[,c("cancer_type", "smoking_status", "age_at_sampling", "gender"), drop = F])

# Wie könnten wir mehrere Informationen gleichzeitig beschriften?

#pheatmap(na.omit(proteome), scale = "row",show_rownames = F, show_colnames = F, color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100), annotation_col = metadata[,"smoking_status", drop = F])

# Statistik:Welche Gene sind unterschiedlich stark exprimiert?
# z.B. SCLC vs. NSCLC
# 
# Wir wollen den "fold-change" und den "p-value" für jedes Gen herausfinden
# --> Wie viel mal mehr/weniger ist jedes Gen exprimiert und ist dies signifikant?

rownames(metadata)
head(rna_seq[,c("ABC1","DMS114")])
head(rna_seq[,rownames(metadata)[1:4]])

rna_seq_sclc_nsclc = rna_seq[,c(rownames(metadata)[metadata$cancer_type == "Small Cell Lung Carcinoma"],rownames(metadata)[metadata$cancer_type == "Non-Small Cell Lung Carcinoma"])]

pheatmap(na.omit(rna_seq_sclc_nsclc), scale = "row", show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-5,5,10/100),
         annotation_col = metadata[,c("cancer_type"), drop = F])



# t.test() Funktion:
#

# fold-change = mean(group1)/mean(group2)
# 1. Berechne für jedes Gen den Mittelwert in Gruppe 1 (SCLC)
# 2. Berechne für jedes Gen den Mittelwert in Gruppe 2 (NSCLC)
# 3. Teile die beiden Mittelwerte für jedes Gen durcheinander
# Wie sieht die Verteilung dieser "fold-changes" aus?

# Volcano Plot

# Gene Set Enrichment


## PROTEOM:
# Sehen wir ähnliche Unterschiede im Proteom?
# Wie sehr korrelieren Transkriptom und Proteom?
# Gibt es Gene, die zwar stark exprimiert, aber anscheinend nicht translatiert werden? Oder andersherum?








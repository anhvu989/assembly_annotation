import pandas as pd
min_reads = 300000
species = "Lactobacillus_crispatus"

df = pd.read_csv("AllSamples_Kraken_NoNorm_March2023.csv", sep=',')
count = 0
for i in df.columns:
    if df.loc[species,i] >= min_reads:# for all columns (samples)
        with open("list_samples.txt", "a") as outfile:
	        print(i + "\t " + str(df.loc[species, i] / df[i].sum()), file = outfile)
        count += 1
print(count)
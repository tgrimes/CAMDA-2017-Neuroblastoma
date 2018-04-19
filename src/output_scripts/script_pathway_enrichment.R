# KEGG pathway enrichment downloaded from DAVID.
# Save pathway enrichment tables, formatted for Latex.

df <- read.table("output/tables/gene_sets/DAVID/OS_all_homo_sapiens.txt",
                 sep = "\t", header = TRUE) %>%
  select(Term, Count, Pop.Hits, PValue, Benjamini) %>%
  mutate(Term = gsub("hsa.*:", "", Term),
         PValue = round(PValue, 3),
         Benjamini = round(Benjamini, 3)) %>%
  rename("Pathway" = Term, "Count" = Count, "Size" = Pop.Hits, 
         "P-value" = PValue, "BH" = Benjamini) %>%
  filter(BH < 0.05)
write.table(df, file = paste("output/tables/OS_enrichment_homo_sapiens.txt", sep = ""), 
            quote = FALSE, sep = " \t& ", eol = " \\\\\n", row.names = FALSE)

df2 <- read.table("output/tables/gene_sets/DAVID/OS_all_gene_cards.txt",
                 sep = "\t", header = TRUE) %>%
  select(Term, Count, Pop.Hits, PValue, Benjamini) %>%
  mutate(Term = gsub("hsa.*:", "", Term),
         PValue = round(PValue, 3),
         Benjamini = round(Benjamini, 3)) %>%
  rename("Pathway" = Term, "Count" = Count, "Size" = Pop.Hits, 
         "P-value" = PValue, "BH" = Benjamini) %>%
  filter(`P-value` < 0.05)
write.table(df2, file = paste("output/tables/OS_enrichment_gene_cards.txt", sep = ""), 
            quote = FALSE, sep = " \t& ", eol = " \\\\\n", row.names = FALSE)

df3 <- read.table("output/tables/gene_sets/DAVID/EFS_all_homo_sapiens.txt",
                 sep = "\t", header = TRUE) %>%
  select(Term, Count, Pop.Hits, PValue, Benjamini) %>%
  mutate(Term = gsub("hsa.*:", "", Term),
         PValue = round(PValue, 3),
         Benjamini = round(Benjamini, 3)) %>%
  rename("Pathway" = Term, "Count" = Count, "Size" = Pop.Hits, 
         "P-value" = PValue, "BH" = Benjamini) %>%
  filter(BH < 0.05)
# write.table(df3, file = paste("output/tables/EFS_enrichment_homo_sapiens.txt", sep = ""), 
#             quote = FALSE, sep = " \t& ", eol = " \\\\\n", row.names = FALSE)

df4 <- read.table("output/tables/gene_sets/DAVID/EFS_all_gene_cards.txt",
                  sep = "\t", header = TRUE) %>%
  select(Term, Count, Pop.Hits, PValue, Benjamini) %>%
  mutate(Term = gsub("hsa.*:", "", Term),
         PValue = round(PValue, 3),
         Benjamini = round(Benjamini, 3)) %>%
  rename("Pathway" = Term, "Count" = Count, "Size" = Pop.Hits, 
         "P-value" = PValue, "BH" = Benjamini) %>%
  filter(`P-value` < 0.05)
write.table(df4, file = paste("output/tables/EFS_enrichment_gene_cards.txt", sep = ""), 
            quote = FALSE, sep = " \t& ", eol = " \\\\\n", row.names = FALSE)


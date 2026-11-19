# Top genes in the somatic main set: samples with a call per gene.
library(dplyr)
library(ggplot2)
library(readxl)
library(here)

TOP_GENES <- 15
theme_set(theme_bw(base_size = 12))

calls <- read_excel(here("data/20261116_somatic_main_calls_review.xlsx"), guess_max = 5000)
genes <- calls %>%
  group_by(Gene) %>%
  summarise(n_samples = n_distinct(Sample.ID), n_calls = n(), n_variants = n_distinct(ID)) %>%
  arrange(desc(n_samples), desc(n_calls))
print(genes, n = Inf)

fig <- genes %>% head(TOP_GENES) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene))) %>%
  ggplot(aes(Gene, n_samples)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n_samples), hjust = -0.3, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Top genes in the somatic main set",
       subtitle = paste0(nrow(calls), " calls in ", n_distinct(calls$Sample.ID), " samples"),
       x = NULL, y = "Samples with a call") +
  theme(panel.grid.major.y = element_blank())
ggsave(here("figures/clinical/top_genes_main.pdf"), width = 6, height = 5)

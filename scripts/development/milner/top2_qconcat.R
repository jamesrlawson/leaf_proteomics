# QconCAT processing

require(tidyverse)

rawdat <- read_tsv('ignore/milner/Milner MultiQuant results SVS 180628.txt')

# extrac only labelled ions

labelled_dat <- rawdat[grep("(\\[\\+08\\]|\\[\\+10\\]|\\[\\+06\\])", rawdat$`Component Name`),]

# find mean (across all samples) area ranks for labelled ions

mean_ion_ranks  <- labelled_dat %>% 
                   select(`Sample Name`, `Component Group Name`, `Component Name`,`Area`) %>%
                   group_by(`Sample Name`, `Component Group Name`) %>%
                   arrange(desc(Area)) %>%
                   mutate(within_component_rank = 1:length(Area)) %>%
                   group_by(`Component Name`) %>%
                   summarise(mean_ion_rank = mean(within_component_rank)) %>%
                   full_join(labelled_dat, by = "Component Name") %>%
                   arrange(`Sample Name`, `Component Group Name`)


# select top2 ions for each peptide

top2 <- mean_ion_ranks %>% 
        group_by(`Sample Name`, `Component Group Name`) %>%
        top_n(2, mean_ion_rank) # %>% 
        # select(`Sample Name`, `Component Group Name`, `Component Name`, `Area`, mean_ion_rank)
  


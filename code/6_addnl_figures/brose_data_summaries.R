#Brose figure
#Ana Miller-ter kuile
#April 16

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "ggforce", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------


data <- read.csv(here("data",
                      "raw_data",
                      "5_brose_data",
                      "data.csv"))


# Tidy data ---------------------------------------------------------------

# select predation only of inverts
predate <- data %>%
  filter(interaction.type == "predacious") %>%
  filter(con.metabolic.type == "invertebrate")

#give diet assignment higher categories
predate_1 <- predate %>%
  mutate(con.size.method = case_when(con.size.method %in% c("measurement", "regression") ~ "Direct measurement",
                                     con.size.method == "measurement  published account  regression" ~ "Measured and literature",
                                     con.size.method == "published accounts" ~ "Literature", 
                                     TRUE ~ "Not reported")) %>%
  mutate(res.size.method = case_when(res.size.method %in% c("measurement",
                                                            "field measurement",
                                                            "regression",
                                                            "measurement: individuals are field-sampled, then masses are derived by linear regressions (weight-length regression with measured lengths)") ~ "Direct measurement*",
                                     res.size.method == "measurement  published account  regression" ~ "Measured* and literature",
                                     res.size.method == "published accounts" ~ "Literature",
                                     TRUE ~ "Not reported")) %>%
  mutate(res.taxonomy.level = case_when(res.taxonomy.level == "Class" ~ "class",
                                        TRUE ~ res.taxonomy.level)) %>%
  mutate(con.taxonomy.level = case_when(con.taxonomy.level == "Class" ~ "class",
                                        TRUE ~ con.taxonomy.level)) %>%
  mutate(link.methodology = case_when(link.methodology %in% c("direct observation ",
                                                              "Feeding trial",
                                                              "Gut content") ~ "Observed",
                                      link.methodology %in% c("published account internet field natural history",
                                                              "field;  gut/stomach analysis; expert;  feeding trial;  natural history",
                                                              "published literature; direct observation; expert knowledge",
                                                              "gut analysis, published account",
                                                              "gut content, extrapolated from similar taxa",
                                                              "Field; gut/stomach analysis; feeding trial; published account") ~ "Inferred and observed",
                                      link.methodology %in% c("published account", "expert ", "expert", " expert", "Published account","published account, expert") ~ "Inferred",
                                      TRUE ~ "Not reported"))


# Figures ------------------------------------------------------------------

#link assignment
 links <- predate_1 %>%
   mutate(link.methodology = factor(link.methodology, 
                                    levels = c("Not reported", 
                                               "Inferred",
                                               "Inferred and observed",
                                               "Observed"))) %>%
ggplot(aes(x = link.methodology)) +
  geom_bar() +
  labs(title = "A) Link assignment method", y = "Number of Interactions") +
  scale_y_continuous(limits = c(0,120000), breaks=c(0,25000, 50000, 75000, 100000)) +
  coord_flip() +
  theme_bw() +
   theme(axis.text = element_text(size = 15),
         axis.title = element_blank(),
         plot.title = element_text(size = 20),
         axis.text.x = element_blank())

#predator size assignment
preds <- predate_1 %>%
  mutate(con.size.method = factor(con.size.method, 
                                  levels = c("Not reported",
                                             "Literature",
                                             "Measured and literature",
                                             "Direct measurement"))) %>%
ggplot(aes(x = con.size.method)) +
  geom_bar() +
  labs(title = "B) Predator body size method", y = "Number of Interactions") +
  scale_y_continuous(limits = c(0,120000), breaks=c(0,25000, 50000, 75000, 100000)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_blank(),
        plot.title = element_text(size = 20),
        axis.text.x = element_blank())

#prey assignment method
prey <- predate_1 %>%
  mutate(res.size.method = factor(res.size.method,
                                  levels = c("Not reported", 
                                             "Literature",
                                             "Measured* and literature",
                                             "Direct measurement*"))) %>%
ggplot(aes(x = res.size.method)) +
  geom_bar() +
  labs(title = "C) Prey body size method", y = "Number of Interactions") +
  scale_y_continuous(limits = c(0,120000), breaks=c(0,25000, 50000, 75000, 100000)) +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 20))

links /
  preds + prey

# summary stats -----------------------------------------------------------

predate_1 %>%
  group_by(link.methodology) %>%
  summarise(total = n())

predate_1 %>%
  tally()

(15504+2089)/131025

(2089)/131025

predate_1 %>%
  group_by(res.size.method) %>%
  tally()

(112622+352)/131025

predate_1 %>%
  group_by(con.size.method) %>%
  tally()



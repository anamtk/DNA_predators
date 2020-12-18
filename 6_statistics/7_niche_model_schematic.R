
# Load packages -----------------------------------------------------------

package.list <- c("tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# create data #

x<-as.data.frame(c(1,2,3,4))

x <- x %>%
  rename("x" = "c(1, 2, 3, 4)")

#plot wedge of sizes for niche model schematic
ggplot(x, aes(x = x, y = x)) +
  geom_line(aes(y = 0.43*x), size = 0.75) +
  geom_line(aes(y = 0.73*x), linetype = "dashed", size = 0.75) +
  geom_line(aes(y = 0.13*x), linetype = "dashed", size = 0.75) +
  geom_ribbon(aes(ymin = 0.13*x, ymax = 0.73*x), alpha = 0.6) +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

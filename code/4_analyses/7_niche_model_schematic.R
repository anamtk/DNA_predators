
# Load packages -----------------------------------------------------------

package.list <- c("tidyverse", "patchwork",
                  "calecopal")

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
  geom_line(aes(y = 0.73*x), color = "#bdbdbd", size = 0.75) +
  geom_line(aes(y = 0.13*x), color = "#bdbdbd", size = 0.75) +
  #geom_ribbon(aes(ymin = 0.13*x, ymax = 0.73*x), fill = "#bdbdbd", alpha = 0.6) +
  geom_line(aes(y = 0.43*x), size = 0.75, color = "#bdbdbd") +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 25))


#plot wedge of sizes for niche model schematic
bigger <- ggplot(x, aes(x = x, y = x)) +
  geom_line(aes(y = 0.73*x), color = "#bdbdbd", linetype = "dashed", size = 0.75) +
  geom_line(aes(y = 0.83*x), color = "#bdbdbd", size = 0.75) +
  geom_line(aes(y = 0.13*x), color = "#bdbdbd", size = 0.75) +
  #geom_ribbon(aes(ymin = 0.13*x, ymax = 0.83*x), fill = "#bdbdbd", alpha = 0.6) +
  geom_line(aes(y = 0.43*x), size = 0.75, linetype = "dashed") +
  geom_line(aes(y = 0.43*x + 0.1), size = 0.75) +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))



#plot wedge of sizes for niche model schematic
smaller <- ggplot(x, aes(x = x, y = x)) +
  geom_line(aes(y = 0.73*x), color = "#bdbdbd", size = 0.75) +
  geom_line(aes(y = 0.13*x), color = "#bdbdbd", linetype = "dashed", size = 0.75) +
  geom_line(aes(y = 0.03*x), color = "#bdbdbd", size = 0.75) +
  #geom_ribbon(aes(ymin = 0.03*x, ymax = 0.73*x), fill = "#bdbdbd", alpha = 0.6) +
  geom_line(aes(y = 0.43*x), size = 0.75, linetype = "dashed") +
  geom_line(aes(y = 0.43*x - 0.1), size = 0.75) +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

basic | (bigger / smaller) +
  plot_layout(widths = c(2, 1))


#plot wedge of sizes for niche model schematic
ggplot(x, aes(x = x, y = x)) +
  geom_line(aes(y = 0.73*x), 
            color = "#bdbdbd",
            linetype= "dashed",
            size = 1) +
  geom_line(aes(y = 0.13*x), 
            color = "#bdbdbd", 
            linetype = "dashed",
            size = 1) +
  #geom_ribbon(aes(ymin = 0.13*x, ymax = 0.73*x), fill = "#bdbdbd", alpha = 0.6) +
  geom_line(aes(y = 0.43*x), 
            size = 1) +
  geom_line(aes(y = 0.43*x - 0.1), 
            size = 1, 
            color = "#067D8D") +
  geom_line(aes(y = 0.63*x), 
            size = 1, 
            color = "#067D8D",
            linetype = "dashed",
            alpha = 0.6) +
  geom_line(aes(y = 0.23*x), 
            size = 1, 
            color = "#067D8D",
            linetype = "dashed",
            alpha = 0.6) +
  geom_line(aes(y = 0.03*x), 
            color = "#EEB00C", 
            linetype = "dashed", 
            alpha = 0.6, 
            size = 1) +
  geom_line(aes(y = 0.43*x + 0.1), 
            size = 1, 
            color = "#EEB00C") +
  geom_line(aes(y = 0.83*x), 
            color = "#EEB00C",
            linetype = "dashed",
            alpha = 0.6,
            size = 1) +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 25))

#LRS: "#C70000" 
# SCY: "#EA7700" 
#NEO: "#EEB00C" 
#CEN: "#C68A2C" 
#SME: "#89742F" 
#EUB: "#496C3C" 
#PHH: "#158D8E"
#PAN: "#067D8D" 
#HEV: "#114C54"

#plot wedge of sizes for niche model schematic
basic <- ggplot(x, aes(x = x, y = x)) +
  geom_line(aes(y = 0.73*x), 
            color = "#bdbdbd",
            linetype= "dashed",
            size = 1) +
  geom_line(aes(y = 0.13*x), 
            color = "#bdbdbd", 
            linetype = "dashed",
            size = 1) +
  #geom_ribbon(aes(ymin = 0.13*x, ymax = 0.73*x), fill = "#bdbdbd", alpha = 0.6) +
  geom_line(aes(y = 0.43*x), 
            size = 1) +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 25))

#plot wedge of sizes for niche model schematic
smaller <- ggplot(x, aes(x = x, y = x)) +
  geom_line(aes(y = 0.73*x), 
            color = "#bdbdbd",
            linetype= "dashed",
            size = 1) +
  geom_line(aes(y = 0.13*x), 
            color = "#bdbdbd", 
            linetype = "dashed",
            size = 1) +
  #geom_ribbon(aes(ymin = 0.13*x, ymax = 0.73*x), fill = "#bdbdbd", alpha = 0.6) +
  geom_line(aes(y = 0.43*x), 
            size = 1) +
  geom_line(aes(y = 0.43*x - 0.1), 
            size = 1, 
            color = "#067D8D") +
  geom_line(aes(y = 0.63*x), 
            size = 1, 
            color = "#067D8D",
            linetype = "dashed",
            alpha = 0.6) +
  geom_line(aes(y = 0.23*x), 
            size = 1, 
            color = "#067D8D",
            linetype = "dashed",
            alpha = 0.6) +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 25))

#plot wedge of sizes for niche model schematic
bigger <- ggplot(x, aes(x = x, y = x)) +
  geom_line(aes(y = 0.73*x), 
            color = "#bdbdbd",
            linetype= "dashed",
            size = 1) +
  geom_line(aes(y = 0.13*x), 
            color = "#bdbdbd", 
            linetype = "dashed",
            size = 1) +
  #geom_ribbon(aes(ymin = 0.13*x, ymax = 0.73*x), fill = "#bdbdbd", alpha = 0.6) +
  geom_line(aes(y = 0.43*x), 
            size = 1) +
  geom_line(aes(y = 0.03*x), 
            color = "#EEB00C", 
            linetype = "dashed", 
            alpha = 0.6, 
            size = 1) +
  geom_line(aes(y = 0.43*x + 0.1), 
            size = 1, 
            color = "#EEB00C") +
  geom_line(aes(y = 0.83*x), 
            color = "#EEB00C",
            linetype = "dashed",
            alpha = 0.6,
            size = 1) +
  labs(x = expression(~log[10]~"predator mass"),
       y = expression(~log[10]~"prey mass")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 25))

basic + plot_spacer() + smaller + plot_spacer() + bigger + plot_spacer()

basic + smaller + bigger

plot_spacer() + bigger + basic + plot_spacer() + plot_spacer() + smaller +
  plot_layout(ncol = 2)

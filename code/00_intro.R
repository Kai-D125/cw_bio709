#' DESCRIPTION:
#' Script for introductory work

# in-class ----------------------------------------------------------------

#non-parametric tests

library(tidyverse)

x <- c(3.2, 5.0, 10.0, 100, 50)
y <- c(1.0, 3.0, 2.0, 2.1, 1.2)

df_xy <- tibble(group = c(rep("x", length(x)),
                 rep("y", length(y))),
                 value = c(x,y)
       )

df_xy %>% 
  ggplot(aes(x = group,
             y = value)) +
  geom_boxplot()

t.test(x, y)

#non-parametric t-test
wilcox.test(x, y)

#ANOVA for >2 groups assuming normal distribution
aov(weight ~ group, 
    data = PlantGrowth)

#non-parametric ANOVA = kruskal-wallis

kruskal.test(weight ~ group,
             data = PlantGrowth)

#confidence interval
m <- lm(Petal.Length ~ Petal.Width,
        data = iris)

summary(m)

confint(m)

#correlation
x <- rnorm(100, mean = 0, sd = 1)
y <- rnorm(100, mean = 0.8 * x, sd = 1)

plot(x, y)

#parametric correlation = pearson
cor.test(x, y)

#non-parametric correlation = spearman
cor.test(x, y, method = "spearman")

#covariance
cov(x, y)

cov(x, y)/sd(x) * sd(y)

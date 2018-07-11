## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  y ~ x, data = my_data

## ----eval=FALSE----------------------------------------------------------
#  y ~ x|z, data = my_data

## ---- message=FALSE------------------------------------------------------
library(effects, warn.conflicts = FALSE)
library(kableExtra, warn.conflicts = FALSE)
library(knitr, warn.conflicts = FALSE)
library(multcomp, warn.conflicts = FALSE)
library(papeR, warn.conflicts = FALSE)
library(pubh, warn.conflicts = FALSE)

## ------------------------------------------------------------------------
data(Hodgkin)
Hodgkin$Ratio <- Hodgkin$CD4/Hodgkin$CD8
kable(head(Hodgkin))

## ------------------------------------------------------------------------
estat(~ CD4, data = Hodgkin)

## ------------------------------------------------------------------------
kable(estat(~ CD4, data = Hodgkin, label = "CD4^+^ T cells (cells / mm^3^)"))

## ------------------------------------------------------------------------
kable(estat(~ Ratio|Group, data = Hodgkin, label = "CD4^+^/CD8^+^ T cells"))

## ------------------------------------------------------------------------
kable(estat(Ratio ~ Group, data = Hodgkin, label = "CD4^+^/CD8^+^ T cells"))

## ------------------------------------------------------------------------
var.test(Ratio ~ Group, data = Hodgkin)

## ------------------------------------------------------------------------
qq_plot(~ Ratio|Group, data = Hodgkin, ylab = "CD4+ / CD8+ T cells", aspect = 1)

## ------------------------------------------------------------------------
wilcox.test(Ratio ~ Group, data = Hodgkin)

## ------------------------------------------------------------------------
strip_error(Ratio ~ Group, data = Hodgkin, ylab = "CD4+ / CD8+ T cells")

## ------------------------------------------------------------------------
strip_error(Ratio ~ Group, data = Hodgkin, ylab = "CD4+ / CD8+ T cells", ylim = c(0, 4.5)) +
  layer(panel.segments(1, 4, 2, 4, lwd=1.5)) +
  layer(panel.text(1.5, 4.1, "**"))

## ------------------------------------------------------------------------
data(birthwt)
birthwt$smoke <- factor(birthwt$smoke, labels = c("Non-smoker", "Smoker"))
birthwt$race <- factor(birthwt$race, labels = c("White", "African American", "Other"))

## ------------------------------------------------------------------------
kable(gen_bst_df(bwt ~ smoke, data = birthwt))

## ------------------------------------------------------------------------
bar_error(bwt ~ smoke, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3500))

## ------------------------------------------------------------------------
qq_plot(~ bwt|smoke, data = birthwt, ylab = "Birth weight (g)", aspect = 1)

## ------------------------------------------------------------------------
t.test(bwt ~ smoke, data = birthwt)

## ------------------------------------------------------------------------
bar_error(bwt ~ smoke, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3500)) +
  layer(panel.segments(1, 3300, 2, 3300, lwd=1.5)) +
  layer(panel.text(1.5, 3350, "**"))

## ------------------------------------------------------------------------
bar_error(bwt ~ smoke|race, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3800),
          col = c("gray80", "gray30"))

## ------------------------------------------------------------------------
strip_error(bwt ~ smoke|race, data = birthwt, ylab = "Birth weight (g)", cex = 0.3)

## ------------------------------------------------------------------------
model_bwt <- lm(bwt ~ smoke + race, data = birthwt)
glm_coef(model_bwt)

## ------------------------------------------------------------------------
unadj <- glm_coef(model_bwt, labels = c("Constant",
                                      "Smoking: smoker - non-smoker",
                                      "Race: African American - White",
                                      "Race: Other - White"))
kable(unadj)

## ------------------------------------------------------------------------
model_glht <- glht(model_bwt, linfct  = mcp(race = "Tukey"))
xymultiple(model_glht)

## ------------------------------------------------------------------------
race_glht <- xymultiple(model_glht, plot = FALSE)
kable(race_glht)

## ------------------------------------------------------------------------
final <- unadj
final[, 5] <- as.character(final[, 5])
race_glht[, 5] <- as.character(race_glht[, 5])
final[3:4, 3:5] <- race_glht[1:2, 3:5]

## ------------------------------------------------------------------------
kable(final)

## ------------------------------------------------------------------------
plot(Effect(c("race", "smoke"), model_bwt), main = NULL, aspect = 3/4, 
     multiline = TRUE, ylab = "Birth weight (g)", xlab = "Race/Ethnicity")

## ------------------------------------------------------------------------
data(Bernard)
kable(head(Bernard))

## ------------------------------------------------------------------------
kable(summarise(Bernard, type = "fac", variables = "treat", group = "fate",
                test = FALSE, labels = "Treatment"))  %>%
  add_header_above(c(" " = 3, "Alive" = 2, " " = 1, "Dead" = 2))

## ------------------------------------------------------------------------
dat <- matrix(c(84, 140 , 92, 139), nrow = 2, byrow = TRUE)
epiR::epi.2by2(as.table(dat))

## ------------------------------------------------------------------------
contingency(fate ~ treat, data = Bernard)

## ----eval=FALSE----------------------------------------------------------
#  outcome ~ stratum/exposure, data = my_data

## ------------------------------------------------------------------------
data(oswego, package = "epitools")
oswego$ill <- factor(oswego$ill)
oswego$sex <- factor(oswego$sex)
oswego$chocolate.ice.cream <- factor(oswego$chocolate.ice.cream)
mhor(ill ~ sex/chocolate.ice.cream, data = oswego)

## ------------------------------------------------------------------------
odds_trend(Cancer ~ Age, data = Macmahon)

## ------------------------------------------------------------------------
Freq <- c(1739, 8, 51, 22)
BCG <- gl(2, 1, 4, labels=c("Negative", "Positive"))
Xray <- gl(2, 2, labels=c("Negative", "Positive"))
tb <- data.frame(Freq, BCG, Xray)
tb <- expand_df(tb)

diag_test(BCG ~ Xray, data=tb)

## ------------------------------------------------------------------------
diag_test2(22, 51, 8, 1739)

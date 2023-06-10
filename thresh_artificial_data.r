library("ggplot2")
library("EbayesThresh")
library("wavethresh")

#Artificial data
# Donoho e Johnstone (1994) introduced four artificial functions to represent different
# features of inhomogeneous signals. We add white noise to these functions to 
#illustrate the performance of two Bayesian wavelet shrinkage techniques.
# We use these functions sampled at 1024 points equally spaced on [0, 1].

v <- DJ.EX()
t <- (1:1024)/1024

blocks <- v$blocks
bumps <- v$bumps
heavi <- v$heavi
doppler <- v$doppler

df <- data.frame(blocks, bumps, heavi, doppler, t)

plot_blocks <- ggplot(df, aes(x = t, y = blocks)) + 
  geom_line(color="#00AFB5", size = 1) + labs(x = "", y = "Blocks") +
  theme_classic()

plot_bumps <- ggplot(df, aes(x = t, y = bumps)) + 
  geom_line(color="#F46036", size = 1) + labs(x = "", y = "Bumps") +
  theme_classic()

plot_heavi <- ggplot(df, aes(x = t, y = heavi)) + 
  geom_line(color="#FF0F80", size = 1) + labs(x = "", y = "Heavisine") +
  theme_classic()

plot_doppler <- ggplot(df, aes(x = t, y = doppler)) + 
  geom_line(color="#04724D", size = 1) + labs(x = "", y = "Doppler") +
  theme_classic()


# Adding independent normally distributed noise to the functions. 
# We establish a signal-to-noise ratio of 4.

SNR <- 4 #defining our signal-to-noise ratio

dopsig <- sd(v$doppler)
bumsig <- sd(v$bumps)
blosig <- sd(v$blocks)
heasig <- sd(v$heavi)

noisedopsig <- dopsig/SNR
noisebumsig <- bumsig/SNR
noiseblosig <- blosig/SNR 
noiseheasig <- heasig/SNR 
noiseexasig <- exasig/SNR


edop <- rnorm(length(v$doppler), mean = 0, sd = noisedopsig)
ebum <- rnorm(length(v$bumps), mean = 0, sd = noisebumsig)
eblo <- rnorm(length(v$blocks), mean = 0, sd = noiseblosig)
ehea <- rnorm(length(v$heavi), mean = 0, sd = noiseheasig)
eexa <- rnorm(length(example.1()$y), mean = 0, sd = noiseexasig)

#Creating our noisy data
ndop <- v$doppler + edop
nbum <- v$bumps + ebum
nblo <- v$blocks + eblo
nhea <- v$heavi + ehea

df <- data.frame(blocks, bumps, heavi, doppler, t, ndop, nbum, nblo, nhea)

plot_nblo <- ggplot(df, aes(x = t, y = nblo)) + 
  geom_line(color="#00AFB5") + labs(x = "", y = "Noisy Blocks") +
  theme_classic()

plot_nbum <- ggplot(df, aes(x = t, y = nbum)) + 
  geom_line(color="#F46036") + labs(x = "", y = "Noisy Bumps") +
  theme_classic()

plot_nhea <- ggplot(df, aes(x = t, y = nhea)) + 
  geom_line(color="#FF0F80") + labs(x = "", y = "Noisy Heavisine") +
  theme_classic()

plot_ndop <- ggplot(df, aes(x = t, y = ndop)) + 
  geom_line(color="#04724D") + labs(x = "", y = "Noisy Doppler") +
  theme_classic()

#Transforming the data with the Wavelet transform
family_choice <- "DaubLeAsymm" 

#Number of vanishing moments
number_vm <- 10

wdndop <- wd(ndop, filter.number = number_vm, family = family_choice)
wdnbum <- wd(nbum, filter.number = number_vm, family = family_choice)
wdnblo <- wd(nblo, filter.number = number_vm, family = family_choice)
wdnhea <- wd(nhea, filter.number = number_vm, family = family_choice)


#Using the threshold approaches to denoise the data
#BayesThresh arguments
type_Bayesthreshold <- "soft" #"soft" or "hard"

by_level_Bayesthreshold <- FALSE # F:global threshold computed.
#V: a threshold is computed and applied separately to each scale level.

#See Wavelet thresholding via Bayesian approach (1998) for details on specifying 
#hyperparameters.
alpha_Bayesthreshold <- 0.5 
beta_Bayesthreshold <- 1


#Doppler
thrwdndopun <- threshold(wdndop,levels= 3:(nlevelsWT(wdndop)-1),
                         type = type_Bayesthreshold, policy = "BayesThresh",
                         by.level = by_level_Bayesthreshold, dev = madmad,
                         boundary = FALSE, alpha = alpha_Bayesthreshold,
                         beta = beta_Bayesthreshold, C1 = NA, C2 = NA, 
                         C1.start = 100)
thrwdndopeb <- ebayesthresh.wavelet(wdndop, prior = "laplace",
                                    threshrule = "median")

#Bumps
thrwdnbumun <- threshold(wdnbum,levels= 3:(nlevelsWT(wdnbum)-1),
                         type = type_Bayesthreshold, policy = "BayesThresh",
                         by.level = by_level_Bayesthreshold, dev = madmad,
                         boundary = FALSE, alpha = alpha_Bayesthreshold,
                         beta = beta_Bayesthreshold, C1 = NA, C2 = NA, 
                         C1.start = 100)
thrwdnbumeb <- ebayesthresh.wavelet(wdnbum, prior = "laplace",
                                    threshrule = "median")
#Blocks
thrwdnbloun <- threshold(wdnblo, levels= 3:(nlevelsWT(wdnblo)-1),
                         type = type_Bayesthreshold, policy = "BayesThresh",
                         by.level = by_level_Bayesthreshold, dev = madmad,
                         boundary = FALSE, alpha = alpha_Bayesthreshold,
                         beta = beta_Bayesthreshold, C1 = NA, C2 = NA, 
                         C1.start = 100)
thrwdnbloeb <- ebayesthresh.wavelet(wdnblo, prior = "laplace",
                                    threshrule = "median")
#Heavisine
thrwdnheaun <- threshold(wdnhea, levels= 3:(nlevelsWT(wdnhea)-1),
                         type = type_Bayesthreshold, policy = "BayesThresh",
                         by.level = by_level_Bayesthreshold, dev = madmad,
                         boundary = FALSE, alpha = alpha_Bayesthreshold,
                         beta = beta_Bayesthreshold, C1 = NA, C2 = NA, 
                         C1.start = 100)

thrwdnheaeb <- ebayesthresh.wavelet(wdnhea, prior = "laplace",
                                    threshrule = "median")

##Taking the inverse DWT

#Doppler
inthrwdndopun <- wr(thrwdndopun)
inthrwdndopeb <- wr(thrwdndopeb)

#Bumps
inthrwdnbumun <- wr(thrwdnbumun)
inthrwdnbumeb <- wr(thrwdnbumeb)

#Blocks
inthrwdnbloun <- wr(thrwdnbloun)
inthrwdnbloeb <- wr(thrwdnbloeb)

#Heavisine
inthrwdnheaun <- wr(thrwdnheaun)
inthrwdnheaeb <- wr(thrwdnheaeb)

#Plots


#Doppler
dops1 <- c(inthrwdndopun,doppler)
dops2 <- c(inthrwdndopeb,doppler)
type1 <- c(rep("Wavethashed",length(inthrwdndopun)),
           rep("Doppler",length(doppler)))
type2 <- c(rep("Empirical",length(inthrwdndopeb)),
           rep("Doppler",length(doppler)))

type1 <- as.factor(type1)
type2 <- as.factor(type2)

df_dops <- data.frame(dops1, dops2, type1, type2, t)

g_dops1 <- ggplot(df_dops, aes(x = t, y = dops1, group = type1, 
                              color=type1)) + 
  geom_line(aes(linetype=type1), size = 1) + labs(x = "",
                                                  y = "Denoised doppler") +
  theme_classic() + 
  scale_color_manual(values=c('black','#04724D')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))

g_dops2 <- ggplot(df_dops, aes(x = t, y = dops2, group = type2, 
                               color=type2)) + 
  geom_line(aes(linetype=type2), size = 1) + labs(x = "", 
                                                  y = "Denoised doppler") +
  theme_classic() + 
  scale_color_manual(values=c('black','#04724D')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))


# Bumps
bumps1 <- c(inthrwdnbumun,bumps)
bumps2 <- c(inthrwdnbumeb,bumps)
type1 <- c(rep("Wavethashed",length(inthrwdnbumun)),
           rep("Bumps",length(bumps)))
type2 <- c(rep("Empirical",length(inthrwdnbumeb)),
           rep("Bumps",length(bumps)))

type1 <- as.factor(type1)
type2 <- as.factor(type2)

df_bumps <- data.frame(bumps1, bumps2, type1, type2, t)

g_bumps1 <- ggplot(df_bumps, aes(x = t, y = bumps1, group = type1, 
                                 color=type1)) + 
  geom_line(aes(linetype=type1), size = 1) + labs(x = "", 
                                                  y = "Denoised bumps") +
  theme_classic() + 
  scale_color_manual(values=c('black','#F46036')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))

g_bumps2 <- ggplot(df_bumps, aes(x = t, y = bumps2, group = type2, 
                                 color=type2)) + 
  geom_line(aes(linetype=type2), size = 1) + labs(x = "", 
                                                  y = "Denoised bumps") +
  theme_classic() + 
  scale_color_manual(values=c('black','#F46036')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))



# Blocks
block1 <- c(inthrwdnbloun,blocks)
block2 <- c(inthrwdnbloeb,blocks)
type1 <- c(rep("Wavethashed",length(inthrwdnbloun)),
           rep("Blocks",length(blocks)))
type2 <- c(rep("Empirical",length(inthrwdnbloeb)),
           rep("Blocks",length(blocks)))

type1 <- as.factor(type1)
type2 <- as.factor(type2)

df_block <- data.frame(block1, block2, type1, type2, t)

g_block1 <- ggplot(df_block, aes(x = t, y = block1, group = type1, 
                              color=type1)) + 
  geom_line(aes(linetype=type1), size = 1) + labs(x = "",
                                                  y = "Denoised blocks") +
  theme_classic() + 
  scale_color_manual(values=c('black','#00AFB5')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))

g_block2 <- ggplot(df_block, aes(x = t, y = block2, group = type2, 
                                 color=type2)) + 
  geom_line(aes(linetype=type2), size = 1) + labs(x = "",
                                                  y = "Denoised blocks") +
  theme_classic() + 
  scale_color_manual(values=c('black','#00AFB5')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))



# Heavisine
heavi1 <- c(inthrwdnheaun,heavi)
heavi2 <- c(inthrwdnheaeb,heavi)
type1 <- c(rep("Wavethashed",length(inthrwdnheaun)),
           rep("Blocks",length(heavi)))
type2 <- c(rep("Empirical",length(inthrwdnheaeb)),
           rep("Blocks",length(heavi)))

type1 <- as.factor(type1)
type2 <- as.factor(type2)

df_heavi <- data.frame(heavi1, heavi2, type1, type2, t)

g_heavi1 <- ggplot(df_heavi, aes(x = t, y = heavi1, group = type1, 
                                 color=type1)) + 
  geom_line(aes(linetype=type1), size = 1) + labs(x = "",
                                                  y = "Denoised Heavisine") +
  theme_classic() + 
  scale_color_manual(values=c('black','#FF0F80')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))

g_heavi2 <- ggplot(df_heavi, aes(x = t, y = heavi2, group = type2, 
                                 color=type2)) + 
  geom_line(aes(linetype=type2), size = 1) + labs(x = "",
                                                  y = "Denoised Heavisine") +
  theme_classic() + 
  scale_color_manual(values=c('black','#FF0F80')) + 
  guides(color = "none", linetype = "none") +
  scale_linetype_manual(values=c("dotted","solid"))









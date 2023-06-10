library("ggplot2")
library("EbayesThresh")
library("wavethresh")


#Application to the Inductance plethysmography data
#The inductance plethysmography data set is a section of plethysmograph recording, 
#made by the Department of Anaesthesia, Bristol University, to measure the flow of 
#air during breathing of patients under general anaesthetic (Nason, 1996). 
#The data contains 4096 observations.

data(ipd)
plot(ipd)

data <- as.numeric(ipd)
t <- seq(1,length(ipd),1)

#Transforming the data with the Wavelet transform
family_choice <- "DaubLeAsymm" #"DaubLeAsymm", "DaubExPhase", "Coiflets"

#Number of vanishing moments
number_vm <- 10

wddata <- wd(data, filter.number = number_vm, family = family_choice)

#Using the threshold approaches to denoise the data
#wavethresh
wavethresh <- threshold(wddata, levels= 3:(nlevelsWT(wddata)-1),
                         type = "soft", policy = "BayesThresh",
                         by.level = FALSE, dev = madmad,
                         boundary = FALSE, alpha = 0.5,
                         beta = 1, C1 = NA, C2 = NA, 
                         C1.start = 100)
#Ebayes
ebayes <- ebayesthresh.wavelet(wddata, prior = "laplace", 
                               threshrule = "median")

##Taking the inverse DWT
#Doppler
inve_wavethresh <- wr(wavethresh)
inve_ebayes <- wr(ebayes)

df <- data.frame(t, data, inve_wavethresh, inve_ebayes)

g_plethysmography <- ggplot(df, aes(x = t, y = data)) + 
  geom_line(linetype="solid", color="#00AFB5", size = 0.5) + 
  labs(x = "Index", y = "Voltage") +
  theme_classic() 

g_wavethresh <- ggplot(df, aes(x = t, y = inve_wavethresh)) + 
  geom_line(linetype="solid", color="#04724D", size = 0.7) + 
  labs(x = "Index", y = "Voltage") +
  theme_classic() 

g_ebayes <- ggplot(df, aes(x = t, y = inve_ebayes)) + 
  geom_line(linetype="solid", color="#F46036", size = 0.7) + 
  labs(x = "Index", y = "Voltage") +
  theme_classic() 




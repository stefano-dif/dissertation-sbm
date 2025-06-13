library(igraph)
library(tidyverse)
library(aricode)

load("nmi_data.RData") 

nmi_SBM_sim1 <- tibble(lmb = c(rep(0.9, 50), rep(0.7, 50), rep(0.5, 50), rep(0.3, 50), rep(0.1, 50)),
                       gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 5),
                       nmi = nmisC)
                         
nmi_SBM_sim2 <- tibble(lmb = c(rep(0.9, 50), rep(0.7, 50), rep(0.5, 50), rep(0.3, 50), rep(0.1, 50)),
                       gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 5),
                       nmi = nmisD)

nmi_SBM_sim3 <- tibble(lmb = c(rep(1, 50), rep(0.8, 50), rep(0.6, 50), rep(0.4, 50), rep(0.2, 50), rep(0, 50)),
                       gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 6),
                       nmi = nmisG)

nmi_SBM_sim4 <- tibble(lmb = c(rep(1, 50), rep(0.8, 50), rep(0.6, 50), rep(0.4, 50), rep(0.2, 50), rep(0, 50)),
                       gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 6),
                       nmi = nmisH)

nmi_SBM <- rbind(nmi_SBM_sim1, nmi_SBM_sim2, nmi_SBM_sim3, nmi_SBM_sim4)

nmi_SBM %>% 
  group_by(lmb) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = lmb, y = totalmean))

nmi_SBM %>% 
  group_by(gammas) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = gammas, y = totalmean))

nmi_SBM %>% 
  group_by(lmb) %>% 
  summarise(mx = max(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = lmb, y = mx))

nmi_SBM %>% 
  group_by(gammas) %>% 
  summarise(mx = max(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = gammas, y = mx))

nmi_SBM %>% 
  group_by(gammas, lmb) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot(aes(x = lmb, y = totalmean, color = gammas, group = gammas))+
  geom_point() + geom_line()+ 
  facet_wrap(~ gammas)+
  theme(legend.position = "none")+
  ylab("Mean NMI") + xlab(expression(lambda))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

nmi_DCSBM %>% 
  group_by(gammas, lmb) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot(aes(x = lmb, y = totalmean, color = gammas, group = gammas))+
  geom_point() + geom_line()+ 
  facet_wrap(~ gammas)+
  theme(legend.position = "none")+
  ylab("Mean NMI") + xlab(expression(lambda))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

nmi_DCSBM %>% 
  group_by(gammas, lmb) %>% 
  summarise(totalmean = mean(nmi)) -> adcsbm

nmi_SBM %>% 
  group_by(gammas, lmb) %>% 
  summarise(totalmean = mean(nmi)) -> asbm

ggplot(data = adcsbm, aes(x = lmb, y = totalmean, group = gammas))+
  geom_point(aes(x = lmb, y = totalmean, group = gammas), color = "blue") + 
  geom_line(data = adcsbm, aes(x = lmb, y = totalmean, group = gammas), color = "blue")+
  geom_point(data = asbm, aes(x = lmb, y = totalmean, group = gammas), color = "red") +
  geom_line(data = asbm, aes(x = lmb, y = totalmean, group = gammas), color = "red")+ 
  facet_wrap(~ gammas)+
  theme(legend.position = "none")+
  ylab("Mean NMI") + xlab(expression(lambda))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

ggplot(data = adcsbm, aes(x = gammas, y = totalmean, group = lmb))+
  geom_point(aes(x = gammas, y = totalmean, group = lmb), color = "blue") + 
  geom_line(data = adcsbm, aes(x = gammas, y = totalmean, group = lmb), color = "blue")+
  geom_point(data = asbm, aes(x = gammas, y = totalmean, group = lmb), color = "red") +
  geom_line(data = asbm, aes(x = gammas, y = totalmean, group = lmb), color = "red")+ 
  facet_wrap(~ lmb)+
  theme(legend.position = "none")+
  ylab("Mean NMI") + xlab(expression(gamma))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

ggplot(data = nmi_SBM)+
  geom_boxplot(aes(x = lmb, group = lmb, color = lmb, y = nmi))+
  facet_wrap(~ gammas)

nmi_SBM %>% 
  group_by(gammas, lmb) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot(aes(x = gammas, y = totalmean, color = lmb, group = lmb))+
  geom_point() + geom_line()+
  facet_wrap(~lmb)+
  ylab("Mean NMI") + xlab("Degree Exponent")


nmi_DCSBM_sim1 <- tibble(lmb = c(rep(0.9, 50), rep(0.7, 50), rep(0.5, 50), rep(0.3, 50), rep(0.1, 50)),
                         gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 5),
                         nmi = nmisA)  

nmi_DCSBM_sim2 <- tibble(lmb = c(rep(0.9, 50), rep(0.7, 50), rep(0.5, 50), rep(0.3, 50), rep(0.1, 50)),
                       gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 5),
                       nmi = nmisB) 

nmi_DCSBM_sim3 <- tibble(lmb = c(rep(1, 50), rep(0.8, 50), rep(0.6, 50), rep(0.4, 50), rep(0.2, 50), rep(0, 50)),
                       gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 6),
                       nmi = nmisE)

nmi_DCSBM_sim4 <- tibble(lmb = c(rep(1, 50), rep(0.8, 50), rep(0.6, 50), rep(0.4, 50), rep(0.2, 50), rep(0, 50)),
                       gammas = rep(c(rep(3, 10), rep(2.75, 10), rep(2.5, 10), rep(2.33, 10), rep(2.25, 10)), 6),
                       nmi = nmisF)

nmi_DCSBM <- rbind(nmi_DCSBM_sim1, nmi_DCSBM_sim2, nmi_DCSBM_sim3, nmi_DCSBM_sim4)

nmi_DCSBM %>% 
  group_by(lmb) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = lmb, y = totalmean))

nmi_DCSBM %>% 
  group_by(gammas) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = gammas, y = totalmean))

nmi_DCSBM %>% 
  group_by(gammas, lmb) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = gammas, y = totalmean, color = lmb))


ggplot(data = nmi_DCSBM)+
  geom_boxplot(aes(x = lmb, group = lmb, fill = lmb, y = nmi))+
  theme(legend.position = "none")+
  facet_wrap(~ gammas)+
  ylab("Mean NMI") + xlab(expression(lambda))


ggplot(data = nmi_DCSBM)+
  geom_boxplot(aes(x = lmb, group = lmb, fill = lmb, y = nmi))+
  theme(legend.position = "none")+
  facet_wrap(~ gammas)+
  ylab("NMI") + xlab(expression(lambda))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

ggplot(data = nmi_SBM)+
  geom_boxplot(aes(x = lmb, group = lmb, fill = lmb, y = nmi))+
  theme(legend.position = "none")+
  facet_wrap(~ gammas)+
  ylab("NMI") + xlab(expression(lambda))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))


ggplot(data = nmi_DCSBM)+
  geom_boxplot(aes(x = gammas, group = gammas, fill = gammas, y = nmi))+
  theme(legend.position = "none")+
  facet_wrap(~ lmb)+
  ylab("Mean NMI") + xlab(expression(gamma))

nmi_DCSBM %>% 
  group_by(lmb) %>% 
  summarise(mx = max(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = lmb, y = mx))

nmi_DCSBM %>% 
  group_by(gammas) %>% 
  summarise(mx = max(nmi)) %>% 
  ggplot()+
  geom_point(aes(x = gammas, y = mx))

nmi_DCSBM %>% 
  group_by(gammas, lmb) %>% 
  summarise(totalmean = mean(nmi)) %>% 
  ggplot(aes(x = lmb, y = totalmean, color = gammas, group = gammas))+
  geom_point() + geom_line()+ 
  facet_wrap(~ gammas)+
  ylab("Mean NMI") + xlab("Degree Exponent")


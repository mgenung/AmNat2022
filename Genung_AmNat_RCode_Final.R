library(ggplot2)
library(magrittr)
library(scales)
library(colorspace)

# PART 1 - WATERMELON ANALYSIS
#####

# Winfree et al. 2015 data, from watermelon farms observed in 2012

# replace the file names as needed
a = read.csv("watermelon_com.csv", row.names=1)
# single-visit pollen deposition data (per-visit)
z = read.csv("watermelon_fun.csv", row.names=1)
# for both data sets, sites are in columns, species in rows

# Main text: abundance-based analysis

# derived data
S = nrow(a)
N = ncol(a)
p = apply(a, 2, function(x) x/sum(x)) # relative species frequencies at each site
A = colSums(a) # total abundance at each site
t = colSums(a*z, na.rm=T) # total function provided to each site
z.bar = colSums(p*z, na.rm=T) # mean per-visit (CWM) function provided to a site

# preliminary summaries and plots
sum(a); dim(a); dim(z)
plot(A, t); abline(h=mean(t))
plot(A, z.bar); abline(h=mean(z.bar))

# population covariance function
pop.cov = function(x, y) {
  x = as.numeric(x)
  y = as.numeric(y)
  
  mean(x*y) - mean(x)*mean(y)}

# goal of abundance-based Price equation analysis: 
# partition difference between each site's function and grand means:
E.t = mean(t)
E.A = mean(A)
E.z.bar = mean(z.bar)

# Equation 5 test plot (not in paper, should be a perfect fit)
plot(t - E.t,

    (A - E.A)*E.z.bar +
    A*(z.bar - E.z.bar) -
    pop.cov(A, z.bar)
 ); abline(0,1)

# we want to study how expected total function changes with overall site abundance, simply because, unlike the unconditional mean expected function, the conditional mean can be imagined as function gain/loss per individual.
abund.mod = lm((A - E.A)*E.z.bar + E.t ~ A)
cwm.mod = lm(A*(z.bar - E.z.bar) + E.t ~ A)
t.mod = lm(t ~ A)

###
# numeric summaries
# check: the slopes of the two partitions sum exactly to the slope of the overall abundance-function trend.
coef(abund.mod)[2] + coef(cwm.mod)[2]; coef(t.mod)[2]
# correlation between abundance and CWM function
round(cor(A, z.bar),2)
# population coveriance term in equation 5 is what percent of grand mean function
round(100*pop.cov(A, z.bar)/E.t, 2)
# total function partitions v. abundance.
sapply(list(t.mod, abund.mod, cwm.mod), 
       function(m) {
         c(coef(m)[2], confint(m, 'A', level=0.95))})
# warning about perfect fit is expected

# setting up first R-generated figure from published manuscript
mytheme <- theme(panel.background = element_rect(fill = 'white'), 
                 axis.title.x = element_text(vjust = -0.5, size = 22), 
                 axis.text.x = element_text(color = "black", size = 16),
                 axis.title.y = element_text(vjust = 0.75, size = 22), 
                 axis.text.y = element_text(color = "black", size = 16),
                 plot.title = element_text(face="bold", size=20),
                 axis.line.x = element_line(color="black", size = 1.2),
                 axis.line.y = element_line(color="black", size = 1.2))

# Fig2a
#####
get = 16 # select one site to illustrate
Fig2a = ggplot()+
  
  # dashed lines
  geom_hline(yintercept = mean(z.bar), lty=2, col="blue", size=1.3, alpha=0.7) +
  geom_vline(xintercept = mean(A), lty=2, col="red", size=1.3, alpha=0.7) +
  
  # arrows
  geom_segment(aes(x = A[get]-2, y = z.bar[get], xend = mean(A), yend = z.bar[get]),
               col="red", size=1.2, arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = A[get], y = z.bar[get]-2, xend = A[get], yend = mean(z.bar)),
               col="blue", size=1.2, arrow = arrow(length = unit(0.5, "cm"))) +
  
  # gray points
  geom_point(aes(x=A, y=z.bar), 
             alpha=0.7, stroke=1, shape=21, col="black", fill="gray", size=t/1200) +
  
  # mean point
  geom_point(aes(x=mean(A), y=mean(z.bar)), 
             stroke=2, shape=21, col="black", fill="white", size=6) +
  
  scale_x_continuous("Site Abundance (A)") +
  scale_y_continuous(expression(paste("per-cap function (", bar("z"), ")"))) +
  annotate("text", x=0, y=170, label= "A", size=14) + 
  mytheme
#####

# Fig2b
#####
Fig2b = ggplot()+
  
  ## gray
  # points
  geom_point(aes(x=A, y=t), 
             stroke=1, shape=21, col="black", fill="gray", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=A, y=t), stat="smooth", 
            col="gray", alpha=0.9, size=2, method="lm", se=F) +
  
  ## red
  # points
  geom_point(aes(x=A, y=(A - E.A)*E.z.bar + E.t), 
             stroke=1, shape=21, col="black", fill="red", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=A, y=(A - E.A)*E.z.bar + E.t), stat="smooth",
            col="red", alpha=0.9, size=2, method="lm", se=F) +
  
  ## blue
  # points
  geom_point(aes(x=A, y=A*(z.bar - E.z.bar) + E.t), 
             stroke=1, shape=21, col="black", fill="blue", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=A, y=A*(z.bar - E.z.bar) + E.t), stat="smooth",
            col="blue", alpha=0.9, size=2, method="lm", se=F) +
  
  scale_x_continuous("Site Abundance (A)") +
  scale_y_continuous("Total Function (T)") +
  annotate("text", x=15, y=18000, label= "B", size=14) + 
  mytheme
#####

# Fig2c
#####
Fig2c = ggplot()+
  
  ## gray
  # points
  geom_point(aes(x=z.bar, y=t), 
             stroke=1, shape=21, col="black", fill="gray", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=z.bar, y=t), stat="smooth", 
            col="gray", alpha=0.9, size=2, method="lm", se=F) +
  
  ## red
  # points
  geom_point(aes(x=z.bar, y=(A - E.A)*E.z.bar + E.t), 
             stroke=1, shape=21, col="black", fill="red", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=z.bar, y=(A - E.A)*E.z.bar + E.t), stat="smooth",
            col="red", alpha=0.9, size=2, method="lm", se=F) +
  
  ## blue
  # points
  geom_point(aes(x=z.bar, y=A*(z.bar - E.z.bar) + E.t), 
             stroke=1, shape=21, col="black", fill="blue", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=z.bar, y=A*(z.bar - E.z.bar) + E.t), stat="smooth",
            col="blue", alpha=0.9, size=2, method="lm", se=F) +
  
  scale_x_continuous(expression(paste("Per-capita Function (", bar("z"), ")"))) +
  scale_y_continuous("Total Function (T)") +
  annotate("text", x=75, y=17600, label= "C", size=14) + 
  mytheme
#####

Fig2pdf = list(Fig2a, Fig2b, Fig2c)
ggsave(
  filename = "/Users/Mark/Downloads/AmNat_Fig2_pdf.pdf", 
  plot = marrangeGrob(Fig2pdf, nrow=1, ncol=3), 
  width = 20, height = 7
)

# Fig A1 uses the same data in a different way
# set up Fig A1 - richness based analysis
#####
# goal of richness-based Price equation analysis: 
# partition difference between each site's function and grand means

# new terms
x = a*z
b = ifelse(a>0, 1, 0)
w = as.data.frame.matrix(apply(b, 2, function(x) x/sum(x)))
colSums(w); dim(w) # nine genera (rows), sixteen sites (columns)
s = colSums(b)
bar.x = colSums(x*w) # mean function per species
t = s*bar.x # total function

# get means across sites
E.s = mean(s)
E.bar.x = mean(bar.x)
E.w = rowMeans(w) # for each sp
E.x = rowMeans(x) # for each sp
E.t = mean(t)

# partition terms for each site
rich = (s - E.s)*E.bar.x
comp = colSums(x*(apply(w, 1, function(x) {x - mean(x)}) %>% t))
perc = colSums(rowMeans(w)*(apply(x, 1, function(x) {x - mean(x)}) %>% t))
# NOTE perc is calculated relative to E(p) not p
# covariances
sp.cov = sum(sapply(1:nrow(a), function(i) pop.cov(c(w[i,]), c(x[i,]))))
site.cov = pop.cov(s, bar.x)
#####

# Fig A1a
#####
get = 15 # select one site to illustrate
ggplot()+
  
  # dashed lines
  geom_hline(yintercept = mean(bar.x), lty=2, col="blue", size=1.3, alpha=0.7) +
  geom_vline(xintercept = mean(s), lty=2, col="red", size=1.3, alpha=0.7) +
  
  # arrows
  geom_segment(aes(x = s[get], y = bar.x[get], xend = mean(s), yend = bar.x[get]),
               col="red", size=1.2, arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = s[get], y = bar.x[get]-2, xend = s[get], yend = mean(bar.x)),
               col="blue", size=1.2, arrow = arrow(length = unit(0.5, "cm"))) +
  
  # gray points
  geom_point(aes(x=s, y=bar.x), 
             alpha=0.7, stroke=1, shape=21, col="black", fill="gray", size=6) +
  
  # mean point
  geom_point(aes(x=mean(s), y=mean(bar.x)), 
             stroke=2, shape=21, col="black", fill="white", size=6) +
  
  scale_x_continuous("Site Richness (s)") +
  scale_y_continuous(expression(paste("Mean Per-Species Function (", bar("x"), ")"))) +
  mytheme
#####

# Fig A1b
#####
ggplot()+
  
  ## gray
  # points
  geom_point(aes(x=s, y=t), 
             stroke=1, shape=21, col="black", fill="gray", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=s, y=t), stat="smooth", 
            col="gray", alpha=0.9, size=2, method="lm", se=F) +
  
  ## red
  # points
  geom_point(aes(x=s, y=E.t + rich), 
             stroke=1, shape=21, col="black", fill="red", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=s, y=E.t + rich), stat="smooth",
            col="red", alpha=0.9, size=2, method="lm", se=F) +
  
  ## blue
  # points
  geom_point(aes(x=s, y=E.t + s*comp), 
             stroke=1, shape=21, col="black", fill="blue", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=s, y=E.t + s*comp), stat="smooth",
            col="blue", alpha=0.9, size=2, method="lm", se=F) +
  
  ## green
  # points
  geom_point(aes(x=s, y=E.t + s*perc), 
             stroke=1, shape=21, col="black", fill="darkgreen", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=s, y=E.t + s*perc), stat="smooth",
            col="darkgreen", alpha=0.9, size=2, method="lm", se=F) +
  
  scale_x_continuous("Site Richness (s)") +
  scale_y_continuous("Total Function (T)") +
  mytheme
#####

# Fig A1c
#####
ggplot()+
  
  ## gray
  # points
  geom_point(aes(x=bar.x, y=t), 
             stroke=1, shape=21, col="black", fill="gray", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=bar.x, y=t), stat="smooth", 
            col="gray", alpha=0.9, size=2, method="lm", se=F) +
  
  ## red
  # points
  geom_point(aes(x=bar.x, y=E.t + rich), 
             stroke=1, shape=21, col="black", fill="red", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=bar.x, y=E.t + rich), stat="smooth",
            col="red", alpha=0.9, size=2, method="lm", se=F) +
  
  ## blue
  # points
  geom_point(aes(x=bar.x, y=E.t + s*comp), 
             stroke=1, shape=21, col="black", fill="blue", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=bar.x, y=E.t + s*comp), stat="smooth",
            col="blue", alpha=0.9, size=2, method="lm", se=F) +
  
  ## green
  # points
  geom_point(aes(x=bar.x, y=E.t + s*perc), 
             stroke=1, shape=21, col="black", fill="darkgreen", alpha=0.5, size=6) +
  # line
  geom_line(aes(x=bar.x, y=E.t + s*perc), stat="smooth",
            col="darkgreen", alpha=0.9, size=2, method="lm", se=F) +
  
  scale_x_continuous(expression(paste("Mean Per-Species Function (", bar("x"), ")"))) +
  scale_y_continuous("Total Function (T)") +
  mytheme
#####

# set up Figure 3
#####
# load in stream invert data
comdat = read.csv("stream_inverts.csv", row.names=1)
sitedat = read.csv("stream_data.csv", row.names=1)
fdat = read.csv("stream_invert_size.csv", row.names=1)

# define Price equation variables for the new analysis
# to maintain clear links with the text, this overwrites some previous variables
# if you run the code in order everything will be fine
a = comdat
p = (a/rowSums(a))
z = fdat
dim(a); dim(p); dim(z)
sum(is.na(z)); prod(dim(z))

A = rowSums(a)
t = rowSums(a*z, na.rm=T)
e = sitedat$dist

# run model
mod.A = lm(A ~ e)
mod.p = apply(p, 2, function(y) lm(y ~ e))
mod.z = apply(z, 2, function(y) lm(y ~ e))
#####

# Fig 3
# stream inverts - make sure to run stream_invert_analysis.R first

T.delta.A = predict(mod.A)*sum(sapply(mod.p, coef)[1,]*sapply(mod.z, coef)[1,])
T.delta.p = coef(mod.A)[1]*rowSums(sapply(mod.p, predict)*
                                     matrix(sapply(mod.z, coef)[1,], nrow=nrow(p), ncol=6, byrow=T))
T.delta.z = coef(mod.A)[1]*rowSums(sapply(mod.z, predict, newdata = data.frame(e=e))*matrix(sapply(mod.p, coef)[1,], nrow=nrow(p), ncol=6, byrow=T))

# New Fig3a
#####
z.bar = A/t
e.std = e - min(e) + 3

f3.colors = rainbow_hcl(6)

az = a*z
az[is.na(az)] = 0

Fig3a = ggplot()+
  
  # All 
  geom_point(aes(x=e, y=rowSums(az)), 
             alpha=0.3, stroke=1, shape=21, col="black", fill="black", size=6) +
  geom_line( aes(x=e, y=rowSums(az)), stat="smooth",
             alpha=0.9, size=2, col="black", method="lm", se=F) +
  
  scale_x_continuous("Pollution (e)") +
  scale_y_continuous("Total Biomass (T)") +
  annotate("text", x=-3.5, y=1, label= "A", size=14) + 
  mytheme
#####

# New Fig 3b
#####
Fig3b = ggplot()+
  
  # Fn group 1
  geom_point(aes(x=e, y=az[,1]), 
             alpha=0.6, stroke=1, shape=21, col="black", fill=f3.colors[1], size=6) +
  geom_line( aes(x=e, y=az[,1]), stat="smooth",
             alpha=0.8, size=1.5, col=f3.colors[1], method="lm", se=F) +
  
  # Fn group 2
  geom_point(aes(x=e, y=az[,2]), 
             alpha=0.6, stroke=1, shape=21, col="black", fill=f3.colors[2], size=6) +
  geom_line( aes(x=e, y=az[,2]), stat="smooth",
             alpha=0.8, size=1.5, col=f3.colors[2], method="lm", se=F) +
  
  # Fn group 3
  geom_point(aes(x=e, y=az[,3]), 
             alpha=0.6, stroke=1, shape=21, col="black", fill=f3.colors[3], size=6) +
  geom_line( aes(x=e, y=az[,3]), stat="smooth",
             alpha=0.8, size=1.5, col=f3.colors[3], method="lm", se=F) +
  
  # Fn group 4
  geom_point(aes(x=e, y=az[,4]), 
             alpha=0.6, stroke=1, shape=21, col="black", fill=f3.colors[4], size=6) +
  geom_line( aes(x=e, y=az[,4]), stat="smooth",
             alpha=0.8, size=1.5, col=f3.colors[4], method="lm", se=F) +
  
  # Fn group 5
  geom_point(aes(x=e, y=az[,5]), 
             alpha=0.6, stroke=1, shape=21, col="black", fill=f3.colors[5], size=6) +
  geom_line( aes(x=e, y=az[,5]), stat="smooth",
             alpha=0.8, size=1.5, col=f3.colors[5], method="lm", se=F) +
  
  # Fn group 6
  geom_point(aes(x=e, y=az[,6]), 
             alpha=0.6, stroke=1, shape=21, col="black", fill=f3.colors[6], size=6) +
  geom_line( aes(x=e, y=az[,6]), stat="smooth",
             alpha=0.8, size=1.5, col=f3.colors[6], method="lm", se=F) +
  
  scale_x_continuous("Pollution (e)") +
  scale_y_continuous("Functional Group\nBiomass", limits=c(0,0.5)) +
  annotate("text", x=-3.5, y=0.5, label= "B", size=14) + 
  mytheme
#####

# New Fig3c
#####
Fig3c = ggplot()+
  
  ## gray
  # points
  geom_point(aes(x=e, y=t), 
             stroke=1, shape=21, col="black", fill="black", alpha=0.3, size=6) +
  # line
  geom_line(aes(x=e, y=t), stat="smooth", 
            col="black", alpha=0.9, size=2, method="lm", se=F) +
  
  ## red
  # line
  geom_line(aes(x=e[order(e)], y=T.delta.A[order(e)]), stat="smooth",
            col="red", alpha=0.85, size=2, method="lm", se=F) +
  
  ## blue
  # line
  geom_line(aes(x=e[order(e)], y=T.delta.p[order(e)]), stat="smooth",
            col="blue", alpha=0.85, size=2, method="lm", se=F) +
  
  ## green
  # line
  geom_line(aes(x=e[order(e)], y=T.delta.z[order(e)]), stat="smooth",
            col="darkgreen", alpha=0.85, size=2, method="lm", se=F) +
  
  scale_x_continuous(expression(paste("Pollution (e)"))) +
  scale_y_continuous("Total Biomass (T)") +
  annotate("text", x=-3.5, y=1, label= "C", size=14) + 
  mytheme
#####

Fig3pdf = list(Fig3a, Fig3b, Fig3c)
ggsave(
  filename = "/Users/Mark/Downloads/AmNat_Fig3_pdf.pdf", 
  plot = marrangeGrob(Fig3pdf, nrow=1, ncol=3), 
  width = 20, height = 7
)

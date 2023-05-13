

###### check gamma distribution ######
 
x11() 
shape = 1
scale = 1/2
sims.1<- rgamma(100, shape =shape,   scale = scale)
p.try1<- ggplot(data.frame(x = sims.1), aes(x=x)) + geom_density() +
ggtitle(paste("gamma)", shape, " scale ", scale,  
               'mean', mean(sims.1), sep = ""))+ xlim(c(0, 3))

p.try1

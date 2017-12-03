
# source our functions
source('select.R')

# Get data set
dat<-mtcars

# specify lm or glm
mod <- lm 

# use AIC function in stats package because it matches the step() functions
library('stats')
f<-function(fit,...){return(extractAIC(fit,...)[2])} 

# Run GA 
num_gens=40
n=30
results<-genetic(mod,dat,f,num_gens=num_gens,n=n)

# inspect fittest model 
summary(results$fittest_model)

# plot results 
X<-matrix(rep(seq(num_gens),nrow(results$fitness),num_gens),nrow=n,ncol=num_gens)
Y<-t(results$fitness)
plot(X,Y,xlab='generation',ylab='AIC',pch=19,cex=0.5)

lines(seq(num_gens),apply(results$fitness,FUN=min,MARGIN=2),lty=1,col='green')


# Comparing against forward selection 
null=lm(dat[,1]~1, data=dat[,-1])
full=lm(dat[,1]~., data=dat[,-1])
step(null, scope=list(lower=null, upper=full), direction="forward")

# Comparing against backwards selection
step(full, direction="backward")


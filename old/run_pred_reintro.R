rep = 1000
dist_mean = 2

source("pred_reintro_functions.R")
pred_reintro_general(dist_mean,rep)

fname = paste('results_dmean_',toString(dist_mean),'_rep_',toString(rep),'.rds',sep='')
load(fname)


x<- read.csv("interactions.csv")  # x is now a data.frame
sppNames <- colnames(x)


# Added this to adjust marjinsback to default after as it's reset below
par(mar=c(5.1, 4,4.1 ,2.1))

#res_num = ceiling(runif(1,max=rep))
res_num = sample( 1:rep, 1)

n_s = dim(y_array)[2]
y1 = y_array[,2,res_num] # not sure why you left the first one off I've added it back in
y1 = y1 /y1[1]
Years = t_vec

plot(Years,y1, ylim=c(0, max(c(1,outcomes_abund[,res_num]))),type="l",ylab = "Relative abundance", col=1)
for (i in 3:n_s){
  y = y_array[,i,res_num]
  y = y /y[1]
  # specify the colour using the col argument, each colour can be specified by a number
  lines(t_vec,y, col=i)
  print(y[length(y)])
}

# add a legend
vec.of.text.for.legend <- sppNames[-2]
legend('topleft', vec.of.text.for.legend, col=1:n_s, lty=1 )



# Chris, below is an example of how to do this without a for loop using the matplot()
# probably neater to delete the above and use this instead
traj.to.plot <- y_array[,,res_num]

# now do the normalization - this was a but tricky to work out but you can do it
# in one line as follows (e.g. see http://stackoverflow.com/questions/9447801/dividing-columns-by-colsums-in-r)
norm.traj.to.plot <- t(t(traj.to.plot)/traj.to.plot[1,])

norm.traj.to.plot <- norm.traj.to.plot [,-1] # remove the first col as per above

no.spp <- dim(norm.traj.to.plot)[2]

# Make the plot of all trajectories
matplot(norm.traj.to.plot, ylab = "Relative abundance", type = 'l', lty=1, 
	col=c(1:no.spp))#, main='new version with matplot()')
# add a legend
legend('topleft', vec.of.text.for.legend, col=1:no.spp, lty=1 )


# bar plots

# add this so the lables will fit in
par(mar=c(5.1, 8,4.1 ,2.1))

mean_outcomes = rowMeans(outcomes)
barplot(mean_outcomes[length(mean_outcomes):1],horiz=TRUE,xlim=c(0,1),
	names.arg=sppNames[length(sppNames):3],
	las=1,
	main = "Equilibrium outcomes",xlab="Frequency of increase")


mean_outcomes = rowMeans(outcomes_dir)
barplot(mean_outcomes[length(mean_outcomes):1],horiz=TRUE,xlim=c(0,1),
	names.arg=sppNames[length(sppNames):3],las=1,
	main = "Initial changes",xlab="Frequency of increase")



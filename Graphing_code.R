# Code for generating the figures in chapters 3,4 and 5

# Here we utilise several functions from the igraph package the documentation of
# which can be found here:https://igraph.org/r/doc/

# Creating, simulating and plotting the modified beta model
library(igraph) # importing igraph for the relevant functions 

# creating a function for the augmented beta model

Augmented_beta <- function(L, k, phi){
  g_base <- sample_smallworld(dim = 1, size = L, nei = k, p = 0)
  # use the sample_smallworld function here as we are not modifying 
  # the lattice as p = 0 
  
  shortcuts <- round(L*k*phi)
  edge_shortcuts <- 2*shortcuts
  # we need this to be a whole number so will have to settle for rounding it
  # edge_shortcuts is twice as long as we will create a vertex list
  
  edges <- sample(1:L, size = edge_shortcuts, replace = TRUE)
  g_new <- make_empty_graph(n = L)
  g_new <- add_edges(g_new,edges)
  g_new <- as_undirected(g_new)
  
  adj_base <- as_adjacency_matrix(g_base)
  adj_new <- as_adjacency_matrix(g_new)
  adj <- adj_base + adj_new
  g_final <- graph_from_adjacency_matrix(adj, mode = 'undirected')
  
  return(g_final)
}


# Illustrating small worlds behaviour of the augmented beta model
phi <- seq(-5,1,length.out = 50)
phi <- 10**(phi)
plot(phi)
length <- length(phi)
L <- 2000
k <- 5
dim <- 1
sims <- 100

char_paths <- rep(0,length)
clust <- rep(0,length)
set.seed(1)
for (i in 1:length) {
  temp_chr <- rep(0,sims)
  temp_clust <- rep(0,sims)
  for (j in 1:sims){
    g <- Augmented_beta(L,k,phi[i])
    temp_chr[j] <- mean_distance(g, unconnected = FALSE)
    temp_clust[j] <- mean(transitivity(g, type = 'local'))
  }
  char_paths[i] <-mean(temp_chr)
  clust[i] <- mean(temp_clust)
}
char_paths
clust
plot(seq(-5,1,length.out = 50),clust/max(clust), ylab = 'Normalised char and clust'
     , xlab = 'Phi in terms of exponents on 10',
     ylim = c(0,1)) 
points(seq(-5,1,length.out = 50),char_paths/max(char_paths), col = 'orange')
lines(seq(-5,1,length.out = 50),clust/max(clust))
lines(seq(-5,1,length.out = 50),char_paths/max(char_paths), col = 'orange')

# note that there are issues to be discussed here namely loops and the way that 
# the shortcuts are selected as we are using large L we have not mitigated 
# further against this


# Here we create a function to show the percolation threshold in the node percolation only case.
library(polynom) # this is used to solve the polynomial

perculation_threshold <- function(k, phi){
  m <- (1 + 1/(2*k*phi))
  
  # we will have to add an adaptation for when k is 1 and for simplicity 
  # when k = 2
  
  if (k ==1) {
    z_poly <- polynomial(coef = c(2,(-m-2),1))
    z <- solve(z_poly)
  }
  
  if (k ==2) {
    z_poly <- polynomial(coef = c(2,-2,-m,1))
    z <- solve(z_poly)
  }
  
  if (k >2) {
    z_poly <- polynomial(coef = c(1,-m,rep(0,k-2),-2,2))
    z_poly <- polynomial(coef = c(2,-2,rep(0,k-2),-m,1))
    z <- solve(z_poly)
  }
  
  
  pc <- 1-z
  # this gives us multiple roots but we only want the roots where pc is a valid 
  # probability can we guarantee that there will only be one of these 
  #print(z_poly)
  return(pc)
}

# we will use this function to find the percolation thresholds.
# we will vary the suceptibility probability which we use to eliminate nodes 
# from the graph with this probability in this way we create our node 
# percolation graphs. A similar strategy can be used in the bond percolation 
# case and the node and bond percolation case.


# need to discuss how we are manually selecting the correct solution as we know 
# we need a probability so require solutions which are real between 0 and 1.

p <- seq(0,1,length.out = 30)
#p <- 10**(p)
plot(p)
length <- length(p)
L <- 10000
k <- 2
dim <- 1
phi <- 0.01

proportion_in_largest <- list()

pc <- perculation_threshold(k,phi)
set.seed(1)
g <- Augmented_beta(L,k,phi)
################# FIX
for (i in 1:length) {
  to_sample <- floor(L*p[i])
  susceptible <- sample(1:L,replace = FALSE,size = to_sample)
  sub_susc <- subgraph(g,susceptible)
  largest_comp <- largest_component(sub_susc)
  v_in <- vcount(largest_comp) 
  prop_in_largest <- v_in/L
  
  proportion_in_largest[i] <- prop_in_largest
}

proportion_in_largest

# here I manually ran this function for k=1, k=2 and k=5 and saved my results in
# the following sections so that I could then plot each line on the same graph 
# recreating Watt's results 

#prop_large_k_1 <- proportion_in_largest
#prop_large_k_2 <- proportion_in_largest
#prop_large_k_5 <- proportion_in_largest

#pc_k_1 <- 0.9622372
#pc_k_2 <- 0.7574794
#pc_k_5 <- 0.4010075

# these critical probs were manually selected from the percolation threshold 
# results 

critical_probs <- c(pc_k_1, pc_k_2, pc_k_5)

proportion_in_largest <- list()
to_sample <- floor(L*pc_k_2) # we need to ensure the relevent pc is used here
susceptible <- sample(1:L,replace = FALSE,size = to_sample)
sub_susc <- subgraph(g,susceptible)
largest_comp <- largest_component(sub_susc)
v_in <- vcount(largest_comp) 
prop_in_largest <- v_in/L

#prop_pc_k_1 <- prop_in_largest
#prop_pc_k_2 <- prop_in_largest
#prop_pc_k_5 <- prop_in_largest
c(prop_pc_k_1,prop_pc_k_2,prop_pc_k_5)

# can add a separate point for each critical property

plot(p,prop_large_k_1,ylab = 'Proportion of nodes in GCC',
     xlab = 'Probability of susceptibility') 
lines(p,prop_large_k_1)
points(p,prop_large_k_2, pch = 2)
lines(p,prop_large_k_2)
points(p,prop_large_k_5, pch = 5)
lines(p,prop_large_k_5)
points(critical_probs,c(prop_pc_k_1,prop_pc_k_2,prop_pc_k_5), pch = c(1,2,5), col = 'red')


# the red is our critical susceptibility probabilities and the respective 
# proportion of nodes in the GCC

# to check that my functions are correct I will try to recreate the graphs in 
# the paper

# I want to vary the shortcut density phi and then plot the percolation 
# threshold for different values of k and then overlay them

# I will let L = 1000000
# k = 1,2,5
# phi 13 evenly spaced between 10^-4 and 1

# as we are calculating the percolation threshold we don't need to iterate here
phi <- seq(-4,1,length.out = 13)
phi <- 10**(phi)
plot(phi)
length <- length(phi)
L <- 1000000
k <- 2
dim <- 1

percs <- list()

for (i in 1:length) {
  percs[[i]] <- perculation_threshold(k,phi[i])
}
# here we are selecting the relevant values by hand, where we want the values 
# between 0 and 1 which are real alternatively we can take the kth value
percs
percs_k_2 <- c(0.9721182, 0.9553579, 0.9289351, 0.8879597, 0.8261762, 0.7373256, 
               0.6189345, 0.4779169, 0.3325478, 0.2052267, 0.1111962, 0.05303857, 0.02283967)
percs_k_1 <- c(0.9996002, 0.9989576, 0.9972859, 0.9929618, 0.9819338, 0.9547767,
               0.8930736, 0.7735294, 0.592981, 0.3901031, 0.2184846, 0.1055813,
               0.04563561)

percs_k_5 <- c(0.7291632, 0.6767954, 0.6159262, 0.5462249, 
               0.4681901, 0.3837737, 0.2969945, 0.2139208, 0.1413245, 0.08437193,
               0.04491957, 0.02127259, 0.009141008)

plot(seq(-4,1,length.out = 13),percs_k_1,ylab = 'Percolation threshold',
     xlab = 'Phi in terms of exponents on 10') 
lines(seq(-4,1,length.out = 13),percs_k_1)
points(seq(-4,1,length.out = 13),percs_k_2, pch = 2)
lines(seq(-4,1,length.out = 13),percs_k_2)
points(seq(-4,1,length.out = 13),percs_k_5, pch = 5)
lines(seq(-4,1,length.out = 13),percs_k_5)


# consider finding a critical phi in terms of the percolation of the population
shortcut_threshold <- function(k,p){
  phi <- (1-p)**k/(2*k*p*(2-(1-p)**k))
  return(phi)
}
# again note how this is not dependent on L - except that we take L to be large
# I am going to plot the corresponding graphs for each of the degrees above

p <- seq(0.0001,1,length.out = 16)
# most of these are high in the critical case so I'll start with a range from 
# 0 to 1 about 16 long
plot(p)
length <- length(p)
L <- 1000000
k <- 5


phis <- list()

for (i in 1:length) {
  phis[[i]] <- shortcut_threshold(k,p[i])
}

phis
#phis_k_1 <- phis
#phis_k_2 <- phis
#phis_k_5 <- phis

x <- seq(0.0001,1,length.out = 16)

plot(x,phis_k_1,ylab = 'Shortcut threshold',
     xlab = 'Susceptibility probability') 
lines(x,phis_k_1)
points(x,phis_k_2, pch = 2)
lines(x,phis_k_2)
points(x,phis_k_5, pch = 5)
lines(x,phis_k_5)


# these all have a massive jump between the first two points so lets narrow 
# this down as I did before 


p <- seq(-4,0,length.out = 16)
p <- 10**(p)
plot(p)
# most of these are high in the critical case so I'll start with a range from 
# 0 to 1 about 16 long
plot(p)
length <- length(p)
L <- 1000000
k <- 1


phis <- list()

for (i in 1:length) {
  phis[[i]] <- shortcut_threshold(k,p[i])
}

phis
#phis_k_1 <- phis
#phis_k_2 <- phis
#phis_k_5 <- phis

x <- p

plot(x,phis_k_1,ylab = 'Shortcut threshold',
     xlab = 'Susceptibility probability in terms of exponents on 10') 
lines(x,phis_k_1)
points(x,phis_k_2, pch = 2)
lines(x,phis_k_2)
points(x,phis_k_5, pch = 5)
lines(x,phis_k_5)

# here I am going to reduce the axis so that I can better see whats going on 

# also I need to consider the bounds on the shortcuts e.g. what is reasonable 
# perhaps I will scale this relative to L or even the total base edges in the 
# network so this is also dependent on k


p <- seq(-4,-2,length.out = 16)
p <- 10**(p)
plot(p)
# most of these are high in the critical case so I'll start with a range from 
# 0 to 1 about 16 long
plot(p)
length <- length(p)
L <- 1000000
k <- 5


phis <- list()

for (i in 1:length) {
  phis[[i]] <- shortcut_threshold(k,p[i])
}

phis
#phis_k_1 <- phis
#phis_k_2 <- phis
#phis_k_5 <- phis

x <- seq(-4,-2,length.out = 16)

plot(x,phis_k_1,ylab = 'Shortcut threshold',
     xlab = 'Susceptibility probability in terms of exponents on 10') 
lines(x,phis_k_1)
points(x,phis_k_2, pch = 2)
lines(x,phis_k_2)
points(x,phis_k_5, pch = 5)
lines(x,phis_k_5)


# all this behavior is occurring in the region I had issues with when justifying 
# the use of this model

# lets look at a region where the behavior is more sensible

p <- seq(0.5,1,length.out = 16)
#p <- 10**(p)
plot(p)
# most of these are high in the critical case so I'll start with a range from 
# 0 to 1 about 16 long
plot(p)
length <- length(p)
L <- 1000000
k <- 5


phis <- list()

for (i in 1:length) {
  phis[[i]] <- shortcut_threshold(k,p[i])
}

phis
#phis_k_1 <- phis
#phis_k_2 <- phis
#phis_k_5 <- phis

x <- seq(0.5,1,length.out = 16)

plot(x,phis_k_1,ylab = 'Shortcut threshold',
     xlab = 'Susceptibility probability '
)
lines(x,phis_k_1)
points(x,phis_k_2, pch = 2)
lines(x,phis_k_2)
points(x,phis_k_5, pch = 5)
lines(x,phis_k_5)

# This gives the graph for the shortcut threshold for larger susceptibility probabilities
# now I will plot the critical shortcut probability

p <- seq(-4,0,length.out = 16)
p <- 10**(p)
plot(p)
# most of these are high in the critical case so I'll start with a range from 
# 0 to 1 about 16 long
plot(p)
length <- length(p)
L <- 1000000
k <- 1


psi <- list()

for (i in 1:length) {
  psi[[i]] <- 2*k*shortcut_threshold(k,p[i])/L
  
}

psi
#psi_k_1 <- psi
#psi_k_2 <- psi
#psi_k_5 <- psi

x <- seq(-4,0,length.out = 16)

plot(x,psi_k_1,ylab = 'Critical Shortcut Probability',
     xlab = 'Susceptibility probability in terms of exponents on 10'
)
lines(x,psi_k_1)
points(x,psi_k_2, pch = 2)
lines(x,psi_k_2)
points(x,psi_k_5, pch = 5)
lines(x,psi_k_5)


# these basically look the same which makes sense for small p


##############################################################################
# Moving on to the bond percolation case
# first creating a graph for the percolation threshold

bond_percolation_threshold <- function(k, phi){
  m <- (1 + 1/(2*k*phi))
  # we will have two cases k=1 and k=2
  
  if (k ==1) {
    z_poly <- polynomial(coef = c(2,(-m-2),1))
    z <- solve(z_poly)
    pc <- 1-z
  }
  
  if (k ==2) {
    z_poly <- polynomial((coef = c(1, (-4-4*phi), 7, (-7-12*phi), (4 + 12*phi), (-1+8*phi), (-20*phi), (8*phi))))
    pc <- solve(z_poly)[4]
  }
  return(pc)
}



phi <- seq(-4,1,length.out = 13)
phi <- 10**(phi)
plot(phi)
length <- length(phi)
L <- 1000000
k <- 2
dim <- 1

percs <- list()

for (i in 1:length) {
  percs[[i]] <- bond_percolation_threshold(k,phi[i])
}
# here we are selecting the relevant values by hand, where we want the values 
# between 0 and 1 which are real alternatively we can take the kth value
percs
#percs_k_2 <- percs
percs_k_1 <- c(0.9996002, 0.9989576, 0.9972859, 0.9929618, 0.9819338, 0.9547767,
               0.8930736, 0.7735294, 0.592981, 0.3901031, 0.2184846, 0.1055813,
               0.04563561)

plot(seq(-4,1,length.out = 13),percs_k_1,ylab = 'Percolation threshold',
     xlab = 'Phi in terms of exponents on 10') 
lines(seq(-4,1,length.out = 13),percs_k_1)
points(seq(-4,1,length.out = 13),percs_k_2, pch = 2)
lines(seq(-4,1,length.out = 13),percs_k_2)



# consider finding a critical phi in terms of the percolation of the population
shortcut_threshold_bond <- function(k,p){
  if (k==1){
    phi <- (1-p)**k/(2*k*p*(2-(1-p)**k)) 
  }
  
  if (k == 2){
    phi <- (1-p)**3*(1-p+p**2)/(4*p*(1+3*p**2 -3*p**3 - 2*p**4 + 5*p**5 -2*p**6))
  }
  return(phi)
}

p <- seq(-4,-2,length.out = 16)
p <- 10**(p)
#p <- seq(0.5,1,length.out = 16)

plot(p)
# most of these are high in the critical case so I'll start with a range from 
# 0 to 1 about 16 long
plot(p)
length <- length(p)
L <- 1000000
k <- 1


phis <- list()

for (i in 1:length) {
  phis[[i]] <- shortcut_threshold_bond(k,p[i])
}

phis
#phis_k_1 <- phis
#phis_k_2 <- phis


x <- seq(-4,-2,length.out = 16)

plot(x,phis_k_1,ylab = 'Shortcut threshold',
     xlab = 'Transmission probability '
)
lines(x,phis_k_1)
points(x,phis_k_2, pch = 2)
lines(x,phis_k_2)



# Finding the amount in the GCC
p <- seq(0,1,length.out = 30)
plot(p)
length <- length(p)
L <- 10000
k <- 2
dim <- 1
phi <- 0.01

pc_temp <- bond_percolation_threshold(k,phi)
#pc_k_1 <- pc_temp[1] - manually selecting the correct solution here (valid probability)
#pc_k_2 <- Re(pc_temp[1]) - removing the 0 i component

proportion_in_largest <- list()
set.seed(1)
for (i in 1:length) {
  g <- Augmented_beta(L,k,phi)
  edges <- as_edgelist(g)
  n_edge <- ecount(g)
  to_sample <- floor(n_edge*p[i])
  present <- sample(1:n_edge,replace = FALSE,size = to_sample)
  ids <- get_edge_ids(g, as.vector(t(edges[present,])))
  sub_susc <- subgraph_from_edges(g,ids)
  largest_comp <- largest_component(sub_susc)
  v_in <- vcount(largest_comp) 
  prop_in_largest <- v_in/L
  
  proportion_in_largest[i] <- prop_in_largest
}

proportion_in_largest

# here I manually ran this function for k=1 and k=2 and saved my results in
# the following sections so that I could then plot each line on the same graph 
# recreating his results 

#prop_large_k_1 <- proportion_in_largest
#prop_large_k_2 <- proportion_in_largest

# these critical probs were manually selected from the percolation threshold 
# results 

#critical_probs <- c(pc_k_1, pc_k_2)

g <- Augmented_beta(L,k,phi)
edges <- as_edgelist(g)
n_edge <- ecount(g)
to_sample <- floor(n_edge*pc_k_1)
present <- sample(1:n_edge,replace = FALSE,size = to_sample)
ids <- get_edge_ids(g, as.vector(t(edges[present,])))
sub_susc <- subgraph_from_edges(g,ids)
largest_comp <- largest_component(sub_susc)
v_in <- vcount(largest_comp) 
prop_in_largest <- v_in/L

proportion_in_largest <- prop_in_largest

#prop_pc_k_1 <- prop_in_largest
#prop_pc_k_2 <- prop_in_largest
#c(prop_pc_k_1,prop_pc_k_2)

# can add a separate point for each critical property

plot(p,prop_large_k_1,ylab = 'Proportion of nodes in GCC',
     xlab = 'Probability of edges present') 
lines(p,prop_large_k_1)
points(p,prop_large_k_2, pch = 2)
lines(p,prop_large_k_2)
points(critical_probs,c(prop_pc_k_1,prop_pc_k_2), pch = c(1,2), col = 'red')

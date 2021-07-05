
# test freefalling --------------------------------------------------------
set.seed(1)
comm <- matrix(sample(c(0, 1), 25, replace = TRUE), 
               nrow = 5, 
               ncol = 5, 
               dimnames = list(paste("comm", 1:5, sep = ""), paste("sp", 1:5, sep = "")
                               )
               )
phylo <- geiger::sim.bdtree(n = 5, stop = "taxa")

test_autoRegress <- free.falling(comm = comm, 
             phylo = phylo, 
             binary = TRUE, 
             test = TRUE, 
             nperm = 1000, 
             parallel = 4)

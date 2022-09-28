using Parameters, Plots #import the libraries we want
include("PS3_MD.jl") #import the functions that solve our model


prim, res = Initialize() #initialize primitive and results structs
mass = Initialize_3()    #initialize the population distribution
res2 = Initialize_2()    #initialize r, w, b
@elapsed Bellman(prim, res, res2) #solve the Bellman!
Distribution(prim, res, res2, mass)
@unpack val_func, pol_func, labor = res
@unpack a_grid = prim
@unpack mu, mu_dist, mu1 = mass


#graphing the value function
Plots.plot(a_grid, val_func[:, :, 50], title="Value Function", label = ["high productivity" "low productivity"])
#The value function is increasing and concave.

#graphing the savings function
Plots.plot(a_grid, pol_func[:, :, 20], title="Savings Function", label = ["high productivity" "low productivity"])


Plots.plot(a_grid, mu_dist[:, :, 65], title = "Population Distribution at J = 65")

println("All done!")
################################

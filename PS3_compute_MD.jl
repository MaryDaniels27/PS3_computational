using Parameters, Plots #import the libraries we want
include("PS3_MD.jl") #import the functions that solve our model


prim, res = Initialize() #initialize primitive and results structs
mass = Initialize_3()    #initialize the population distribution
res2 = Initialize_2()    #initialize r, w, b
@elapsed Bellman(prim, res, res2) #solve the Bellman!
@elapsed Distribution()  #solve for the stationary distribution
MarketClearing(prim, res, mass)
@unpack val_func, pol_func, labor = res
@unpack a_grid = prim
@unpack mu, mu_dist, mu1 = mass


#graphing the value function
Plots.plot(a_grid, val_func[:, :, 50], title="Value Function", label = ["high productivity" "low productivity"])
#The value function is increasing and concave.

#graphing the savings function
Plots.plot(a_grid, pol_func[:, :, 20], title="Savings Function", label = ["high productivity" "low productivity"])

#just checking the population
Plots.plot(a_grid, mu1[:, :, 65], title = "Population Distribution at J = 65")



MarketClearing(prim, res, mass)
K_1 = MarketClearing(prim,res,mass)[1]
L_1 = MarketClearing(prim,res,mass)[2]

#aggregate capital is 3.9967 which is a little too high
#aggregate labor is 0.325 which is too low


function Solve(prim::Primitives, res::Results, res2::Results_2, mass::mu_results)
    Prices(prim, res2)   #Define prices by solving the firm's problem and GBC
    println("Capital is", res2.K_0)
    println("Labor is", res2.L_0)
    Bellman(prim, res, res2) #solve the Bellman!
    Distribution()   #solve for the stationary distribution
    MarketClearing(prim, res, mass)
    K_1 = MarketClearing(prim,res,mass)[1]
    L_1 = MarketClearing(prim,res,mass)[2]
    return K_1, L_1
end

Solve(prim, res, res2, mass)



##Now we will update the prices to try to get a better aggregate capital and aggregate labor 
function Iterate(prim::Primitives, res::Results, res2::Results_2, mass::mu_results; tol::Float64 = 0.001)
    @unpack K_0, L_0 = res2  #Unpack initial guess of K and L
    Solve(prim, res, res2, mass)   #solve the model with the initial guess


    #Solving aggregate capital and labor
    K_1 = MarketClearing(prim,res,mass)[1]  #define the new aggregate capital which is found from market clearing
    L_1 = MarketClearing(prim,res,mass)[2]   #new aggregate labor from market clearing

    println("Aggregate capital initially", K_1)
    println("Aggregate labor initially", L_1)

    #calculate the errors for capital and labor separately
    err_k = abs(K_1 - res2.K_0)
    err_l = abs(L_1 - res2.L_0)
    n = 0  #counter
    
    while err_k > tol || err_l > tol
        println("on iteration", n)
        if abs(err_k) > tol
            res2.K_0 = 0.9res2.K_0 + 0.1K_1
        end

        if abs(err_l) > tol
            res2.L_0 = 0.9res2.L_0 + 0.1L_1
        end

        #After updating the initial guess of K and L, we now want to re-run the model
        Solve(prim, res, res2, mass)
        K_1 = MarketClearing(prim,res,mass)[1] 
        L_1 = MarketClearing(prim,res,mass)[2]

        #update the errors
        err_k = abs(K_1 - res2.K_0)
        err_l = abs(L_1 - res2.L_0)

        n+=1
        println("Aggregate capital now", K_1)
        println("Aggregate labor now", L_1)
        
    end  #end while loop
end #end iteration function 

Iterate(prim, res, res2, mass)

























println("All done!")
################################

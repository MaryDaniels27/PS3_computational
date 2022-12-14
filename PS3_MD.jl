### PS 3 ###
using Parameters, Plots #import the libraries we want

#creating the capital grid which has more density at lower values
i_max = 75
i_min = 0
nk = 1000
inc = (i_max - i_min)/(nk - 1).^2  #increment
aux = collect(1:1:nk)
assets = zeros(length(aux))
for i = 1:length(aux)
    assets[i] = i_min + inc*(aux[i] - 1).^2
end
assets

@with_kw struct Primitives
    N::Int64 = 66 #live span of agents which is also the number of generations in our model
    n::Float64 = 0.011 #population growth rate 
    R::Int64 = 46 #Retirement age
    tw::Int64 = 45 #length of your working life
    tR::Int64 = N - (R+1) #length of your retirement
    θ::Float64 = 0.11 #labor income tax
    γ::Float64 = 0.42 #weight on consumption
    σ::Float64 = 2 #coefficient of relative risk aversion
    η::Array{Float64, 1} = [0.59923239, 
    0.63885106, 
    0.67846973, 
    0.71808840, 
    0.75699959,
    0.79591079, 
    0.83482198, 
    0.87373318,
    0.91264437, 
    0.95155556, 
    0.99046676, 
    0.99872065, 
    1.0069745, 
    1.0152284, 
    1.0234823, 
    1.0317362, 
    1.0399901, 
    1.0482440, 
    1.0564979, 
    1.0647518, 
    1.0730057, 
    1.0787834, 
    1.0845611, 
    1.0903388, 
    1.0961165, 
    1.1018943, 
    1.1076720, 
    1.1134497, 
    1.1192274, 
    1.1250052, 
    1.1307829, 
    1.1233544, 
    1.1159259, 
    1.1084974, 
    1.1010689, 
    1.0936404, 
    1.0862119, 
    1.0787834, 
    1.0713549, 
    1.0639264,
    1.0519200,
    1.0430000,
    1.0363000,
    1.0200000,
    1.0110000]
    z::Array{Float64, 1} = [3.0, 0.5] #productivity levels which can either be high z = 3 or low z = 0.5
    markov::Array{Float64, 2} = [0.9261 0.0739; 0.0189 0.9811]  #markov matrix for productivity
    α::Float64 = 0.36  #capital share of output
    δ::Float64 = 0.06 #depreciation rate
    β::Float64 = 0.97 # discount factor
    a_grid::Array{Float64, 1} = assets #asset grid
    na::Int64 = length(a_grid) #length of the asset grid
    nz::Int64 = length(z) #length of productivity vector
end


#Note to self: If you want to see something inside of the primitives struct, type "println(object)" into the struct itself then run the code for the struct
#Note to self: Array{Float64, 3} tells Julia this is a three dimensional array

mutable struct Results
    val_func::Array{Float64, 3} #value function, we need three dimensions because we need a value function for each age and productivity state
    pol_func::Array{Float64, 3} #policy function for assets
    labor::Array{Float64, 3}  #optimal labor supply
    
    #prices, capital, labor
    r::Float64 #interest rate
    w::Float64 #wage rate
    b::Float64  #social security benefit
    K_0::Float64 #aggregate capital intial guess
    L_0::Float64 #aggregate labor initial guess
end

function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.na, prim.nz, prim.N) #initial value function guess
    pol_func = zeros(prim.na, prim.nz, prim.N) #initial policy function for asset guess
    labor = zeros(prim.na, prim.nz, prim.tw) #inital guess of the optimal labor supply

    K_0 = 3.64  #initial guess of capital
    L_0 = 0.43  #initial guess of labor
    w = 1.05
    r = 0.05
    b = 0.2
    res = Results(val_func, pol_func, labor, r, w, b, K_0, L_0) #initialize results struct
    prim, res #return deliverables
end

#Initializing the population distribution 

mutable struct mu_results
    mu::Array{Float64, 1}
    mu_dist::Array{Float64, 3}
    mu1::Array{Float64, 3}
end

function Initialize_3()
    @unpack N, n, na, nz = prim
    mu_dist = zeros(prim.na, prim.nz, prim.N) #this will hold our unweighted distribution
    mu1 = zeros(prim.na, prim.nz, prim.N)  #this will hol our final, properly weight distribution
    mu = ones(prim.N)
    for i = 2:N
        mu[i] = mu[i-1]/(1 + n) #finding the relative sizes of each cohort (accounting for population growth)
    end
    mu = mu/sum(mu)  #normalizing mu so that it sums to 1
    mass = mu_results(mu, mu_dist, mu1)
end


function Prices(prim::Primitives, res::Results, mass::mu_results)
    @unpack α, δ, θ, R, N = prim
    @unpack mu = mass
    @unpack K_0, L_0 = res
    w = (1-α)*(L_0)^(-α)*((K_0)^α)  #wage rate
    r_k = α*(K_0)^(α - 1)*(L_0)^(1-α)  #the firm's interest rate
    #Using the government budget constraint, we can solve for b
    b = (θ*w*(L_0))/(sum(mu[R:N]))
    r = r_k - δ  #the household's interest rate which is the firm's interest rate adjusted for depreciation
    res.w = w
    res.r = r
    res.b = b
end


#Solving the dynamic programming problem of retirees and workers

function Bellman(prim::Primitives, res::Results)
    @unpack N, R, tw, θ, γ, σ, η, z, markov, β, a_grid, na, nz = prim
    @unpack val_func, r, w, b = res
    
    #Solving the problem of retirees
    for a_index = 1:na, z_index = 1:nz #looping over assets and states
        c_N = (1 + r)*a_grid[a_index] + b    #last period consumption
        val_N = (c_N^((1-σ)*γ))/(1-σ)          #last period utility
        res.val_func[a_index, z_index, 66] = val_N #storing last period utility in the value function
        res.pol_func[a_index, z_index, 66] = 0.0 #storing the last period policy function
    end
    for j = (N-1):-1:R
        for a_index = 1:na, z_index = 1:nz  #looping over assets today
            a = a_grid[a_index] #setting the value of a
            candidate_max = -Inf  #initial guess of candidate max
            budget = (1 + r)*a + b #calculate the budget
            for ap_index = 1:na #looping over assets tomorrow
                c = budget - a_grid[ap_index] #consumption given a' selection
                if c>0  #check for positivity
                    val = (c^((1-σ)*γ))/(1-σ) + β * res.val_func[ap_index, z_index, j+1] #calculate the value function while looking at next period's value function
                    if val > candidate_max
                        candidate_max = val  #if our value function association with a' is greater than the current candidate max, update candidate max 
                        #to be equal to val. Then repeat the process for the next value of a' and update candidate max if this value of a' yields a higher
                        #value function. If the next value of a' does not yield a higher val, then candidate max will not update and we will just be left with
                        #the candidate max from before.
                        res.pol_func[a_index, z_index, j] = a_grid[ap_index] #update the policy function
                        res.val_func[a_index, z_index, j] = val  #updating the value function
                    end
                end
            end
        end
    end

    #Solving the problem of the working age people
    for i = tw:-1:1
        for a_index = 1:na  #looping over assets today
            a = a_grid[a_index] #setting the value of a
            for z_index = 1:nz #looping for productivity states
                candidate_max = -Inf  #initial guess of candidate max
                for ap_index = 1:na #looping over assets tomorrow
                    l = (γ*(1-θ)*(z[z_index] * η[i])*w - (1-γ)*((1+r)*a_grid[a_index] - a_grid[ap_index]))/((1-θ)*w*(z[z_index] * η[i]))  #for a given combination (z, a, a'), this is the 
                    #optimal labor supply
                    if l > 1  #this if, else loop ensures that labor supply is bounded between 0 and 1
                        l = 1
                    elseif l < 0
                        l = 0
                    end
                    budget = w*(1-θ)*(z[z_index] * η[i])*l + (1 + r)*a  #calculate the budget
                    c = budget - a_grid[ap_index] #consumption given a' selection
                        if c>0  #check for positivity
                            val = (((c^γ)*((1-l)^(1-γ)))^(1-σ))/(1-σ) + β * sum(res.val_func[ap_index, :, i+1].* markov[z_index, :]) #calculate the value 
                            #function while looking at next period's value function
                                if val > candidate_max
                                    candidate_max = val  #update candidate max
                                    res.pol_func[a_index, z_index, i] = a_grid[ap_index] #update the policy function
                                    res.val_func[a_index, z_index, i] = val  #updating the value function
                                    res.labor[a_index, z_index, i] = l  #updating the optimal labor supply vector 
                                end
                        end
                end
            end
        end
    end
end

function Distribution(prim::Primitives, res::Results, mass::mu_results)
    @unpack N, R, tw, z, a_grid, na, nz, markov = prim
    @unpack pol_func = res
    @unpack mu, mu_dist, mu1 = mass

    mu_dist[1, 1, 1] = 0.2037  #initial mass of high productivity people
    mu_dist[1, 2, 1] = 0.7963   #initial mass of low productivity people
    for j = 2:prim.N  #loopig over the ages
        for ap = 1:prim.na  #looping over assets tomorrow
            for z = 1:2  #looping over z from today
                d = findall(x->x == a_grid[ap], res.pol_func[:, z, j-1])
                    for i = 1:length(d) #loop over the indices you find for a specific z
                        for zp = 1:2  #loop over zp 
                            mu_dist[ap, zp, j] = mu_dist[ap, zp, j] + mu_dist[d[i], z, j-1] .* markov[z, zp]
                        end
                    end
            end
        end
    end

    #The final step is weighting each generation according to the relative sizes calculated by mu.
    for i = 1:prim.N
        mu1[:, :, i] = mu_dist[:, :, i] .* mu[i]
    end
    mass.mu_dist = mu_dist  #updating the unweighted distribution
    mass.mu1 = mu1  #updating the weighted distribution 
end

###  Market Clearing  ###

#Calculating aggregate capital
#1. for each generation, sum the two states together in the distribution
#2. multiply the agggregate distribution for each asset by the asset holding
function MarketClearing(prim::Primitives, res::Results, mass::mu_results)
    @unpack labor = res
    @unpack a_grid, N, na, nz, tw, z, η = prim
    @unpack mu1 = mass

    agg_dist = zeros(prim.na, prim.N) #In this array each column corresponds to a different generation. Each row is contains the result of summing the distribution at a given
    #asset holding level then multiplying that by the asset holding amount. 
    #3. sum all the generations together to get aggregate capital
    for i = 1:prim.N
        for j = 1:prim.na
        agg_dist[j, i] = sum(mu1[j, :, i]) .* a_grid[j]
        end
    end
    agg_capital = sum(agg_dist)

    #Calculating aggregate labor 
    #begin by calculating e(z, η_j) as vector
    e = zeros(prim.tw, prim.nz)
    for i = 1:prim.tw
        for j = 1:prim.nz
            e[i, j] = prim.z[j] * prim.η[i]
        end
    end

    agg_l = zeros(prim.na, prim.tw)  #each column corresponds to a different generation and each row corresponds to a different asset holding position. The first element agg_l[1,1] is
    #for generation 1, the sum over productivity states of mu1*e(z, η_j)*labor where this multiplication occurs for each productivity state first then we sum.
    for i = 1:prim.tw  #pick a generation
        for a = 1:prim.na  #pick an asset holding
            for j = 1:prim.nz  #pick a productivity state
                agg_l[a, i]= sum(mu1[a,:,i].*e[i,:].*labor[a,:,i])
            end
        end
    end
    agg_labor = sum(agg_l)

    return agg = [agg_capital, agg_labor]
end

function Solve_model(prim::Primitives, res::Results, mass::mu_results)
    Prices(prim, res, mass)   #Define prices by solving the firm's problem and GBC
    println("Capital is", res.K_0)
    println("Labor is", res.L_0)
    Bellman(prim, res) #solve the Bellman!
    Distribution(prim, res, mass)   #solve for the stationary distribution
    MarketClearing(prim, res, mass)
    K_1 = MarketClearing(prim,res,mass)[1]
    L_1 = MarketClearing(prim,res,mass)[2]
    return K_1, L_1
end

prim, res = Initialize()
mass = Initialize_3()

##Now we will update the prices to try to get a better aggregate capital and aggregate labor 
function Iterate(prim::Primitives, res::Results, mass::mu_results; tol::Float64 = 0.001)
    @unpack K_0, L_0 = res  #Unpack initial guess of K and L
    prim, res = Initialize()
    mass = Initialize_3()
    Solve_model(prim, res, mass)   #solve the model with the initial guess
    #Solving aggregate capital and labor
    K_1 = MarketClearing(prim,res,mass)[1]  #define the new aggregate capital which is found from market clearing
    L_1 = MarketClearing(prim,res,mass)[2]   #new aggregate labor from market clearing

    println("Aggregate capital initially", K_1)
    println("Aggregate labor initially", L_1)

    #calculate the errors for capital and labor separately
    err_k = abs(K_1 - res.K_0)
    err_l = abs(L_1 - res.L_0)
    n = 0  #counter
    
    while err_k > tol || err_l > tol
        println("on iteration", n)
        if abs(err_k) > tol
            res.K_0 = 0.9res.K_0 + 0.1K_1
        end

        if abs(err_l) > tol
            res.L_0 = 0.9res.L_0 + 0.1L_1
        end

        #After updating the initial guess of K and L, we now want to re-run the model
        mass = Initialize_3()   #re-initialize the population distribution
        Solve_model(prim, res, mass)
        K_1 = MarketClearing(prim,res,mass)[1] 
        L_1 = MarketClearing(prim,res,mass)[2]

        #update the errors
        err_k = abs(K_1 - res.K_0)
        err_l = abs(L_1 - res.L_0)

        n+=1
        println("Aggregate capital now", K_1)
        println("Aggregate labor now", L_1)
        
    end  #end while loop
end #end iteration function 

Iterate(prim, res, mass)
##Converged in 51 iterations. The aggregate capital is 3.364317 and the aggregate labor is 0.3432





### PS 3 ###
using Parameters, Plots, Optim

#creating the capital grid which has more density at lower values
max = 14
min = 0
nk = 180
inc = (max - min)/(nk - 1).^2  #increment
aux = collect(1:1:nk)
assets = zeros(length(aux))
for i = 1:length(aux)
    assets[i] = min + inc*(aux[i] - 1).^2
end


@with_kw struct Primitives
    N::Float64 = 66 #live span of agents which is also the number of generations in our model
    n::Float64 = 0.011 #population growth rate 
    R::Float64 = 46 #Retirement age
    tw::Int64 = 45 #length of your working life
    tR::Int64 = N - (R+1) #length of your retirement
    r_age::Array{Float64, 1} = collect(R:1:N)
    w_age::Array{Float64, 2} = collect(1:1:tw)
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
end

mutable struct Results_2
    r::Float64 #interest rate
    w::Float64 #wage rate
    b::Float64  #social security benefit
end


function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.na, prim.nz, prim.N) #initial value function guess
    pol_func = zeros(prim.na, prim.nz, prim.N) #initial policy function for asset guess
    labor = zeros(prim,na, prim.nz, prim.N) #inital guess of the optimal labor supply
    res = Results(val_func, pol_func, labor) #initialize results struct
    prim, res #return deliverables
end

function Initialize_2()
    r = 0.05 #initial guess of the interest rate
    w = 1.05 #initial guess of the wage rate
    b = 0.2 #initial guess of social security benefit
    res2 = Results_2(r, w, b) #initialize our economy's Parameters
end

#Initializing the population distribution 

mutable struct mu_results
    mu::Array{Float64, 1}
end

function Initialize_3()
    @unpack N, n = prim
    mu = ones(prim.N,1)
    for i = 2:N
        mu[i] = mu[i-1]/(1 + n) #finding the relative sizes of each cohort (accounting for population growth)
    end
    mu = mu/sum(mu)  #normalizing mu so that it sums to 1
    mass = mu_results(mu)
end


#Solving the dynamic programming problem of retirees and workers

function Backwards_Induct(prim::Primitives,res::Results)
    #First we fill in the values for the retired people
    for j in reverse(r_age)  #looping over the retirement ages in reverse order
        println("working on retiree", j)  #Tells us which age we are working on
        res.val_func[:, :, j] = Bellman(prim, res, res2, j, p) #filling in the value function using the results of the Bellman below
    end
    #Now we fill in the values for workers
    for p in reverse(w_age) #looping over the working ages in reverse order
        println("working on worker", p) 
        res.val_func[:, :, p] = Bellman(prim, res, res2, j, p)
    end
end

function Bellman(prim::Primitives, res::Results, res2::Results_2, j::Int64, p::Int64)
    @unpack N, n, R, tw, tR, r_age, w_age, θ, γ, σ, η, z, markov, α, δ, β, a_grid, na, nz = prim
    @unpack val_func = res
    @unpack r, w, b = res2

    v_next = zeros(na, nz, N) #filling in the next guess of the value function
    
    #Solving the problem of retirees
    if j == N
        for a_index = 1:na, z_index = 1:nz #looping over assets and states
            c_N = (1 + r)a_grid[a_index] + b    #last period consumption
            val_N = (c_N^((1-σ)γ))/(1-σ)          #last period utility
            res.val_func[a_index, z_index, j] = val_N #storing last period utility in the value function
            res.pol_func[a_index, z_index, j] = 0.0
        end
    end 
    elseif j!=N
        for a_index = 1:na, z_index = 1:nz  #looping over assets today
            obj(ap_index) = -obj_func(prim, res2, ap_index)
            opt = optimize(obj, 0.0, 14.0)  #The boundaries on a' are 0 and 14 from the capital grid
            v_next[a_index, z_index, j] = -opt.minimum  #fill in the maximum value 
            res.pol_func[a_index, z_index, j] = opt.minimizer  #fill in the choice of a' that gives the maximum value
        end
    end

    #Solving the problem of the working age people
    for a_index = 1:na, z_index = 1:nz  #looping over assets today
        obj_2(ap_index) = -obj_func_2(prim, res2, ap_index)
        opt_2 = optimize(obj_2, 0.0, 14.0)  #The boundaries on a' are 0 and 14 from the capital grid
        v_next[a_index, z_index, p] = -opt.minimum  #fill in the maximum value 
        res.pol_func[a_index, z_index, p] = opt.minimizer  #fill in the choice of a' that gives the maximum value
    end

end


#Objective function

function obj_func(prim::Primitives, res2::Results_2, ap_index::Float64)
    @unpack N, n, R, tw, tR, r_age, θ, γ, σ, η, z, markov, α, δ, β, a_grid, na, nz = prim
    @unpack r, w, b = res2
    for a_index = 1:na  #looping over assets today
        a = a_grid[a_index] #setting the value of a
        budget = (1 + r)a + b 
        for ap_index = 1:na #looping over assets tomorrow
            c = budget - a_grid[ap_index] #consumption given a' selection
            if c>0  #check for positivity
                val = (c^((1-σ)γ))/(1-σ) + β * v_next[ap_index, :, i+1] #calculate the value function
            end
        end
    end
end

function obj_func(prim::Primitives, res2::Results_2, ap_index::Float64)
    @unpack N, n, R, tw, tR, r_age, θ, γ, σ, η, z, markov, α, δ, β, a_grid, na, nz = prim
    @unpack r, w, b = res2
    for a_index = 1:na  #looping over assets today
        a = a_grid[a_index] #setting the value of a
        budget = (1 + r)a + b 
        for ap_index = 1:na #looping over assets tomorrow
            c = budget - a_grid[ap_index] #consumption given a' selection
            if c>0  #check for positivity
                val = (c^((1-σ)γ))/(1-σ)  #calculate the value function
            end
        end
    end
end








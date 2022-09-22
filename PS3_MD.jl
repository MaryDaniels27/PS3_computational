### PS 3 ###
using Parameters, Plots

@with_kw struct Primitives
    N::Float64 = 66 #live span of agents which is also the number of generations in our model
    n::Float64 = 0.011 #population growth rate 
    R::Float64 = 46 #Retirement age
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
    a_grid::Array{Float64, 1} = collect(range(0, length = 180, stop = 14)) #asset grid
    na::Int64 = length(a_grid) #length of the asset grid
    nz::Int64 = length(z) #length of productivity vector
end

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



using Parameters, Plots #import the libraries we want
include("PS3_MD.jl") #import the functions that solve our model


prim, res = Initialize() #initialize primitive and results structs
mass = Initialize_3()    #initialize the population distribution
res2 = Initialize_2()    #initialize r, w, b
@elapsed Bellman(prim, res, res2) #solve the Bellman!
Distribution(prim, res, res2, mass)
@unpack val_func, pol_func, labor = res
@unpack a_grid, markov, N = prim
@unpack mu, mu_dist = mass


#graphing the value function
Plots.plot(a_grid, val_func[:, :, 50], title="Value Function", label = ["high productivity" "low productivity"])
#The value function is increasing and concave.

#graphing the savings function
Plots.plot(a_grid, pol_func[:, :, 20], title="Savings Function", label = ["high productivity" "low productivity"])


###### Below, I go through a step by step process of creating the stationary population distribution for this model. I begin with 
# a slower code that is more broken down and easy to understand. Then I speed up the code by adding for loops. I add each for loop
# one at a time, so we can see exactly what each loop is doing and how it is speeding up the process. ######################## 

res.pol_func[1, :, 1]     
findall(x-> x == res.pol_func[1, 1, 1], a_grid)
#the policy function says that an agent with a = 0 today and z = high should hold a' = 2.2026 tomorrow. We can find this value of a' indexed
#in the asset grid at element number 72.
findall(x-> x == res.pol_func[1, 2, 1], a_grid)
#likewise, if you're in the low state today a' = 0.367 which can be found as the 30th element of the asset grid
p = findall(x->x == a_grid[72], res.pol_func[:, 1, 1])
#for the first twelve asset holdings in our asset grid, if z = high and j = 1, then each agent should choose to hold a' = 2.20261 tomorrow
mu_dist[p[1], 1, 1]
#we have a mass of 0.2037 at the first asset holding position in our grid
mu_dist[p[2], 1, 1]
g = 0
for i = 1:length(p)
            g = g + mu_dist[p[i], 1, 1]*0.9261
end
g

h = findall(x->x == a_grid[72], res.pol_func[:, 2, 1]) #z = low people who will choose a' = 2.20261 tomorrow
#these are people at the 65 index in the asset grid. However, we know there's no one in the initial distribution at the 65th index, so 
#this shouldn't add any mass to our current population
v = zeros(1, 2)
for i = 1:length(p)
    #z = high today, a' = 2.20261 tomorrow
    for zp = 1:2  #state tomorrow
        v[1, zp] = v[1, zp] + mu_dist[p[i], 1, 1] .* markov[1, zp]  #populate the row vector for a' = 2.20261 tomorrow with the people coming for z = high today
    end
end
for i = 1:length(h)
    #z = low today, a' = 2.20261 tomorrow
    for zp = 1:2 #state tomorrow
        v[1, zp] = v[1, zp] + mu_dist[h[i], 2, 1] .* markov[2, zp]  #add to the row vector for a' = 2.20261 the people who are coming from z = low today
    end
end

v
#as predicted, adding the people coming from the low state didn't change the mass at all. So the above code is doing what I want. Now, I want to
#try to speed up the process a bit 

output = zeros(1, 2)
for z = 1:2  #looping over z from today
    d = findall(x->x == a_grid[72], res.pol_func[:, z, 1])
        for i = 1:length(d) #loop over the indices you find for a specific z
            for zp = 1:2  #loop over zp 
                output[1, zp] = output[1, zp] + mu_dist[d[i], z, 1] .* markov[z, zp]
            end
        end
end
output
#these are the masses that would be associated with a' = 2.20261 tomorrow which would occur at element 72 in the asset grid. 

#yay! It populated the same vector but with less lines of code!

#So far, we've found that holding a' constant, we can populate a row vector associated with it. Great. Now, we want to loop over a' and populate
#all the row vectors

output_2 = zeros(prim.na, prim.nz)    #this array will hold our results for all asset holdings
for ap = 1:prim.na  #looping over assets tomorrow
    for z = 1:2  #looping over z from today
        d = findall(x->x == a_grid[ap], res.pol_func[:, z, 1])
            for i = 1:length(d) #loop over the indices you find for a specific z
                for zp = 1:2  #loop over zp 
                    output_2[ap, zp] = output_2[ap, zp] + mu_dist[d[i], z, 1] .* markov[z, zp]
                end
            end
    end
end

output_2
#yay! It ran and when I call row 72 which is associated with a' = 2.20261, we get the correct population mass. Recall for early in the 
#script that the policy function recommends anyone in j = 1 who is in the low state to choose a' = 0.367 which is the 30th index in the
#asset grid. Well I checked row 30 of our final distribution and it has the correct population mass! The code is doing what we want it to!

#now we update our mu distribution
mu_dist[:, :, 2] = output_2

#At this point, we've calculated the distribution for period j = 2 using the distribution for j = 1. Now, if we add a for loop that looks back one generation,
#we should be able to repeat the process for all the generations. Before I do this though, I will populate the distribution for j = 3.

#first, checking the policy function, we'll see where it recommends people from j=2 with a = 2.20261 go and where people with a = 0.367 go since this is where our current
#population is located. 

findall(x-> x == res.pol_func[72, 1, 2], a_grid)
#the policy function says that an agent with a = 2.20261 today and z = high should hold a' = 4.546 tomorrow. We can find this value of a' indexed
#in the asset grid at element number 103.
findall(x-> x == res.pol_func[72, 2, 2], a_grid)
#the policy function says that an agent with a = 2.20261 today and z = low should hold a' = 2.658 tomorrow. We can find this value of a' indexed
#in the asset grid at element number 79.

#Thus, we can expect population masses at indices 103 and 79 of our j = 3 grid. 

output_3 = zeros(prim.na, prim.nz)    #this array will hold our results for all asset holdings for j = 3
for ap = 1:prim.na  #looping over assets tomorrow
    for z = 1:2  #looping over z from today
        d = findall(x->x == a_grid[ap], res.pol_func[:, z, 2])
            for i = 1:length(d) #loop over the indices you find for a specific z
                for zp = 1:2  #loop over zp 
                    output_3[ap, zp] = output_3[ap, zp] + mu_dist[d[i], z, 2] .* markov[z, zp]
                end
            end
    end
end
output_3
output_3[103, :]     #there's a mass! A little math on the calculator confirms this is the correct mass too.

#Okay, now I feel pretty confident about my code. What I want to do now is wrap this up by adding the for loop for ages

#out mu_distribution is initialized at zero already so we don't need to do that again
#filling in the initial mass
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

mu_dist

#Checking that we have masses at indices 72 and 30 for j = 2
mu_dist[72, :, 2]   #it has the correct mass!
mu_dist[30, :, 2]   #it has the correct mass!

#checking that we have (correct) masses at indices 103 and 79 for j = 3
mu_dist[103, :, 3]   #it has the right mass!
mu_dist[79, :, 3]    #you already know what it is, ah!

#The final step is weighting each generation according to the relative sizes calculated by mu. I'm just going to add the array mu1 into my mutable struct
#and initialize it to hold my final, weighted population distribution.
mu1 = zeros(prim.na, prim.nz, prim.N)
for i = 1:prim.N
    mu1[:, :, i] = mu_dist[:, :, i] .* mu[i]
end

mu1
mu1[72, :, 2]    #yup, this is the correct weighted mass
sum(mu1)     #now the distribution sums to 1 which is what want. We're done! Take a bow



#Just checking the distribution to see where the mass is at in the final period
Plots.plot(a_grid, mu_dist[:, :, 66], title = "Population Distribution at J = 66")






println("All done!")
################################

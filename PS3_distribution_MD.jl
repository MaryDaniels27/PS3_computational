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
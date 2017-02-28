using Distributions

function evsi_an(μ = [0, 0], σ = [1, 1], nsamples = [10000, 10000])
    m   = μ[1] - μ[2]
    s   = σ .* σ
    ps  = nsamples .* s
    θ   = sqrt(s[1] * (ps[1] / (ps[1] + 1)) + s[2] * (ps[2] / (ps[2] + 1)))
    (θ * sqrt(2 / π) * exp((-m^2) / (2 * θ^2)) + m * erf(m / (θ * sqrt(2))) - abs(m)) / 2
end

function evsi_sim(μ = [0, 0], σ = [1, 1], nsamples = [10000, 10000]; nsims = 1000000)
  
    ndists        = length(μ)
    true_values   = Array(Float64, nsims, ndists)
    prior_weight  = Array(Float64, ndists)
    sample_weight = Array(Float64, ndists) 
    sample_sd     = Array(Float64, ndists)
   
    for dist in 1:ndists

        true_values[:, dist] = rand(Normal(μ[dist], σ[dist]), nsims)
        prior_weight[dist]   = (1 / σ[dist]) / (nsamples[dist] + (1 / σ[dist]))
        sample_weight[dist]  = nsamples[dist] / (nsamples[dist] + (1 / σ[dist]))
        sample_sd[dist]      = sqrt(1 / nsamples[dist])
    
    end
    
    benefit = Array(Float64, nsims)
    
    for sim in 1:nsims
        
        posterior_mean = Array(Float64, ndists)
        
        for dist in 1:ndists
            posterior_mean[dist] = μ[dist] * prior_weight[dist]
            if nsamples[dist] > 0 
                posterior_mean[dist] +=
                    rand(Normal(true_values[sim, dist], sample_sd[dist])) *
                        sample_weight[dist]
            end
        end
        
        benefit[sim] = true_values[sim, indmax(posterior_mean)]

    end    
    
    mean(benefit) - maximum(μ)

end
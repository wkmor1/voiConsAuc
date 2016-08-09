using Distributions, Optim, RCall, DataFrames

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
        prior_weight[dist]   = (1 / σ[dist]^2) / (nsamples[dist] + 1)
        sample_weight[dist]  = nsamples[dist] / (nsamples[dist] + 1)
        sample_sd[dist]      = sqrt(1 / nsamples[dist])
    
    end
    
    benefit = Array(Float64, nsims)
    
    for sim in 1:nsims
        
        posterior_mean = Array(Float64, ndists)
        
        for dist in 1:ndists
            posterior_mean[dist] = μ[dist] * prior_weight[dist] +
              rand(Normal(true_values[sim, dist], sample_sd[dist])) * sample_weight[dist]
        end
        
        benefit[sim] = true_values[sim, indmax(posterior_mean)]

    end    
    
    mean(benefit) - maximum(μ)

end

function evsi_opt(b, μ, σ, f)
    function fn(b, μ, σ, f)
        g(p) = -f(μ, σ, [b * p, b * (1 - p)])
        res  = optimize(g, 0, 1, method = GoldenSection())
        return hcat(Optim.minimizer(res), -Optim.minimum(res))
    end
    res = hcat(b, vcat([fn(x, μ, σ, f) for x = b]...), vcat([f(μ, σ, [.5x, .5x]) for x = b]...))
    res = convert(DataFrame, res)
    rename!(res, Dict(:x1 => :n, :x2 => :p, :x3 => :evsi_opt, :x4 => :evsi_naive))
end

R"
plot_evsi <- function(df, evpi) {
        with(df,
             {
               par(mfrow = c(3, 1), las = 1, mar = c(0.5, 4, 0, 0), oma = c(5, 0, 0, 0))
               plot(p ~ n, log = 'x', type = 'l', ylab = 'p', lwd = 2, ylim = 0:1, xaxt = 'n', xlab = '')
               abline(h = .5, lwd = 2, lty = 2, col = 'grey')
               plot(evsi_opt ~ n, log = 'x', type = 'l', lwd = 2, ylim = c(0, evpi), xaxt = 'n', xlab = '')
               points(n, evsi_naive, type = 'l', lwd = 2, lty = 2, col = 'grey')
               plot((evsi_opt - evsi_naive) ~ n, log = 'x', type = 'l', lwd = 2, ylim = c(0, evpi), xlab = '')
               mtext('n', 1, 3, outer = TRUE)
             }
        )
    }
"

function evsi_plot(μ = [0, 0], σ = [1, 1]; f = evsi_an,
                   b = exp(collect(linspace(log(1), log(5000), 20))))
    df   = evsi_opt(b, μ, σ, f)
    evpi = evsi_an(μ, σ)
    R"plot_evsi($df, $evpi)"
end
    
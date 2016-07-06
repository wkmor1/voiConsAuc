using Distributions

function posterior(x, m, s, p, nsp)
    s2   = s^2
    s2p1 = s2 * p + 1
    μ    = (m + x * s2 * p) / s2p1
    σ    = sqrt(s2 / s2p1)
    rand(Normal(μ, σ), nsp)
end

function preposterior(m, s, p, nps, npd)
    hcat(map(posterior, rand(Normal(m, s), npd))...)
end

function preposterior_mult(m, s, p, nps, npd)
    f(m, s, p) = preposterior(m, s, p, nps, npd)
    cat(3, map(f, m, s, p)...)
end

function evsi_sim(m = [0, 0], s = [1, 1], p = [0, 0], npd = 1000, nps = 1000)
    mean(preposterior_mult(m, s, p, npd, nps), 1)
end

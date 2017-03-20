pids = addprocs(10)

@everywhere begin 
    include("voiConsAuc.jl")
    f(n, p) = evsi_sim([1, 0], [1, .2], [n * p, n * (1 - p)], nsims = 10000000)
end

x = @parallel hcat for i = linspace(1, 50, 50)
  [f(i, j) for j = linspace(.99, .01, 50)]
end

rmprocs(pids)

writedlm("evsi2sim.csv", x, ',')

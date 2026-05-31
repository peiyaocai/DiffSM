using KLIEPInference
using ProximalBase, CoordinateDescent
using LinearAlgebra, SparseArrays, Statistics
using Distributions, StatsBase, JLD
using DelimitedFiles
using CSV, DataFrames

function ψlinear(X)
    m, n = size(X)
    p = div(m * (m + 1), 2)
    out = zeros(Float64, p, n)
    for i = 1:n
        ind = 0
        for row = 1:m
            for col = 1:row
                ind += 1
                if row == col
                    out[ind, i] = -0.5 * X[col, i] * X[row, i]
                else
                    out[ind, i] = -1.0 * X[col, i] * X[row, i]
                end
            end
        end
    end
    out
end

function trimap(i, j; diagonal=true)
    if i > j
        return trimap(j, i)
    else
        return i + div(j * (j - 1), 2)
    end
end

function itrimap(k; diagonal=true)
    j = convert(Int, ceil((-1. + sqrt(1. + 8. * k)) / 2.))
    i = k - div(j * (j - 1), 2)
    CartesianIndex(i, j)
end

function unpack(θ; diagonal=true)
    p = length(θ)
    m = convert(Int, ceil((-1. + sqrt(1. + 8. * p)) / 2.))
    out = zeros(m, m)
    for k = 1:p
        out[itrimap(k, diagonal=true)] = θ[k]
    end
    out + out' - diagm(diag(out))
end

n = parse(Int, ARGS[1])
p0 = parse(Int, ARGS[2])
sig = parse(Int, ARGS[3])
seed = parse(Int, ARGS[4])

xdir ="Data/XNull_n$(n)p$(p0)sig$(sig)seed$(seed).csv"
ydir ="Data/YNull_n$(n)p$(p0)sig$(sig)seed$(seed).csv"
x_data,nodesx = readdlm(xdir, ',', Float64, '\n'; header=true)
y_data,nodesy = readdlm(ydir, ',', Float64, '\n'; header=true)

rm(xdir)
rm(ydir)

Ψx = ψlinear(x_data)
Ψy = ψlinear(y_data)

p, ny = size(Ψy)
λ1 = 1.01 * quantile(Normal(), 1. - 0.05 / p)
λ2 = sqrt(2. * log(p) / ny)

θ = spKLIEP(Ψx, Ψy, λ1, CD_KLIEP(); loadings=true)
c1 = unpack(θ)
H = KLIEP_Hessian(spzeros(Float64, p), Ψy)
Hinv = Vector{SparseIterate{Float64}}(undef, p)
σ = Vector{Float64}(undef, p)
for k in 1:p
    print(k)
    ω = Hinv_row(H, k, λ2)

    supp = KLIEPInference._find_supp(k, ω)
    h = view(H, supp, supp)
    δ = (supp .=== k)
    ω[supp] = h\δ

    Hinv[k] = ω

    σ[k] = stderr_SparKLIE(Ψx, Ψy, θ, ω)
end

t_start = time()
boot = boot_SparKLIE(Ψx, Ψy, θ, Hinv; bootSamples=100,debias = 1)
println(time() - t_start)

crit = boot_quantile(boot, 0.95)

idx = findall(sqrt(ny) .* abs.(boot.θhat) .>= crit)
inf_graph = zeros(p)
inf_graph[idx] .= 1
inf_graph = unpack(inf_graph)

θd = boot.θhat
clb = θd - 1.96*σ
cub = θd + 1.96*σ
clb = unpack(clb)
cub = unpack(cub)

# define output file names
outfile_clb = "Res/clb_n$(n)p$(p0)sig$(sig)seed$(seed).csv"
outfile_cub = "Res/cub_n$(n)p$(p0)sig$(sig)seed$(seed).csv"
outfile_inf = "Res/inf_graph_n$(n)p$(p0)sig$(sig)seed$(seed).csv"
outfile_est = "Res/est_n$(n)p$(p0)sig$(sig)seed$(seed).csv"
# make sure directory exists
mkpath("Res")

# convert to DataFrame (safe for vectors/matrices)
CSV.write(outfile_clb, DataFrame(clb, :auto))
CSV.write(outfile_cub, DataFrame(cub, :auto))
CSV.write(outfile_inf, DataFrame(inf_graph, :auto))
CSV.write(outfile_est, DataFrame(c1, :auto))

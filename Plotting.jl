Plotting


# Pkg.add("BenchmarkTools")
# Pkg.add("Distributions")
# Pkg.add("CxxWrap")
# Pkg.add("JLD")
# Pkg.add("ImageIO")
# Pkg.add("CircularArrays")
# Pkg.add("OffsetArrays")
# Pkg.add("TensorOperations")
# Pkg.add("TensorCast")
# Pkg.add("Tullio")
# Pkg.add("DiffEqOperators")
# Pkg.add("BandedMatrices")

#using BenchmarkTools
#using Distributions
#using Random
#using Pkg
#using CxxWrap
#using JLD
#using ImageIO
#using Plots
#using PyPlot
#using Images
#using StaticArrays
#using CircularArrays
#using TensorCast
#using Tullio
#using DiffEqOperators
#using BandedMatrices

using BenchmarkTools
using Distributions
using Random
#using Plots
using PyPlot
#using Images
using CircularArrays
using CxxWrap
#using CosmicStrings

x = randn(1024,1024);
y = randn(1024,1024);
@time x = x.+y;
@time x.= x.+y;


n = 8;
Δx = 1.0;
Δt = 0.1;
η = 1.0;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;
t₀ = 0.1

for n ∈ 9:9

@time begin
N = 2^n;
println(N)

Ȧ₁ = zeros(N,N);
Ȧ₂ = zeros(N,N);

μ,σ = 0, 0.1;
C₁ = rand(Normal(μ,σ),N,N);
C₂ = rand(Normal(μ,σ),N,N);



plotting_2D!(N, t₀, 100, C₁, C₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)

end
end

for n ∈ 6:6

@time begin

N = 2^n;
println(N)
Ȧ₁ = zeros(N,N,N);
Ȧ₂ = zeros(N,N,N);

μ,σ = 0, 0.1;
C₁ = rand(Normal(μ,σ),N,N,N);
C₂ = rand(Normal(μ,σ),N,N,N);

plotting_3D!(N, t₀, 100, C₁, C₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)

end

end

t₀ = 0.1 ;
# t₂ = t/2 ;
# t₅ = t/5 ;
# t₁₀ = t/10 ;
# t₂₀ = t/20 ;
# t₅₀ = t/50 ;


n = 7;
Δx = 1.0;
Δt = 0.1;
η = 1.0;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;
t₀ = 0.1 ;
t₅ = t/5 ;

for n ∈ 1:8

@time begin

N = 2^n;
print(N)
Ȧ₁ = zeros(N,N,N);
Ȧ₂ = zeros(N,N,N);

μ,σ = 0, 0.1;
C₁ = rand(Normal(μ,σ),N,N,N);
C₂ = rand(Normal(μ,σ),N,N,N);

plotting_3D!(N, t₀, 100, C₁, C₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)

end

end
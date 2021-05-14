Plotting


Pkg.add("BenchmarkTools")
Pkg.add("Distributions")
Pkg.add("CxxWrap")
Pkg.add("JLD")
Pkg.add("ImageIO")
Pkg.add("CircularArrays")
Pkg.add("OffsetArrays")
Pkg.add("TensorOperations")
Pkg.add("TensorCast")
Pkg.add("Tullio")
Pkg.add("DiffEqOperators")
Pkg.add("BandedMatrices")

using BenchmarkTools
using Distributions
using Random
using Pkg
using CxxWrap
using JLD
using ImageIO
using Plots
#using PyPlot
using Images
#using StaticArrays
using CircularArrays
#using TensorCast
#using Tullio
#using DiffEqOperators
#using BandedMatrices

x = randn(1024,1024);
y = randn(1024,1024);
@time x = x.+y;
@time x.= x.+y;


n = 5;
Δx = 1.0;
Δt = 0.1;
η = 1.0;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;

Ȧ₁ = zeros(N,N);
Ȧ₂ = zeros(N,N);

μ,σ = 0, 0.1;
C₁ = rand(Normal(μ,σ),N,N);
C₂ = rand(Normal(μ,σ),N,N);
D₁ = C₁;
D₂ = C₂;
P = zeros(N,N);

t₀ = 0.1 ;
t₂ = t/2 ;
t₅ = t/5 ;
t₁₀ = t/10 ;
t₂₀ = t/20 ;
t₅₀ = t/50 ;


@time saving_circle(N, t₀, t₂₀, C₁, C₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)
Ȧ₁ = zeros(N,N);
Ȧ₂ = zeros(N,N);


@time saving_laplace!(N, t₀, t₂₀, D₁, D₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)

Ȧ₁ = zeros(N,N);
Ȧ₂ = zeros(N,N);

@time plotting_2D!(N, t₀, t₂₀, D₁, D₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)




n = 7;
Δx = 1.0;
Δt = 0.1;
η = 1.0;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;

Ȧ₁ = zeros(N,N,N);
Ȧ₂ = zeros(N,N,N);

μ,σ = 0, 0.1;
C₁ = rand(Normal(μ,σ),N,N,N);
C₂ = rand(Normal(μ,σ),N,N,N);
P = zeros(N,N,N);

t₀ = 0.1 ;
t₂ = t/2 ;
t₅ = t/5 ;
t₁₀ = t/10 ;
t₂₀ = t/20 ;
t₅₀ = t/50 ;

@time plotting_slow3D(N, t₀, t₅, C₁, C₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)

Ȧ₁ = zeros(N,N,N);
Ȧ₂ = zeros(N,N,N);

@time plotting_3D!(N, t₀, t₅, C₁, C₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)

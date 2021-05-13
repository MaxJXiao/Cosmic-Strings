Plotting


Pkg.add("BenchmarkTools")
Pkg.add("Distributions")
Pkg.add("CxxWrap")
Pkg.add("JLD")
Pkg.add("Images")
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
using PyPlot
using Images
using StaticArrays
using CircularArrays
#using TensorCast
#using Tullio
using DiffEqOperators
using BandedMatrices




n = 8;
Δx = 1.0;
Δt = 0.1;
η = 1.0;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;

Ȧ₁ = zeros(N,N);
Ȧ₂ = zeros(N,N);

μ,σ = 0, 0.1;
A₁ = rand(Normal(μ,σ),N,N);
A₂ = rand(Normal(μ,σ),N,N);

t₀ = 0.1 ;
t₂ = t/2 ;
t₅ = t/5 ;
t₁₀ = t/10 ;
t₂₀ = t/20 ;

@time begin
A₁ = rand(Normal(μ,σ),N,N);
A₂ = rand(Normal(μ,σ),N,N);
A₁,A₂ = saving_circle(N, t₀, t₂₀, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;
end

@time begin
    A₁ = rand(Normal(μ,σ),N,N,N);
    A₂ = rand(Normal(μ,σ),N,N,N);
    A₁,A₂ = saving_3D(N, t₀, 1, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;
end



@time begin

A₁ = rand(Normal(μ,σ),N,N);
A₂ = rand(Normal(μ,σ),N,N);
A₁,A₂ = saving(N, t₀, t₂₀, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;

end

@time begin
n = 8;
Δx = 1.0;
Δt = 0.1;
η = 1.0;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;
t₂₀ = t/20;

Ȧ₁ = zeros(N,N);
Ȧ₂ = zeros(N,N);

A₁ = rand(Normal(μ,σ),N,N);
A₂ = rand(Normal(μ,σ),N,N);
A₁,A₂ = saving_laplace!(N, t₀, t₂₀, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;
    
end

@time begin
n = 9;
Δx = 1.0;
Δt = 0.1;
η = 1.0;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;
t₂₀ = t/20;

Ȧ₁ = zeros(N,N);-
Ȧ₂ = zeros(N,N);

A₁ = rand(Normal(μ,σ),N,N);
A₂ = rand(Normal(μ,σ),N,N);
saving_circle(N, t₀, t₂₀, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;

end



@time begin
for i = 1:10
    A₁ = rand(Normal(μ,σ),N,N);
    A₂ = rand(Normal(μ,σ),N,N);
    A₁,A₂ = saving_roll(N, t₀, t₁₀, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;
end
end
    
@time begin
    for i = 1:10
        A₁ = rand(Normal(μ,σ),N,N);
        A₂ = rand(Normal(μ,σ),N,N);
        A₁,A₂ = saving_function(N, t₀, t₁₀, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;
    end
    end




@time begin
    Laplacian(10,zeros(10,10),1.0)
end
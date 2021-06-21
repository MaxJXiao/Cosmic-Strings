Plotting

# using Pkg
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
# Pkg.add("PyPlot")
# Pkg.add("Plots")
# Pkg.add("FFTW")
# Pkg.add("SciPy")
#Pkg.add("GR")

import GR
using BenchmarkTools
using Distributions
using Random
using Plots
using PyCall
using PyPlot
#using Images
using CircularArrays
using CxxWrap
using FFTW
using SciPy
stats = pyimport("scipy.stats")

x = randn(64,64,64);
y = randn(512,512);
@time x = x.+1;
@time y= y.+1;


n = 8;
Δx = 1;
Δt = 0.1;
η = 1;
ω = 5.0;
N = 2^n;
t = 0.5N*Δx/Δt;
t₀ = 0.1 ;
t₂ = t/2 ;
t₅ = t/5 ;
t₁₀ = t/10 ;
t₂₀ = t/20 ;
t₅₀ = t/50 ;


for n ∈ 9:9

    @time begin
        t₀ = 1e-10;
        N = 2^n;
        println(N)
        
        rₐ = range(15,stop = 24,length = 10);
        fₐ = 10 .^rₐ;
        
        for i ∈ 2:2
            print(string(i)*" ")
            println(fₐ[i])
        
            Ȧ₁ = zeros(N,N);
            Ȧ₂ = zeros(N,N);
        
            μ,σ = 0, 0.1;
            C₁ = rand(Normal(μ,σ),N,N);
            C₂ = rand(Normal(μ,σ),N,N);
        
        
        
            t₁ = run_2D!(N, t₀, 50, C₁, C₂, Ȧ₁, Ȧ₂,ω,η, Δx, Δt,i)
            plotting_2D!(N,t₀,t₁,250,C₁,C₂,Ȧ₁,Ȧ₂,ω,η, Δx,Δt,fₐ[i],i)
        
        end
        
    end
end


for n ∈ 9:9
    
    @time begin
    t₀ = 1;
    N = 2^n;
    println(N)
    
    rₐ = range(15,stop = 24,length = 10);
    fₐ = 10 .^rₐ;
    Δx = 1;
    Δt = 0.1;
    
    for i ∈ 2:2
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
    
        μ,σ = 0, 1e-10;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
    
    
        t₁ = PQrun_2D!(N, t₀, 100, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)
        #PQplotting_2D!(N,t₀,t₁,t₂₀,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
    end
    #PQ EPOCH
    end
end


for n ∈ 9:9

    @time begin
    b = 1
    N = 2^n;
    println(N)
    
    rₐ = range(15,stop = 24,length = 10);
    fₐ = 10 .^rₐ;
    s = 50
    
    
    for i ∈ 2:2
        t₀ = 1;
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
    
        μ,σ = 0, 0.1;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
        Δx = 1
        Δt = 0.1
    
        t₁ = EQCDrun_2D!(N, t₀, s, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)

        Δx = 1
        Δt = 1e-6

        EQCDplotting_2D!(N,s,t₁,5000Δt,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
    end
    
    end
end


for n ∈ 9:9

    @time begin
    b = 6.68
    N = 2^n;
    println(N)
    
    rₐ = range(15,stop = 24,length = 10);
    fₐ = 10 .^rₐ;
    s = 50;
    
    for i ∈ 2:2
        t₀ = 1
        print(string(i)*" ")
        println(fₐ[i])
        Δx = 1;
        Δt = 0.1;
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
        Ȧ = zeros(N,N);
        
    
        μ,σ = 0, 1e-10;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
    
        t₁ = Lrun_2D!(N, t₀, s, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)
        angle = zeros(N,N);
        angler!(angle,C₁,C₂);
        angler!(Ȧ,Ȧ₁,Ȧ₂)

        Δx = 1;
        Δt = 1e-8;
        #t₁ = 1.61e-10 * (fₐ[i]/1e12)
        Lplotting_2D!(N, s,t₁, 5000Δt, angle, Ȧ, Δx, Δt,fₐ[i],i)
        #PQplotting_2D!(N,t₀,t₁,t₂₀,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
    end
    
    end
end




for n ∈ 9:9

    @time begin
    b = 1
    N = 2^n;
    println(N)
    
    rₐ = range(15,stop = 24,length = 10);
    fₐ = 10 .^rₐ;
    s = 150
    
    
    for i ∈ 2:2
        t₀ = 1;
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
        Ȧ = zeros(N,N);
    
        μ,σ = 0, 0.1;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
        Δx = 1
        Δt = 0.1
    
        t₁ = FPQrun_2D!(N, t₀, s, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)

        Δx = 1
        Δt = 1e-6

        r = 10s

        t₁ = FErun_2D!(N,s,t₁,r*Δt,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)

        angle = zeros(N,N);
        angler!(angle,C₁,C₂);
        angler!(Ȧ,Ȧ₁,Ȧ₂)

        Δx = 1;
        Δt = 1e-8;

        FLrun_2D!(N, s + r/10,t₁, r*Δt, angle, Ȧ, Δx, Δt,fₐ[i],i)
    
    end
    
    end
end






for n ∈ 1:7

    @time begin

    N = 2^n;
    println(N)
    Ȧ₁ = zeros(N,N,N);
    Ȧ₂ = zeros(N,N,N);

    μ,σ = 0, 0.1;
    C₁ = rand(Normal(μ,σ),N,N,N);
    C₂ = rand(Normal(μ,σ),N,N,N);

    plotting_3D!(N, t₀, 400, C₁, C₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt)

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

for n ∈ 1:6

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
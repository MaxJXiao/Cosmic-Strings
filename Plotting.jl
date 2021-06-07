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

using BenchmarkTools
using Distributions
using Random
#using Plots
using PyCall
using PyPlot
#using Images
using CircularArrays
using CxxWrap
#using CosmicStrings

x = randn(64,64,64);
y = randn(512,512);
@time x = x.+1;
@time y= y.+1;


n = 8;
Δx = 1e-9;
Δt = 1e-10;
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
        
        rₐ = range(14.5,stop = 23.5,length = 10);
        fₐ = 10 .^rₐ;
        
        for i ∈ 1:10
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
    
    rₐ = range(9,stop = 18,length = 10);
    fₐ = 10 .^rₐ;
    
    for i ∈ 7:10
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
    
        μ,σ = 0, 1e-10;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
    
    
        t₁ = PQrun_2D!(N, t₀, Δt * 100, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)
        #PQplotting_2D!(N,t₀,t₁,t₂₀,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
    end
    
    end
end


for n ∈ 9:9

    @time begin
    b = 6.68
    N = 2^n;
    println(N)
    
    rₐ = range(9,stop = 18,length = 10);
    fₐ = 10 .^rₐ;
    
    for i ∈ 1:10
        t₀ = 3.01e-7 * (fₐ[i]/1e12)^(4/(4+b));
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
    
        μ,σ = 0, 0.1;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
    
    
        t₁ = EQCDrun_2D!(N, t₀, Δt * 100, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)
        #PQplotting_2D!(N,t₀,t₁,t₂₀,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
    end
    
    end
end




for n ∈ 9:9

    @time begin
    b = 6.68
    N = 2^n;
    println(N)
    
    rₐ = range(9,stop = 18,length = 10);
    fₐ = 10 .^rₐ;
    
    for i ∈ 1:10
        t₀ = 1.61e-10 * (fₐ[i]/1e12)
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
        Ȧ = zeros(N,N);
        
    
        μ,σ = 0, 1e-10;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
    
        t₁ = Lrun_2D!(N, t₀, Δt * 50, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)
        angle = zeros(N,N);
        angler!(angle,C₁,C₂);
        angler!(Ȧ,Ȧ₁,Ȧ₂)
        Lplotting_2D!(N, t₀,t₁, Δt * 100, angle, Ȧ, Δx, Δt,fₐ[i],i)
        #PQplotting_2D!(N,t₀,t₁,t₂₀,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
    end
    
    end
end



for n ∈ 9:9

    @time begin
    b = 6.68
    N = 2^n;
    println(N)
    
    rₐ = range(9,stop = 18,length = 10);
    fₐ = 10 .^rₐ;
    
    for i ∈ 1:10
        t₀ = 1.61e-10 * (fₐ[i]/1e12)
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ = zeros(N,N);

    
        μ,σ = 0, 0.1;
        C = rand(Normal(μ,σ),N,N);
    
    
        t₁ = Lrun_2D!(N, t₀, Δt * 50, C ,Ȧ, Δx, Δt,fₐ[i],i)

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
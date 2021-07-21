
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
# Pkg.add("GR")

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
using JLD
stats = pyimport("scipy.stats")
np = pyimport("numpy")
cv2 = pyimport("cv2")

x = randn(64,64,64);
y = randn(512,512);
@time x = x.+1;
@time y= y.+1;


η = 1;
ω = 5.0;
N = 2^n;





@time begin
    n = 8
    t₀ = 1;
    N = 2^n;
    println(N)

    Δx = 3
    Δt = 1


    t = 0.5N*Δx/Δt;
    t₀ = 0.1 ;
    t₂ = t/2 ;
    t₅ = t/5 ;
    t₁₀ = t/10 ;
    t₂₀ = t/20 ;
    t₅₀ = t/50 ;


    sw = 1
    

    Ȧ₁ = zeros(N,N);
    Ȧ₂ = zeros(N,N);

    μ,σ = 0, 0.1;
    C₁ = rand(Normal(μ,σ),N,N);
    C₂ = rand(Normal(μ,σ),N,N);
    
    i = 1

    t₁ = Gorg_2D!(N,t₀,t₅,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,sw,i)

end
        




for n ∈ 9:9

    @time begin
        t₀ = 1;
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
        
            Δx = 1
            Δt = 0.1
        
            t₁ = run_2D!(N, t₀, 250, C₁, C₂, Ȧ₁, Ȧ₂,ω,η, Δx, Δt,i)
            #plotting_2D!(N,t₀,t₁,250,C₁,C₂,Ȧ₁,Ȧ₂,ω,η, Δx,Δt,fₐ[i],i)
            
        end
        
    end
end





for n ∈ 9:9
    
    @time begin
    t₀ = 1;
    N = 2^n;
    println(N)
    
    fₐ = 10 .^15;
    Δx = 1;
    Δt = 0.1;
    

    
    println(fₐ)

    Ȧ₁ = zeros(N,N);
    Ȧ₂ = zeros(N,N);

    μ,σ = 0, 1e-10;
    C₁ = rand(Normal(μ,σ),N,N);
    C₂ = rand(Normal(μ,σ),N,N);



    t₁ = PQrun_2D!(N, t₀, 1000*Δt, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ,1)
    #PQplotting_2D!(N,t₀,t₁,t₂₀,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)

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
    s = 500
    
    
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
    
        t₁ = EQCDrun_2D!(N, t₀, s*Δt, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)

        Δx = 1
        Δt = 1e-6

        EQCDplotting_2D!(N,s/10,t₁,5000Δt,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
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
    s = 500;
    
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
    
    
        t₁ = Lrun_2D!(N, t₀, s*Δt, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)
        angle = zeros(N,N);
        angler!(angle,C₁,C₂);
        angler!(Ȧ,Ȧ₁,Ȧ₂)

        Δx = 1;
        Δt = 1e-8;
        #t₁ = 1.61e-10 * (fₐ[i]/1e12)
        Lplotting_2D!(N, s/10,t₁, 5000Δt, angle, Ȧ, Δx, Δt,fₐ[i],i)
        #PQplotting_2D!(N,t₀,t₁,t₂₀,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i)
    
    end
    
    end
end




for n ∈ 3:3

    @time begin
    
    N = 2^n;
    println(N)
    
    rₐ = range(15,stop = 24,length = 10);
    fₐ = 10 .^rₐ;
    r = 1
    s = 1
    
    
    for i ∈ 1:1
        
        t₀ = 0.0001;
        print(string(i)*" ")
        println(fₐ[i])
    
        Ȧ₁ = zeros(N,N);
        Ȧ₂ = zeros(N,N);
        Ȧ = zeros(N,N);
    
        μ,σ = 0, 0.1;
        C₁ = rand(Normal(μ,σ),N,N);
        C₂ = rand(Normal(μ,σ),N,N);
    
        Δx = round(8000/N,digits = 2)
        Δt = 0.0004
    
        FPQrun_2D!(N, t₀, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)

        
        save("PQend.jld", "Real", C₁,"Imaginary",C₂)
    
    
    end
    
    end
end


C₁ = load("PQend.jld")["Real"]
C₂ = load("PQend.jld")["Imaginary"]

t₁ = 0.4

Δx = round(4/N,digits = 5)
Δt = 0.001 



t₂ = FErun_2D!(N,t₁,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i,r,s)

angle = zeros(N,N);
angler!(angle,C₁,C₂);
angler!(Ȧ,Ȧ₁,Ȧ₂)


FLrun_2D!(N, t₂, 10, angle, Ȧ, Δx, Δt,fₐ[i],i,s)






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
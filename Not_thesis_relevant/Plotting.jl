
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
# Pkg.add("Statistics")

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
using Statistics
stats = pyimport("scipy.stats")
np = pyimport("numpy")
cv2 = pyimport("cv2")

x = randn(64,64,64);
y = randn(512,512);
@time x = x.+1;
@time y= y.+1;





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
    
η = 1;
ω = 5.0;
N = 2^n;



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


@time begin
    for n ∈ 9:9
    
        
        
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
    
        
            μ,σ = 0, 0.1;
            C₁ = rand(Normal(μ,σ),N,N);
            C₂ = rand(Normal(μ,σ),N,N);
        
      
            Δx = round(8000/N,digits = 2)
            Δt = 0.004*10
        
            FPQrun_2D!(n, N, t₀, C₁, C₂, Ȧ₁, Ȧ₂, Δx, Δt,fₐ[i],i)
    
      
    
        end
        
        
    end
        println(' ')
    end

saxion = load("PQStatstest9.jld")["saxion"]
axion = load("PQStatstest9.jld")["axion"]
tracker = load("PQStringstest9.jld")["time"]
axenergy = load("PQStatstest9.jld")["axenergy"]

n = 9

N = 2^n

k_freq = fftfreq(N)*N
kx,ky = meshgrid(k_freq,k_freq)

knrm = sqrt.( kx.^2 + ky.^2)
knrm = collect(Iterators.flatten(knrm))

kbins = range(0.5, N/2+1, step = 1)
kvals = 0.5 * (kbins[2:end] + kbins[1:end-1])

Plots.plot(kvals,axion[120],xaxis= :log,yaxis =:log,legend = false)
Plots.plot(tracker,axenergy)




tracker = load("PQStrings.jld")["time"]
cores = load("PQStrings.jld")["number"]

using Plots

Plots.plot(tracker,cores)


@time begin

    n = 9
    N = 2^n
    
    
    Δx = round(4/N,digits = 5)
    Δt = 0.001 * 10
    
    rₐ = range(15,stop = 24,length = 10);
    fₐ = 10 .^rₐ;
    
    i = 1
    for r ∈ 1:1
        for s ∈ 1:1
            C₁ = load("Saving/PQ/PQend"*string(n)*".jld")["Real"]
            C₂ = load("Saving/PQ/PQend"*string(n)*".jld")["Imaginary"]
            Ȧ₁ = load("Saving/PQ/PQend"*string(n)*".jld")["Realvel"]
            Ȧ₂ = load("Saving/PQ/PQend"*string(n)*".jld")["RealIm"]
    
            t₁ = 0.4
    
    
            FErun_2D!(n,N,t₁,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i,r,s)
    
        
        end
    end
        
    # for r ∈ 2:5
    #     for s ∈ 1:1
    #         C₁ = load("PQend"*string(n)*".jld")["Real"]
    #         C₂ = load("PQend"*string(n)*".jld")["Imaginary"]
    #         Ȧ₁ = load("PQend"*string(n)*".jld")["Realvel"]
    #         Ȧ₂ = load("PQend"*string(n)*".jld")["RealIm"]
    
    #         t₁ = 0.4
    
    
    #         FErun_2D!(N,t₁,C₁,C₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ[i],i,r,s)
    
    
    #     end
    # end
        
    println(' ')
    
    end




C₁ = load("EQCDend911.jld")["Real"]
C₂ = load("EQCDend911.jld")["Imaginary"]
angle = zeros(2^9,2^9)
angler!(angle,C₁,C₂)
PyPlot.imsave("E.png",angle,vmin=-π,vmax = π,cmap = "twilight")



tracker = load("EQCDStrings911.jld")["time"]
cores = load("EQCDStrings911.jld")["number"]

using Plots

Plots.plot(tracker,cores)





@time begin

    n = 9
    N = 2^n
    
    
    
    Δx = round(4/N,digits = 5)
    Δt = 0.001 
    
    rₐ = range(15,stop = 24,length = 10);
    fₐ = 10 .^rₐ;
    
    i = 1
    
    for s ∈ 1:1
    
        t₂ = 7
    
        C₁ = load("EQCDend911.jld")["Real"]
        C₂ = load("EQCDend911.jld")["Imaginary"]
        Ȧ₁ = load("EQCDend911.jld")["Realvel"]
        Ȧ₂ = load("EQCDend911.jld")["RealIm"]
    
        angle = zeros(N,N);
        Ȧ = zeros(N,N)
    
        angler!(angle,C₁,C₂);
        angler!(Ȧ,Ȧ₁,Ȧ₂)
    
    
        tracker,cores = FLrun_2D!(N, t₂, 10, angle, Ȧ, Δx, Δt,fₐ[i],i,s)
    
    
    
        save("LateStrings"*string(n)*string(s)*".jld","time",tracker,"number",cores)
    
    end
    
    
    end


    

tracker = load("LateStrings91.jld")["time"]
cores = load("LateStrings91.jld")["number"]


using Plots

Plots.plot(tracker,cores)







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

Plotting


Pkg.add("BenchmarkTools")
using BenchmarkTools

Pkg.add("Distributions")
using Distributions
using Random

using Pkg

Pkg.add("CxxWrap")

using CxxWrap

Pkg.add("JLD")

using JLD

Pkg.instantiate()
Pkg.add("ImageIO")
using ImageIO

using Plots
using PyPlot

using Images

function Laplacian(A,Δx)
    😄 = length(A[1,:]);
    P = zeros(😄,😄);
    for i = 1:😄
        for j = 1:😄
            P[i,j] = (-A[mod(i+2-1,😄)+1,j] + 16A[mod(i+1-1,😄)+1,j] 
            + 16A[mod(i-1-1,😄)+1,j] -A[mod(i-2-1,😄)+1,j] 
            -A[i,mod(j+2-1,😄)+1] + 16A[i,mod(j+1-1,😄)+1] 
            + 16A[i,mod(j-1-1,😄)+1] -A[i,mod(j-2-1,😄)+1] 
            - 60A[i,j]) /(12Δx^2);
        end
    end
    return P
end


#a(t) = 4.36000000000006e-18t - 6.78288102293483e-23
#F₁ = Laplacian(A₁,Δx) - a(time)^β * λ * A₁ .* (A₁.^2 + A₂.^2 - η^2) - α * © * Ȧ₁ / time
#F₂ = Laplacian(A₂,Δx) - a(time)^β * λ * A₂ .* (A₁.^2 + A₂.^2 - η^2) - α * © * Ȧ₂ / time

#B₁ = A₁ + Δt * (Ȧ₁ + 0.5Δt * F₁)
#B₂ = A₂ + Δt * (Ȧ₂ + 0.5Δt * F₂)

function update(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time)
    #radiation dominated

    function F(Aₓ,Aₛ,Ȧ,t)
        λ = 2π^2/ω^2;
        β = 0.0;
        α = 3.0;
        © = 1.0;
        r = Laplacian(Aₓ,Δx) .- (4.36000000000006e-18t - 6.78288102293483e-23)^β * 
         λ * Aₓ .* (Aₓ.^2 .+ Aₛ.^2 .- η^2) .- α * © * Ȧ / t;
        return r
    end

    

    F₁ = F(A₁,A₂,Ȧ₁,time);
    F₂ = F(A₂,A₁,Ȧ₂,time);

    time += Δt;

    function G(A,Ȧ,F)
        r = A .+ Δt * (Ȧ .+ 0.5Δt * F);
        return r
    end

    A₁ = G(A₁,Ȧ₁,F₁); #G updates A
    A₂ = G(A₂,Ȧ₂,F₂);

    Ȧ₁ = Ȧ₁ .+ 0.5Δt * (F₁ .+ F(A₁,A₂,Ȧ₁,time) );
    Ȧ₂ = Ȧ₂ .+ 0.5Δt * (F₂ .+ F(A₂,A₁,Ȧ₂,time) );

    return A₁,A₂,Ȧ₁,Ȧ₂,time
end

function update_2(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time)
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = Laplacian(A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = Laplacian(A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ Laplacian(B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ Laplacian(B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)

    return B₁,B₂,Ḃ₁,Ḃ₂,lime

end

function run(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    steps = round(t/Δt,digits = 0);
    
    for i in 1:steps
        time = round(time,digits = 1);
        A₁,A₂,Ȧ₁,Ȧ₂,time = update(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time);
    end

    return A₁,A₂,Ȧ₁,Ȧ₂,time
end



function saving(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    steps = round(t/Δt,digits = 0);


    for i in 1:steps
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2)
            mod(mod .> 1) .= 1
            #save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
            save("gray.png",colorview(Gray,mod))
        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_2(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end


n = 5;
Δx = 1;
Δt = 0.1;
η = 1;
ω = 5;
N = 2^n;
t = 0.5N*Δx/Δt;

Ȧ₁ = zeros(N,N);
Ȧ₂ = zeros(N,N);

μ,σ = 0, 0.1;
A₁ = rand(Normal(μ,σ),N,N);
A₂ = rand(Normal(μ,σ),N,N);
#sparse array
t₀ = 0.1 ;
t₂ = t/2 ;
t₅ = t/5 ;
t₁₀ = t/10 ;

print(t)

A₁,A₂ = saving(N, t₀, t, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt) ;




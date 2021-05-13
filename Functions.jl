function Laplacian(😄,A,Δx)
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

function Laplacian_circle(😄,A,Δx)
    B = CircularArray(A);
    P = zeros(😄,😄);
    for i = 1:😄
        for j = 1:😄
            P[i,j] = (-B[i+2,j] + 16B[i+1,j] 
            + 16B[i-1,j] -B[i-2,j] 
            -B[i,j+2] + 16B[i,j+1] 
            + 16B[i,j-1] -B[i,j-2] 
            - 60A[i,j]) /(12Δx^2);
        end
    end
    return P
end

function Laplacian_3D(😄::Int16,A::AbstractArray{Float16},Δx::Float16,P::AbstractArray{Float16})

    B = CircularArray(A);
    for i = 1:😄
        for j = 1:😄
            for k = 1:😄
                P[i,j,k] = (-B[i+2,j,k] + 16B[i+1,j,k] + 16B[i-1,j,k] - B[i-2,j,k]
                -B[i,j+2,k] + 16B[i,j+1,k] + 16B[i,j-1,k] - B[i,j-2,k]
                -B[i,j,k+2] + 16B[i,j,k+1] + 16B[i,j,k-1] - B[i,j,k-2] - 90B[i,j,k])/(12Δx^2)
            end
        end
    end
    return P
end

function Laplacian_roll(😄,A,Δx)
    P = zeros(😄,😄)
    P = (-circshift(A,(-2,0)) + 16circshift(A,(-1,0)) + 16circshift(A,(1,0)) - circshift(A,(2,0))
    -circshift(A,(0,-2)) + 16circshift(A,(0,-1)) + 16circshift(A,(0,1)) - circshift(A,(0,2)) - 60A)/(12Δx^2)
end


function F(Aₓ,Aₛ,Ȧ,t)
    λ = 2π^2/ω^2;
    β = 0.0;
    α = 3.0;
    © = 1.0;
    r = Laplacian(N,Aₓ,Δx) .- (4.36000000000006e-18t - 6.78288102293483e-23)^β * 
     λ * Aₓ .* (Aₓ.^2 .+ Aₛ.^2 .- η^2) .- α * © * Ȧ / t;
    return r
end

function G(A,Ȧ,F)
    r = A .+ Δt * (Ȧ .+ 0.5Δt * F);
    return r
end


function update(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time)
    #radiation dominated

    F₁ = F(A₁,A₂,Ȧ₁,time);
    F₂ = F(A₂,A₁,Ȧ₂,time);

    time += Δt;

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

    F₁ = Laplacian(N,A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = Laplacian(N,A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ Laplacian(N,B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ Laplacian(N,B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)

    return B₁,B₂,Ḃ₁,Ḃ₂,lime

end

function update_3(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time)
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = Laplacian_roll(N,A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = Laplacian_roll(N,A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ Laplacian_roll(N,B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ Laplacian_roll(N,B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)

    return B₁,B₂,Ḃ₁,Ḃ₂,lime

end

function update_4(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time)
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = Laplacian_circle(N,A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = Laplacian_circle(N,A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ Laplacian_circle(N,B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ Laplacian_circle(N,B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)

    return B₁,B₂,Ḃ₁,Ḃ₂,lime

end

function update_6(N,A₁::AbstractArray{Float16},A₂::AbstractArray{Float16},Ȧ₁::AbstractArray{Float16},Ȧ₂::AbstractArray{Float16},ω,η,Δx,Δt,time,P::AbstractArray{Float16})
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = Laplacian_3D(N,A₁,Δx,P) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = Laplacian_3D(N,A₂,Δx,P) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ Laplacian_3D(N,B₁,Δx,P) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ Laplacian_3D(N,B₂,Δx,P) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)

    return B₁,B₂,Ḃ₁,Ḃ₂,lime

end

function update_7(N,A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = M₁ .- α * © .* Ȧ₁ ./ time
    F₂ = M₂ .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    L₁ = Laplacian_circle(N,B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2)
    L₂ = Laplacian_circle(N,B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ L₁ .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ L₂ .- α * © .* Ȧ₂ ./ lime)

    return B₁,B₂,Ḃ₁,Ḃ₂,L₁,L₂,lime

end



function run(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    steps = round(t/Δt,digits = 0);
    
    for _ ∈ 1:steps
        time = round(time,digits = 1);
        A₁,A₂,Ȧ₁,Ȧ₂,time = update(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time);
    end

    return A₁,A₂,Ȧ₁,Ȧ₂,time
end

function saving_function(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for _ ∈ 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));

        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for _ ∈ 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_2(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_roll(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for _ ∈ 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));

        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_3(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_circle(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for _ ∈ 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));

        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_4(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_laplace(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);

    M₁ = Laplacian_circle(N,A₁,Δx)
    M₂ = Laplacian_circle(N,A₂,Δx)

    for _ ∈ 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));

        end
    A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,time = update_7(N, A₁, A₂, Ȧ₁, Ȧ₂, M₁, M₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_3D(N::Int16,t₀::Float16,t,A₁::AbstractArray{Float16},A₂::AbstractArray{Float16},Ȧ₁::AbstractArray{Float16},Ȧ₂::AbstractArray{Float16},ω::Float16,η::Float16,Δx::Float16,Δt::Float16)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);

    P = zeros(N,N,N)
    P = convert(Array{Float16},P)
    for _ ∈ 1:steps
        time = round(time,digits = 1);
    #     # if time % 1 == 0

    #     #     #save("plottting_3D/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,sqrt.(A₁.^2 .+ A₂.^2)[sqrt.(A₁.^2 .+ A₂.^2) .> 1]));
    #     #     #save("gray.png",colorview(Gray,mod))
    #     # end
        A₁,A₂,Ȧ₁,Ȧ₂,time = update_6(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time,P);
        
    end
    return A₁,A₂
end

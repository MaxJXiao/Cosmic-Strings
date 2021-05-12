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

function Laplacian_tensor(😄,A,Δx)
    B = CircularArray(A);
    P = zeros(😄,😄);
    @tullio P[i,j] := (-B[i+2,j] + 16B[i+1,j] 
    + 16B[i-1,j] -B[i-2,j] 
    -B[i,j+2] + 16B[i,j+1] 
    + 16B[i,j-1] -B[i,j-2] 
    - 60A[i,j]) /(12Δx^2);
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

function Laplacian_3D(😄,A,Δx)
    B = CircularArray(A);
    P = zeros(😄,😄,😄);
    for i = 1: 😄
        for j = 1: 😄
            for k = 1: 😄
                P[i,j,k] = (-B[i+2,j,k] + 16B[i+1,j,k] + 16B[i-1,j,k] - B[i-2,j,k]
                -B[i,j+2,k] + 16B[i,j+1,k] + 16B[i,j-1,k] - B[i,j-2,k]
                -B[i,j,k+2] + 16B[i,j,k+1] + 16B[i,j,k-1] - B[i,j,k-2] - 90A[i,j,k])/(12Δx^2)
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

function update_5(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time)
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = Laplacian_tensor(N,A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = Laplacian_tensor(N,A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ Laplacian_tensor(N,B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ Laplacian_tensor(N,B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)

    return B₁,B₂,Ḃ₁,Ḃ₂,lime

end

function update_6(N,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,time)
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = Laplacian_3D(N,A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = Laplacian_3D(N,A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time
    
    lime = time + Δt

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ Laplacian_3D(N,B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ Laplacian_3D(N,B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)

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

function saving_function(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for i in 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
            #save("gray.png",colorview(Gray,mod))
        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for i in 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
            #save("gray.png",colorview(Gray,mod))
        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_2(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_roll(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for i in 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
            #save("gray.png",colorview(Gray,mod))
        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_3(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_circle(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for i in 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
            #save("gray.png",colorview(Gray,mod))
        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_4(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_tensor(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for i in 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
            #save("gray.png",colorview(Gray,mod))
        end
    A₁,A₂,Ȧ₁,Ȧ₂,time = update_5(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

function saving_3D(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);


    for i in 1:steps
        time = round(time,digits = 1);
    #     # if time % 1 == 0
    #     #     mod = sqrt.(A₁.^2 .+ A₂.^2);
    #     #     mod[mod .> 1] .= 1;
    #     #     #save("plottting_3D/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));
    #     #     #save("gray.png",colorview(Gray,mod))
    #     # end
        A₁,A₂,Ȧ₁,Ȧ₂,time = update_6(N, A₁, A₂, Ȧ₁, Ȧ₂, ω, η, Δx, Δt, time);
        
    end
    return A₁,A₂
end

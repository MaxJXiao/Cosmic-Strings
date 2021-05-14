function Laplacian_circle(😄,A,Δx)
    B = CircularArray(A);
    for i ∈ 1: 😄
        for j ∈ 1: 😄
            @fastmath P[i,j] = (-B[i+2,j] + 16B[i+1,j] 
                + 16B[i-1,j] -B[i-2,j] 
                -B[i,j+2] + 16B[i,j+1] 
                + 16B[i,j-1] -B[i,j-2] 
                - 60A[i,j]) /(12Δx^2);
        end
    end
    return P
end

function Laplacian_circle!(P,A,Δx)
    B = CircularArray(A);
    for 😄 ∈ CartesianIndices(P)
        (i,j) = Tuple(😄)
        @fastmath P[i,j] = (-B[i+2,j] + 16B[i+1,j] 
            + 16B[i-1,j] -B[i-2,j] 
            -B[i,j+2] + 16B[i,j+1] 
            + 16B[i,j-1] -B[i,j-2] 
            - 60A[i,j]) /(12Δx^2);
    end
    return nothing
end


function update_4(N,A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
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


function update_7!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = M₁ .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ time
    F₂ = M₂ .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ time

    B₁ = A₁ .+ Δt * (Ȧ₁ .+ 0.5Δt * F₁)
    B₂ = A₂ .+ Δt * (Ȧ₂ .+ 0.5Δt * F₂)

    lime = time + Δt

    Laplacian_circle!(M₁,B₁,Δx)
    Laplacian_circle!(M₂,B₂,Δx)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2) .- α * © .* Ȧ₂ ./ lime)


    return B₁,B₂,Ḃ₁,Ḃ₂,lime

end



function saving_circle(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    time = round(time,digits = 1)

    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0

    λ = 2π^2/ω^2
    
    steps = round(t/Δt,digits = 0);
    M₁ = Laplacian_circle(N,A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2)
    M₂ = Laplacian_circle(N,A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2)


    for _ ∈ 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));

        end
        A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,time = update_4(N, A₁, A₂, Ȧ₁, Ȧ₂,M₁,M₂, ω, η, Δx, Δt, time);
    end

    
    return A₁,A₂

end

function saving_laplace!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    
    steps = round(t/Δt,digits = 0);

    M₁ = zeros(N,N)
    M₂ = zeros(N,N)

    Laplacian_circle!(M₁,A₁,Δx) 
    Laplacian_circle!(M₂,A₂,Δx)

    for _ ∈ 1:steps
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plotting/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));

        end
        A₁,A₂,Ȧ₁,Ȧ₂,time = update_7!(A₁, A₂, Ȧ₁, Ȧ₂, M₁, M₂, ω, η, Δx, Δt, time);
    end
    return nothing
end

function Laplacian_2D!(P,A,Δx)
    B = CircularArray(A);
    for 😄 ∈ CartesianIndices(P)
        (i,j) = Tuple(😄)
        @fastmath P[i,j] = (-B[i+2,j] + 16B[i+1,j] 
            + 16B[i-1,j] -B[i-2,j] 
            -B[i,j+2] + 16B[i,j+1] 
            + 16B[i,j-1] -B[i,j-2] 
            - 60A[i,j]) /(12Δx^2);
    end
    return nothing
end

function update_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
    a = 4.36000000000006e-18*time - 6.78288102293483e-23
    a₁ = 4.36000000000006e-18*(time + Δt) - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    F₂ = M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    Laplacian_2D!(M₁,A₁,Δx)
    Laplacian_2D!(M₂,A₂,Δx)

    Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    return nothing
end

function plotting_2D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)

    time = t₀

    M₁ = zeros(N,N)
    M₂ = zeros(N,N)

    Laplacian_2D!(M₁,A₁,Δx)
    Laplacian_2D!(M₂,A₂,Δx)

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if time % 1 == 0
            mod = sqrt.(A₁.^2 .+ A₂.^2);
            mod[mod .> 1] .= 1;
            save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,mod));

        end
        update_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
        time = time + Δt

    end

    return nothing
end




function Laplacian_3D!(P,B,Δx)
    A = CircularArray(B);
    for 😄 ∈ CartesianIndices(P)
        (i,j,k) = Tuple(😄)
        @fastmath P[i,j,k] = (A[i-1,j-1,k-1] + 3A[i-1,j,k-1] + A[i-1,j+1,k-1] 
                + 3A[i-1,j-1,k] + 14A[i-1,j,k] + 3A[i-1,j+1,k] 
                + A[i-1,j-1,k+1] + 3A[i-1,j,k+1] + A[i-1,j+1,k+1] 
                + 3A[i,j-1,k-1] + 14A[i,j,k-1] + 3A[i,j+1,k-1] 
                + 14A[i,j-1,k] - 128A[i,j,k] + 14A[i,j+1,k] 
                + 3A[i,j-1,k+1] + 14A[i,j,k+1] + 3A[i,j+1,k+1] 
                + A[i+1,j-1,k-1] + 3A[i+1,j,k-1] + A[i+1,j+1,k-1] 
                + 3A[i+1,j-1,k] + 14A[i+1,j,k] + 3A[i+1,j+1,k] 
                + A[i+1,j-1,k+1] + 3A[i+1,j,k+1] + A[i+1,j+1,k+1])/(30Δx^2)
    end
    return nothing
end

function update_3D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
    a = 4.36000000000006e-18*time - 6.78288102293483e-23
    a₁ = 4.36000000000006e-18*(time + Δt) - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    F₁ = M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    F₂ = M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    Laplacian_3D!(M₁,A₁,Δx)
    Laplacian_3D!(M₂,A₂,Δx)

    Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    return nothing
end

function plotting_3D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)

    time = t₀

    M₁ = zeros(N,N,N)
    M₂ = zeros(N,N,N)

    Laplacian_3D!(M₁,A₁,Δx)
    Laplacian_3D!(M₂,A₂,Δx)

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        update_3D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
        time = time + Δt
    end

    return nothing
end



function Laplacian_slow3D(😄,B,Δx)
    A = CircularArray(B);
    for i ∈ 1:😄
        for j ∈ 1:😄
            for k ∈ 1:😄
                @fastmath P[i,j,k] = (A[i-1,j-1,k-1] + 3A[i-1,j,k-1] + A[i-1,j+1,k-1] 
                + 3A[i-1,j-1,k] + 14A[i-1,j,k] + 3A[i-1,j+1,k] 
                + A[i-1,j-1,k+1] + 3A[i-1,j,k+1] + A[i-1,j+1,k+1] 
                + 3A[i,j-1,k-1] + 14A[i,j,k-1] + 3A[i,j+1,k-1] 
                + 14A[i,j-1,k] - 128A[i,j,k] + 14A[i,j+1,k] 
                + 3A[i,j-1,k+1] + 14A[i,j,k+1] + 3A[i,j+1,k+1] 
                + A[i+1,j-1,k-1] + 3A[i+1,j,k-1] + A[i+1,j+1,k-1] 
                + 3A[i+1,j-1,k] + 14A[i+1,j,k] + 3A[i+1,j+1,k] 
                + A[i+1,j-1,k+1] + 3A[i+1,j,k+1] + A[i+1,j+1,k+1])/(30Δx^2)
            end
        end
    end
    return P
end


function update_slow3D(N,A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
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

    L₁ = Laplacian_slow3D(N,B₁,Δx) .- a(lime).^β * λ .* B₁ .* (B₁.^2 .+ B₂.^2 .- η^2)
    L₂ = Laplacian_slow3D(N,B₂,Δx) .- a(lime).^β * λ .* B₂ .* (B₁.^2 .+ B₂.^2 .- η^2)

    Ḃ₁ = Ȧ₁ .+ 0.5Δt .* (F₁ .+ L₁ .- α * © .* Ȧ₁ ./ lime)
    Ḃ₂ = Ȧ₂ .+ 0.5Δt .* (F₂ .+ L₂ .- α * © .* Ȧ₂ ./ lime)
    
    return B₁,B₂,Ḃ₁,Ḃ₂,L₁,L₂,lime

end


function plotting_slow3D(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)
    time = t₀;
    time = round(time,digits = 1)

    a(t) = 4.36000000000006e-18*t - 6.78288102293483e-23

    β = 0

    λ = 2π^2/ω^2

    M₁ = Laplacian_slow3D(N,A₁,Δx) .- a(time).^β * λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η^2)
    M₂ = Laplacian_slow3D(N,A₂,Δx) .- a(time).^β * λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η^2)


    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,time = update_slow3D(N,A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,ω,η,Δx,Δt,time)
    end

    
    return A₁,A₂

end

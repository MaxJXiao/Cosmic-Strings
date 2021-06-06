#Early Times in QCD Epoch
#T ∼ Λ 
#High Temperature mass of the axion is neglected.
#λ = 1, free parameter

function Laplacian_2D!(P₁,P₂,A₁,A₂,Δx)
    B₁ = CircularArray(A₁);
    B₂ = CircularArray(A₂);
    Threads.@threads for 😄 ∈ CartesianIndices(P₁)
        (i,j) = Tuple(😄)
        P₁[i,j] = @fastmath (-B₁[i+2,j] + 16B₁[i+1,j] 
            + 16B₁[i-1,j] -B₁[i-2,j] 
            -B₁[i,j+2] + 16B₁[i,j+1] 
            + 16B₁[i,j-1] -B₁[i,j-2] 
            - 60B₁[i,j]) /(12Δx^2);
        P₂[i,j] = @fastmath (-B₂[i+2,j] + 16B₂[i+1,j] 
            + 16B₂[i-1,j] -B₂[i-2,j] 
            -B₂[i,j+2] + 16B₂[i,j+1] 
            + 16B₂[i,j-1] -B₂[i,j-2] 
            - 60B₂[i,j]) /(12Δx^2);  
    end
    return nothing
end


function EQCDfupdate_2D!(F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,λ,fₐ,ηₓ)
    n = 6.68;
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F₁)
        (i,j) = Tuple(😄)
        F₁[i,j] = @fastmath M₁[i,j] - λ * C₁[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + ηₓ^n *η^2 * C₂[i,j]^2/((A₁[i,j]^2 .+ A₂[i,j]^2)^1.5) - 2/η * Ȧ₁[i,j];
        F₂[i,j] = @fastmath M₂[i,j] - λ * C₂[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) - ηₓ^n *η^2 * C₁[i,j]*C₂[i,j]/((A₁[i,j]^2 .+ A₂[i,j]^2)^1.5) - 2/η * Ȧ₂[i,j];
    end
    return nothing
end


function EQCDAupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(A₁)
        (i,j) = Tuple(😄)
        A₁[i,j] = @fastmath A₁[i,j] .+ Δt .* (Ȧ₁[i,j] .+ 0.5Δt .* F₁[i,j])
        A₂[i,j] = @fastmath A₂[i,j] .+ Δt .* (Ȧ₂[i,j] .+ 0.5Δt .* F₂[i,j])
    end
    return nothing
end


function EQCDvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,η,λ,fₐ,ηₓ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ₁)
        (i,j) = Tuple(😄)
        Ȧ₁[i,j] = @fastmath Ȧ₁[i,j] .+ 0.5Δt .* (F₁[i,j] .+ M₁[i,j] - λ * C₁[i,j] * (η^2 .* (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + 8.4e5 * 1e12/(3fₐ)) - 2/η * Ȧ₁[i,j])
        Ȧ₂[i,j] = @fastmath Ȧ₂[i,j] .+ 0.5Δt .* (F₂[i,j] .+ M₂[i,j] - λ * C₂[i,j] * (η^2 .* (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + 8.4e5 * 1e12/(3fₐ)) - 2/η * Ȧ₂[i,j])
    end
    return nothing
end

function mass(fₐ,T)
    fₐ = fₐ*1e3
    n = 6.68;
    Λ = 400;
    αₐ = 1.68e-7;
    mᵤ = 1.7; #1.7 - 3.3MeV
    m₍d₎ = 4.1; #4.1 - 5.8MeV
    m₍π₎ = 135;
    f₍π₎ = 130;

    mass = αₐ * Λ^(4+n) / (fₐ^2 * T^n)
    mₐ = sqrt( m₍π₎^2 * f₍π₎^2 / fₐ^2 * mᵤ * m₍d₎ / (mᵤ + m₍d₎)^2 )
    if mass > mₐ
        mass = mₐ
    end
    return mass
end

function ηtime(time,fₐ,T)
    n = 6.68
    if T > 103
        t₁ = 3.01e-7 * (fₐ/1e12)^(4/(4+n));
    else
        t₁ = 1.61e-10 * (fₐ/1e12)^(4/(4+n));
    end
    η = (time/t₁)^0.5;

    ηₓ = 2e3/1.15e2;

    if ηₓ > η
        ηₓ = η
    end

    return ηₓ,η
end

function EQCDupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₀,time,fₐ)

    λ = (fₐ/mass(fₐ,400))^2;
    ηₓ,η = ηtime(time,fₐ,400);
    

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    EQCDfupdate_2D!(F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,Ȧ₁,Ȧ₂,η,λ,fₐ,ηₓ)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    EQCDAupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    ηₓ,η = ηtime(time+Δt,fₐ,400);

    EQCDvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,η,λ,fₐ,ηₓ)

    return nothing
end

function EQCDrun_2D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)

    time = t₀

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);


    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    angle = zeros(N,N);
    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if time % 1 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
            angler!(angle,A₁,A₂);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            PyPlot.imsave("EQCD/"*string(i)*"/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end
        EQCDupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₀,time,fₐ)
        time = time + Δt

    end

    return time
end

function EQCDplotting_2D!(N,t₀,t₁,t,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)

    time = t₁

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    # moo = zeros(N,N);
    angle = zeros(N,N);
    angler!(angle,A₁,A₂);

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if time % 1 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
            angler!(angle,A₁,A₂);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            PyPlot.imsave("EQCD/"*string(i)*"/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end
        EQCDupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₁,time,fₐ)
        time = time + Δt

    end

    return nothing
end
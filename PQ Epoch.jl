#PQ Epoch
#T ≫ Λ 
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


function PQfupdate_2D!(F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,λ,fₐ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F₁)
        (i,j) = Tuple(😄)
        F₁[i,j] = @fastmath M₁[i,j] - λ * C₁[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + λ * C₁[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₁[i,j];
        F₂[i,j] = @fastmath M₂[i,j] - λ * C₂[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + λ * C₂[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₂[i,j];
    end
    return nothing
end


function PQAupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(A₁)
        (i,j) = Tuple(😄)
        A₁[i,j] = @fastmath A₁[i,j] .+ Δt .* (Ȧ₁[i,j] .+ 0.5Δt .* F₁[i,j])
        A₂[i,j] = @fastmath A₂[i,j] .+ Δt .* (Ȧ₂[i,j] .+ 0.5Δt .* F₂[i,j])
    end
    return nothing
end


function PQvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,η,λ,fₐ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ₁)
        (i,j) = Tuple(😄)
        Ȧ₁[i,j] = @fastmath Ȧ₁[i,j] .+ 0.5Δt .* (F₁[i,j] .+ M₁[i,j] - λ * C₁[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + λ * C₁[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₁[i,j]);
        Ȧ₂[i,j] = @fastmath Ȧ₂[i,j] .+ 0.5Δt .* (F₂[i,j] .+ M₂[i,j] - λ * C₂[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + λ * C₂[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₂[i,j]);
    end
    return nothing
end

function PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₀,time,fₐ)

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    PQfupdate_2D!(F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,Ȧ₁,Ȧ₂,(time/t₀)^0.5,1,fₐ)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    PQAupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    PQvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,((time + Δt)/t₀)^0.5,1,fₐ)

    return nothing
end

function PQrun_2D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)

    time = t₀

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);


    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    angle = zeros(N,N);

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 10);
        #if time % 1 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        angler!(angle,A₁,A₂);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
        PyPlot.imsave("PQEpoch/"*string(i)*"/"*lpad( string(trunc(Int,(time-t₀)/Δt)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        #end
        PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₀,time,fₐ)
        time = time + Δt

    end

    return time
end

function PQplotting_2D!(N,t₀,t₁,t,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)

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
            PyPlot.imsave("PQEpoch/"*string(i)*"/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end
        PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₁,time,fₐ)
        time = time + Δt

    end

    return nothing
end
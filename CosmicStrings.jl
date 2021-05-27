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

function fupdate_2D!(F₁,F₂,M₁,M₂,a,©,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,time,β,α,λ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F₁)
        (i,j) = Tuple(😄)
        F₁[i,j] = @fastmath M₁[i,j] - a^β * λ * C₁[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - α * © * Ȧ₁[i,j] /time;
        F₂[i,j] = @fastmath M₂[i,j] - a^β * λ * C₂[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - α * © * Ȧ₂[i,j] /time;
    end
    return nothing
end

function Aupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(A₁)
        (i,j) = Tuple(😄)
        A₁[i,j] = @fastmath A₁[i,j] .+ Δt .* (Ȧ₁[i,j] .+ 0.5Δt .* F₁[i,j])
        A₂[i,j] = @fastmath A₂[i,j] .+ Δt .* (Ȧ₂[i,j] .+ 0.5Δt .* F₂[i,j])
    end
    return nothing
end

function velupdate_2D!(Ȧ₁,Ȧ₂,Δt,a,©,F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,η,time,β,α,λ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ₁)
        (i,j) = Tuple(😄)
        Ȧ₁[i,j] = @fastmath Ȧ₁[i,j] .+ 0.5Δt .* (F₁[i,j] .+ M₁[i,j] .- a.^β .* λ .* C₁[i,j] .* (A₁[i,j].^2 .+ A₂[i,j].^2 .- η.^2) .- α .* © .* Ȧ₁[i,j] ./ time)
        Ȧ₂[i,j] = @fastmath Ȧ₂[i,j] .+ 0.5Δt .* (F₂[i,j] .+ M₂[i,j] .- a.^β .* λ .* C₂[i,j] .* (A₁[i,j].^2 .+ A₂[i,j].^2 .- η.^2) .- α .* © .* Ȧ₂[i,j] ./ time)
    end
    return nothing
end

function update_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,ω,η,Δx,Δt,time)
    a = 4.36000000000006e-18*time - 6.78288102293483e-23
    a₁ = 4.36000000000006e-18*(time + Δt) - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    fupdate_2D!(F₁,F₂,M₁,M₂,a,©,A₁,A₂,A₁,A₂,Ȧ₁,Ȧ₂,η,time,β,α,λ)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    Aupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    velupdate_2D!(Ȧ₁, Ȧ₂, Δt, a₁, ©, F₁, F₂, M₁, M₂, A₁, A₂, A₁, A₂, η, time + Δt, β, α, λ)

    return nothing
end

function setting!(moo)
    @inbounds Threads.@threads for t in eachindex(moo)
        if moo[t] > 1
            moo[t] = 1
        end
    end
    return nothing
end

function mooing!(moo,A₁,A₂)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(moo)
        (i,j) = Tuple(😄)
        @fastmath moo[i,j] = sqrt(A₁[i,j]^2 + A₂[i,j]^2)
    end
    return nothing
end

function angler!(angle,A₁,A₂)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(angle)
        (i,j) = Tuple(😄)
        angle[i,j] = @fastmath atan(A₂[i,j],A₁[i,j]);
    end
    return nothing
end

function plotting_2D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)

    time = t₀

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    # moo = zeros(N,N);
    # angle = zeros(N,N);

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        # if time % 1 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        #     angler!(angle,A₁,A₂);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
        #     PyPlot.imsave("plottting_angle/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        # end
        update_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,ω,η,Δx,Δt,time)
        time = time + Δt

    end

    return nothing
end




function Laplacian_3D!(P₁,P₂,B₁,B₂,Δx)
    A₁ = CircularArray(B₁);
    A₂ = CircularArray(B₂);
    Threads.@threads for 😄 ∈ CartesianIndices(P₁)
        (i,j,k) = Tuple(😄)
        P₁[i,j,k] = @fastmath (A₁[i-1,j-1,k-1] + 3A₁[i-1,j,k-1] + A₁[i-1,j+1,k-1] 
                + 3A₁[i-1,j-1,k] + 14A₁[i-1,j,k] + 3A₁[i-1,j+1,k] 
                + A₁[i-1,j-1,k+1] + 3A₁[i-1,j,k+1] + A₁[i-1,j+1,k+1] 
                + 3A₁[i,j-1,k-1] + 14A₁[i,j,k-1] + 3A₁[i,j+1,k-1] 
                + 14A₁[i,j-1,k] - 128A₁[i,j,k] + 14A₁[i,j+1,k] 
                + 3A₁[i,j-1,k+1] + 14A₁[i,j,k+1] + 3A₁[i,j+1,k+1] 
                + A₁[i+1,j-1,k-1] + 3A₁[i+1,j,k-1] + A₁[i+1,j+1,k-1] 
                + 3A₁[i+1,j-1,k] + 14A₁[i+1,j,k] + 3A₁[i+1,j+1,k] 
                + A₁[i+1,j-1,k+1] + 3A₁[i+1,j,k+1] + A₁[i+1,j+1,k+1])/(30Δx^2);
        P₂[i,j,k] = @fastmath (A₂[i-1,j-1,k-1] + 3A₂[i-1,j,k-1] + A₂[i-1,j+1,k-1] 
                + 3A₂[i-1,j-1,k] + 14A₂[i-1,j,k] + 3A₂[i-1,j+1,k] 
                + A₂[i-1,j-1,k+1] + 3A₂[i-1,j,k+1] + A₂[i-1,j+1,k+1] 
                + 3A₂[i,j-1,k-1] + 14A₂[i,j,k-1] + 3A₂[i,j+1,k-1] 
                + 14A₂[i,j-1,k] - 128A₂[i,j,k] + 14A₂[i,j+1,k] 
                + 3A₂[i,j-1,k+1] + 14A₂[i,j,k+1] + 3A₂[i,j+1,k+1] 
                + A₂[i+1,j-1,k-1] + 3A₂[i+1,j,k-1] + A₂[i+1,j+1,k-1] 
                + 3A₂[i+1,j-1,k] + 14A₂[i+1,j,k] + 3A₂[i+1,j+1,k] 
                + A₂[i+1,j-1,k+1] + 3A₂[i+1,j,k+1] + A₂[i+1,j+1,k+1])/(30Δx^2);
    end
    return nothing
end


function fupdate_3D!(F₁,F₂,M₁,M₂,a,©,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,time,β,α,λ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F₁)
        (i,j,k) = Tuple(😄)
        F₁[i,j,k] = @fastmath M₁[i,j,k] - a^β * λ * C₁[i,j,k] *(A₁[i,j,k]^2 .+ A₂[i,j,k]^2 .- η^2) - α * © * Ȧ₁[i,j,k] /time;
        F₂[i,j,k] = @fastmath M₂[i,j,k] - a^β * λ * C₂[i,j,k] *(A₁[i,j,k]^2 .+ A₂[i,j,k]^2 .- η^2) - α * © * Ȧ₂[i,j,k] /time;
    end
    return nothing
end

function Aupdate_3D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(A₁)
        (i,j,k) = Tuple(😄)
        A₁[i,j,k] = @fastmath A₁[i,j,k] .+ Δt .* (Ȧ₁[i,j,k] .+ 0.5Δt .* F₁[i,j,k])
        A₂[i,j,k] = @fastmath A₂[i,j,k] .+ Δt .* (Ȧ₂[i,j,k] .+ 0.5Δt .* F₂[i,j,k])
    end
    return nothing
end

function velupdate_3D!(Ȧ₁,Ȧ₂,Δt,a,©,F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,η,time,β,α,λ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ₁)
        (i,j,k) = Tuple(😄)
        Ȧ₁[i,j,k] = @fastmath Ȧ₁[i,j,k] .+ 0.5Δt .* (F₁[i,j,k] .+ M₁[i,j,k] .- a.^β .* λ .* C₁[i,j,k] .* (A₁[i,j,k].^2 .+ A₂[i,j,k].^2 .- η.^2) .- α .* © .* Ȧ₁[i,j,k] ./ time)
        Ȧ₂[i,j,k] = @fastmath Ȧ₂[i,j,k] .+ 0.5Δt .* (F₂[i,j,k] .+ M₂[i,j,k] .- a.^β .* λ .* C₂[i,j,k] .* (A₁[i,j,k].^2 .+ A₂[i,j,k].^2 .- η.^2) .- α .* © .* Ȧ₂[i,j,k] ./ time)
    end
    return nothing
end


function update_3D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,ω,η,Δx,Δt,time)
    a = 4.36000000000006e-18*time - 6.78288102293483e-23
    a₁ = 4.36000000000006e-18*(time + Δt) - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    fupdate_3D!(F₁,F₂,M₁,M₂,a,©,A₁,A₂,A₁,A₂,Ȧ₁,Ȧ₂,η,time,β,α,λ)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    Aupdate_3D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)

    Laplacian_3D!(M₁,M₂,A₁,A₂,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    velupdate_3D!(Ȧ₁, Ȧ₂, Δt, a₁, ©, F₁, F₂, M₁, M₂, A₁, A₂, A₁, A₂, η, time + Δt, β, α, λ)

    return nothing
end

function plotting_3D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt)

    time = t₀

    M₁ = zeros(N,N,N)
    M₂ = zeros(N,N,N)

    F₁ = zeros(N,N,N);
    F₂ = zeros(N,N,N);

    Laplacian_3D!(M₁,M₂,A₁,A₂,Δx)

    #moo = zeros(N,N,N);

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
 
        # if time % 5 == 0
        #     mooing(moo,A₁,A₂)
        #     mood = Tuple.(findall(<(0.5),moo));
        #     p = plot(mood,seriestype = :scatter,xlims = (0,N),ylims = (0,N),zlims = (0,N))
        #     display(p)
        # end
        update_3D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,ω,η,Δx,Δt,time)
        time = time + Δt
    end

    return nothing
end



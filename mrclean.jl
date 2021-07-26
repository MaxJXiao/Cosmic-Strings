function Laplacian_2D!(P₁,P₂,A₁,A₂,Δx)
    B₁ = CircularArrays.CircularArray(A₁);
    B₂ = CircularArrays.CircularArray(A₂);
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



function Cores_2D!(N,angle,thr)
    s = []
    count = 0
    accept = 0.5 - 0.5*thr/100
    for 😄 ∈ 1:(N-1)
        for 🥪 ∈ 1:(N-1)
            norm₁ = (angle[😄,🥪] + π)/(2π)
            norm₂ = (angle[😄+1,🥪] + π)/(2π)
            norm₃ = (angle[😄+1,🥪+1] + π)/(2π)
            norm₄ = (angle[😄,🥪+1] + π)/(2π)

            θ₁ = min(abs(norm₂ - norm₁), 1 - abs(norm₂ - norm₁))
            θ₂ = min(abs(norm₃ - norm₂), 1 - abs(norm₃ - norm₂))
            θ₃ = min(abs(norm₄ - norm₃), 1 - abs(norm₄ - norm₃))
            θₛ = θ₁ + θ₂ + θ₃

            if θₛ > accept 
                append!(s,[[😄,🥪]])
            end
        end
    end

    if length(s) > 0
        for 🇸🇦 ∈ 1:(length(s)-1)
  
            diffᵣ = s[🇸🇦 + 1][1] - s[🇸🇦][1]
            diffₛ = s[🇸🇦 + 1][2] - s[🇸🇦][2]

            if diffᵣ == 0 && diffₛ == 1
                count += 1
            end
            if diffᵣ == 1 && diffₛ == 0
                count += 1
            end
        end
    end

    num = length(s) - count

    return num
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

function meshgrid(xin,yin)
    nx=length(xin)
    ny=length(yin)
    xout=zeros(ny,nx)
    yout=zeros(ny,nx)
    for jx=1:nx
        for ix=1:ny
            xout[ix,jx]=xin[jx]
            yout[ix,jx]=yin[ix]
        end
    end
    return (x=xout, y=yout)
end


function PQfupdate_2D!(F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,λ,fₐ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F₁)
        (i,j) = Tuple(😄)
        F₁[i,j] = @fastmath M₁[i,j] - λ * C₁[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) - λ * C₁[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₁[i,j];
        F₂[i,j] = @fastmath M₂[i,j] - λ * C₂[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) - λ * C₂[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₂[i,j];
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
        Ȧ₁[i,j] = @fastmath Ȧ₁[i,j] .+ 0.5Δt .* (F₁[i,j] .+ M₁[i,j] - λ * C₁[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) - λ * C₁[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₁[i,j]);
        Ȧ₂[i,j] = @fastmath Ȧ₂[i,j] .+ 0.5Δt .* (F₂[i,j] .+ M₂[i,j] - λ * C₂[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) - λ * C₂[i,j]*8.4e5 * 1e12/(3fₐ) - 2/η * Ȧ₂[i,j]);
    end
    return nothing
end

function PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,time,fₐ)

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    PQfupdate_2D!(F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,Ȧ₁,Ȧ₂,time,1,fₐ)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    PQAupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))
    if time < 250
        PQvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,time + Δt,1,fₐ)
    else 
        PQvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,time + Δt * 250/time,1,fₐ)
    end

    return nothing
end



function EQCDfupdate_2D!(F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,λ,ηₓ,n)
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


function EQCDvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,η,λ,ηₓ,n)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ₁)
        (i,j) = Tuple(😄)
        Ȧ₁[i,j] = @fastmath Ȧ₁[i,j] .+ 0.5Δt .* (F₁[i,j] .+ M₁[i,j] - λ * C₁[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) + ηₓ^n *η^2 * C₂[i,j]^2/((A₁[i,j]^2 .+ A₂[i,j]^2)^1.5) - 2/η * Ȧ₁[i,j])
        Ȧ₂[i,j] = @fastmath Ȧ₂[i,j] .+ 0.5Δt .* (F₂[i,j] .+ M₂[i,j] - λ * C₂[i,j] * η^2 * (A₁[i,j]^2 .+ A₂[i,j]^2 .- 1) - ηₓ^n *η^2 * C₁[i,j]*C₂[i,j]/((A₁[i,j]^2 .+ A₂[i,j]^2)^1.5) - 2/η * Ȧ₂[i,j])
    end
    return nothing
end


function EQCDupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,time,fₐ,r,s)

    n = 1

    #λ = (fₐ/mass(fₐ))^2;

    λᵪ = [1024 1448 3072 3584 5504];
    #ηₓ,η = ηtime(time,fₐ);
    λ = λᵪ[r]

    ηᵪ = range(2.8,stop = 3.6,length = 5)
    ηₓ = ηᵪ[s]

    if time < ηₓ
        ηₓ = time
    end
    

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    EQCDfupdate_2D!(F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂,Ȧ₁,Ȧ₂,time,λ,ηₓ,n)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    EQCDAupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    #ηₓ,η = ηtime(time+Δt,fₐ);
    if time < 1.8
        EQCDvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂, time + Δt,λ,ηₓ,n)
    else
        EQCDvelupdate_2D!(Ȧ₁,Ȧ₂,Δt,F₁,F₂,M₁,M₂,A₁,A₂,A₁,A₂, time + Δt* (1.8/time)^3.34,λ,ηₓ,n)
    end

    return nothing
end



function LLaplacian_2D!(P,A,Δx)
    B = CircularArray(A);
    Threads.@threads for 😄 ∈ CartesianIndices(P)
        (i,j) = Tuple(😄)
        P[i,j] = @fastmath (-B[i+2,j] + 16B[i+1,j] 
            + 16B[i-1,j] -B[i-2,j] 
            -B[i,j+2] + 16B[i,j+1] 
            + 16B[i,j-1] -B[i,j-2] 
            - 60B[i,j]) /(12Δx^2);
    end
    return nothing
end


function Lfupdate_2D!(F,M,A,Ȧ,η,ηₓ,n)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F)
        (i,j) = Tuple(😄)
        F[i,j] = @fastmath M[i,j] - ηₓ^n * η^2 * sin(A[i,j]) - 2/η * Ȧ[i,j];
    end
    return nothing
end


function LAupdate_2D!(A,Δt,Ȧ,F)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(A)
        (i,j) = Tuple(😄)
        A[i,j] = @fastmath A[i,j] .+ Δt .* (Ȧ[i,j] .+ 0.5Δt .* F[i,j])
    end
    return nothing
end


function Lvelupdate_2D!(Ȧ,Δt,F,M,A,η,ηₓ,n)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ)
        (i,j) = Tuple(😄)
        Ȧ[i,j] = @fastmath Ȧ[i,j] .+ 0.5Δt .* (F[i,j] .+ M[i,j] - ηₓ^n * η^2 * sin(A[i,j]) - 2/η * Ȧ[i,j])
    end
    return nothing
end


function Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ,s)

    #ηₓ,η = Lηtime(time,fₐ);
    n = 1

    ηᵪ = range(2.8,stop = 3.6,length = 5)
    ηₓ = ηᵪ[s]

    if time < ηₓ
        ηₓ = time
    end

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    Lfupdate_2D!(F,M,A,Ȧ,time,ηₓ,n)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    LAupdate_2D!(A,Δt,Ȧ,F)

    LLaplacian_2D!(M,A,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    #ηₓ,η = Lηtime(time+Δt,fₐ)

    Lvelupdate_2D!(Ȧ,Δt,F,M,A,time + Δt,ηₓ,n)

    return nothing
end

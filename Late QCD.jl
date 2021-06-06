#Late QCD
#T << Λ


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


function Lfupdate_2D!(F,M,A,Ȧ,η,ηₓ)
    n = 6.68;
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


function Lvelupdate_2D!(Ȧ,Δt,F,M,A,η,ηₓ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ)
        (i,j) = Tuple(😄)
        Ȧ[i,j] = @fastmath Ȧ[i,j] .+ 0.5Δt .* (F[i,j] .+ M[i,j] - ηₓ^n * η^2 * sin(A[i,j]) - 2/η * Ȧ[i,j])
    end
    return nothing
end



function Lηtime(time,fₐ,T)
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

function Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ)

    ηₓ,η = Lηtime(time,fₐ,50);
    

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    Lfupdate_2D!(F,M,A,Ȧ,η,ηₓ)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    LAupdate_2D!(A,Δt,Ȧ,F)

    LLaplacian_2D!(M,A,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    ηₓ,η = Lηtime(time+Δt,fₐ,50);

    Lvelupdate_2D!(Ȧ,Δt,F,M,A,η,ηₓ)

    return nothing
end

function Lrun_2D!(N,t₀,t,A,Ȧ,Δx,Δt,fₐ,i)

    time = t₀

    M = zeros(N,N);


    F = zeros(N,N);


    LLaplacian_2D!(M,A,Δx)



    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if time % 1 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            PyPlot.imsave("Late/"*string(i)*"/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",A,vmin=-π,vmax = π,cmap = "twilight")
        end
        Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ)
        time = time + Δt

    end

    return time
end

function Lplotting_2D!(N,t₀,t₁,t,A,Ȧ,Δx,Δt,fₐ,i)

    time = t₁

    M = zeros(N,N);

    F = zeros(N,N);

    LLaplacian_2D!(M,A,Δx)

    # moo = zeros(N,N);

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if time % 1 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            PyPlot.imsave("Late/"*string(i)*"/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",A,vmin=-π,vmax = π,cmap = "twilight")
        end
        Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ)
        time = time + Δt

    end

    return nothing
end
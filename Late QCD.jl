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
    n = 6.68
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ)
        (i,j) = Tuple(😄)
        Ȧ[i,j] = @fastmath Ȧ[i,j] .+ 0.5Δt .* (F[i,j] .+ M[i,j] - ηₓ^n * η^2 * sin(A[i,j]) - 2/η * Ȧ[i,j])
    end
    return nothing
end



function Lηtime(time,fₐ)
    t₁ = 1.61e-10 * (fₐ/1e12)
    η = (time/t₁)^0.5
    T = 42.3e3 * (fₐ/1e12)^(-0.5)

    ηₓ = T/103;

    if ηₓ > η
        ηₓ = η
    end

    return ηₓ,η
end

function Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ)

    ηₓ,η = Lηtime(time,fₐ);
    

    #F₁ .= M₁ .- a.^β .* λ .* A₁ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./time
    #F₂ .= M₂ .- a.^β .* λ .* A₂ .*(A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./time

    Lfupdate_2D!(F,M,A,Ȧ,η,ηₓ)

    #A₁ .= A₁ .+ Δt .* (Ȧ₁ .+ 0.5Δt .* F₁)
    #A₂ .= A₂ .+ Δt .* (Ȧ₂ .+ 0.5Δt .* F₂)

    LAupdate_2D!(A,Δt,Ȧ,F)

    LLaplacian_2D!(M,A,Δx)

    #Ȧ₁ .= Ȧ₁ .+ 0.5Δt .* (F₁ .+ M₁ .- a₁.^β .* λ .* A₁ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₁ ./ (time + Δt))
    #Ȧ₂ .= Ȧ₂ .+ 0.5Δt .* (F₂ .+ M₂ .- a₁.^β .* λ .* A₂ .* (A₁.^2 .+ A₂.^2 .- η.^2) .- α .* © .* Ȧ₂ ./ (time + Δt))

    ηₓ,η = Lηtime(time+Δt,fₐ)

    Lvelupdate_2D!(Ȧ,Δt,F,M,A,η,ηₓ)

    return nothing
end

function Lrun_2D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)
    time = t₀

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);

    moo = zeros(N,N);
    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    angle = zeros(N,N);

    #fₛ = Δx
    #freq_s = fₛ/2 * range(-1,stop = 1,length = N)
    #freq = fₛ/2 * range(0,stop = 1, length = N/2)
    k_freq = fftfreq(N)*N
    kx,ky = meshgrid(k_freq,k_freq)

    knrm = sqrt.( kx.^2 + ky.^2)
    knrm = collect(Iterators.flatten(knrm))

    kbins = range(0.5, N/2+1, step = 1)
    kvals = 0.5 * (kbins[2:end] + kbins[1:end-1])


    for lo ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if lo % 10 == 0
            mooing!(moo,A₁,A₂);
            setting!(moo);
            angler!(angle,A₁,A₂);

            PyPlot.imsave("Late/1"*string(i)*"/"*lpad( string(trunc(Int,t₀ + lo/10 - 1)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            
            f_image = FFTW.fft(moo)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotd = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false, ylims = (1e1,1e8))
            Plots.savefig(plotd,"Late/Fourier/"*lpad( string(trunc(Int,t₀ + lo/10 - 1)) ,3,"0")*".png")
            
            f_image = FFTW.fft(angle)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false, ylims = (1e6,1e11))
            Plots.savefig(plotc,"Late/Angle/"*lpad( string(trunc(Int,t₀ + lo/10 - 1)) ,3,"0")*".png")
            

            PyPlot.imsave("Late/"*string(i)*"/"*lpad( string(trunc(Int,t₀ + lo/10 - 1)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end
        PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₀,time,fₐ)
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
    k_freq = fftfreq(N)*N
    kx,ky = meshgrid(k_freq,k_freq)

    knrm = sqrt.( kx.^2 + ky.^2)
    knrm = collect(Iterators.flatten(knrm))

    kbins = range(0.5, N/2+1, step = 1)
    kvals = 0.5 * (kbins[2:end] + kbins[1:end-1])

    for lo ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 10);
        if lo % 10 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            
            f_image = FFTW.fft(A)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))

            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
        
            plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false, ylims = (1e6,1e11))
            Plots.savefig(plotc,"Late/Angle/"*lpad( string(trunc(Int,t₀ + lo/10 - 1)) ,3,"0")*".png")
            
        
            PyPlot.imsave("Late/"*string(i)*"/"*lpad( string(trunc(Int,t₀+lo/10 - 1)) ,3,"0")*".png",A,vmin=-π,vmax = π,cmap = "twilight")
        end
        Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ)
        time = time + Δt

    end

    return nothing
end
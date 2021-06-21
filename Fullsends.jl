# Full sends



function FPQrun_2D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)

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
            PyPlot.imsave("Full/"*string(i)*"/"*lpad( string(trunc(Int,(time-t₀))) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end
        PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,t₀,time,fₐ)
        time = time + Δt

    end

    return time
end



function FErun_2D!(N,t₀,t₁,t,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)

    time = t₁

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    # moo = zeros(N,N);
    angle = zeros(N,N);
    angler!(angle,A₁,A₂);

    for lo ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 10);
        if lo % 10 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
            angler!(angle,A₁,A₂);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            PyPlot.imsave("Full/"*string(i)*"/"*lpad( string(trunc(Int,t₀ + lo/10 - 1)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end
        EQCDupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,time,fₐ)
        time = time + Δt
        
    end

    return time
end


function FLrun_2D!(N,t₀,t₁,t,A,Ȧ,Δx,Δt,fₐ,i)

    time = t₁

    M = zeros(N,N);

    F = zeros(N,N);

    LLaplacian_2D!(M,A,Δx)

    # moo = zeros(N,N);

    for lo ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 10);
        if lo % 10 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            PyPlot.imsave("Full/"*string(i)*"/"*lpad( string(trunc(Int,t₀+lo/10 - 1)) ,3,"0")*".png",A,vmin=-π,vmax = π,cmap = "twilight")
        end
        Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ)
        time = time + Δt

    end

    return nothing
end
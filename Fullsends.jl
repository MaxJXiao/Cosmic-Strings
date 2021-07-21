# Full sends




function FPQrun_2D!(N,t₀,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i)
    #700 images on per 5000 updates
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
    x = 5000


    for lo ∈ 1:round((250- t₀)/Δt,digits = 0)
        time = round(time,digits = 4);


        if lo % x == 0
            #125 pictures
            mooing!(moo,A₁,A₂);
            #println(moo[1,1])
            setting!(moo);
            angler!(angle,A₁,A₂);

            PyPlot.imsave("PQEpoch/1"*string(i)*"/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            
            f_image = FFTW.fft(moo)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotd = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e1,1e8))
            Plots.savefig(plotd,"PQEpoch/Fourier/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png")
            
            f_image = FFTW.fft(angle)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e6,1e11))
            Plots.savefig(plotc,"PQEpoch/Angle/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png")
            

            PyPlot.imsave("PQEpoch/"*string(i)*"/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end


        PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,time,fₐ)
        time = time + Δt

    end
    print(time)
    tim = time
    lo = round((250-t₀)/Δt,digits = 0) + 1

    while tim < 800
        #make a note a 280 where we get PQ transition
        tim = round(tim,digits = 5);


        # if lo % x == 0
        #     mooing!(moo,A₁,A₂);
        #     #println(moo[1,1])
        #     setting!(moo);
        #     angler!(angle,A₁,A₂);

        #     PyPlot.imsave("PQEpoch/1"*string(i)*"/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            
        #     f_image = FFTW.fft(moo)
        #     f_images = (abs.(f_image)).^2
        #     f_images = collect(Iterators.flatten(f_images))
   
        #     Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
        #     Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
        #     plotd = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e1,1e8))
        #     Plots.savefig(plotd,"PQEpoch/Fourier/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png")
            
        #     f_image = FFTW.fft(angle)
        #     f_images = (abs.(f_image)).^2
        #     f_images = collect(Iterators.flatten(f_images))
   
        #     Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
        #     Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
        #     plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e6,1e11))
        #     Plots.savefig(plotc,"PQEpoch/Angle/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png")
            

        #     PyPlot.imsave("PQEpoch/"*string(i)*"/"*lpad( string(trunc(Int,lo/x - 1)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        # end


        PQupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt * 250/tim,tim,fₐ)
        tim = tim + Δt * 250/tim
        lo += 1
    end
        

    return nothing
end



function FErun_2D!(N,t₁,A₁,A₂,Ȧ₁,Ȧ₂,Δx,Δt,fₐ,i,r,s)

    #

    time = t₁

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    moo = zeros(N,N);
    angle = zeros(N,N);
    angler!(angle,A₁,A₂);

    k_freq = fftfreq(N)*N
    kx,ky = meshgrid(k_freq,k_freq)

    knrm = sqrt.( kx.^2 + ky.^2)
    knrm = collect(Iterators.flatten(knrm))

    kbins = range(0.5, N/2+1, step = 1)
    kvals = 0.5 * (kbins[2:end] + kbins[1:end-1])
    
    x = 100

    for lo ∈ 1:round((1.8-t₁)/Δt,digits = 0)
        time = round(time,digits = 3);
        if lo % x == 0
            #18 pictures
            mooing!(moo,A₁,A₂);
            #println(moo[1,1])
            setting!(moo);
            angler!(angle,A₁,A₂);
            PyPlot.imsave("EQCD/1"*string(i)*"/"*lpad( string(trunc(Int,lo/x - 1)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")

            f_image = FFTW.fft(moo)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotd = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e1,1e8))
            Plots.savefig(plotd,"EQCD/Fourier/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png")
            
            f_image = FFTW.fft(angle)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e6,1e11))
            Plots.savefig(plotc,"EQCD/Angle/"*lpad( string(trunc(Int,lo/x - 1)) ,3,"0")*".png")
            

            PyPlot.imsave("EQCD/"*string(i)*"/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")

        end
        EQCDupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt,time,fₐ,r,s)
        time = time + Δt
        
    end

    lo = round((1.8-t₁) /Δt,digits = 0) + 1
    print(time)
    tim = time
    while tim < 7
        
        tim = round(tim,digits = 6);
        if lo % x == 0
            #13 pictures
            mooing!(moo,A₁,A₂);
            #println(moo[1,1])
            setting!(moo);
            angler!(angle,A₁,A₂);
            PyPlot.imsave("EQCD/1"*string(i)*"/"*lpad( string(trunc(Int,lo/x - 1)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")

            f_image = FFTW.fft(moo)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotd = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e1,1e8))
            Plots.savefig(plotd,"EQCD/Fourier/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png")
            
            f_image = FFTW.fft(angle)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))
   
            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e6,1e11))
            Plots.savefig(plotc,"EQCD/Angle/"*lpad( string(trunc(Int,lo/x - 1)) ,3,"0")*".png")
            

            PyPlot.imsave("EQCD/"*string(i)*"/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
            
        end
        EQCDupdate_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,Δx,Δt* (1.8/tim)^3.34,tim,fₐ,r,s)
        tim = tim + Δt* (1.8/tim)^3.34
        lo += 1
    end

    return tim
end






function FLrun_2D!(N,t₁,t,A,Ȧ,Δx,Δt,fₐ,i,s)

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

    x = 100

    for lo ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 3);
        if lo % x == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            
            f_image = FFTW.fft(A)
            f_images = (abs.(f_image)).^2
            f_images = collect(Iterators.flatten(f_images))

            Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
        
            plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)#, ylims = (1e6,1e11))
            Plots.savefig(plotc,"Late/Angle/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png")
            
        
            PyPlot.imsave("Late/"*string(i)*"/"*lpad( string(trunc(Int, lo/x - 1)) ,3,"0")*".png",A,vmin=-π,vmax = π,cmap = "twilight")
        end
        Lupdate_2D!(A,Ȧ,M,F,Δx,Δt,time,fₐ,s)
        time = time + Δt

    end

    return nothing
end
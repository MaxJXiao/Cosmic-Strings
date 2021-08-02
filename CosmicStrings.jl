
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


function fupdate_2D!(F₁,F₂,M₁,M₂,a,©,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,time,β,α,λ)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F₁)
        (i,j) = Tuple(😄)
        F₁[i,j] = @fastmath M₁[i,j] - a^β * λ * C₁[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - α * © * Ȧ₁[i,j] /time;
        F₂[i,j] = @fastmath M₂[i,j] - a^β * λ * C₂[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - α * © * Ȧ₂[i,j] /time;
    end
    return nothing
end

function fupdate_mass2D!(F₁,F₂,M₁,M₂,a,©,C₁,C₂,A₁,A₂,Ȧ₁,Ȧ₂,η,time,β,α,λ,mass,angle,fₐ)
    T = 100/(1e3fₐ);
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(F₁)
        (i,j) = Tuple(😄)
        F₁[i,j] = @fastmath M₁[i,j] - a^β * λ * C₁[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - λ*T^2/3 * C₁[i,j] - α * © * Ȧ₁[i,j] /time
        + mass*η^2 * sin(angle[i,j])* C₂[i,j]/ (A₁[i,j]^2 .+ A₂[i,j]^2);
        F₂[i,j] = @fastmath M₂[i,j] - a^β * λ * C₂[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - λ*T^2/3 * C₂[i,j] - α * © * Ȧ₂[i,j] /time
        - mass * η^2 * sin(angle[i,j]) * C₁[i,j]/ (A₁[i,j]^2 .+ A₂[i,j]^2);
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

function velupdate_mass2D!(Ȧ₁,Ȧ₂,Δt,a,©,F₁,F₂,M₁,M₂,C₁,C₂,A₁,A₂,η,time,β,α,λ,mass,angle,fₐ)
    T = 100/(1e3fₐ);
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(Ȧ₁)
        (i,j) = Tuple(😄)
        Ȧ₁[i,j] = Ȧ₁[i,j] .+ 0.5Δt .* (F₁[i,j] .+ M₁[i,j] - a^β * λ * C₁[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - λ*T^2/3 * C₁[i,j] - α * © * Ȧ₁[i,j] /time
        + mass*η^2 * sin(angle[i,j])* C₂[i,j]/ (A₁[i,j]^2 .+ A₂[i,j]^2) );
        Ȧ₂[i,j] = Ȧ₂[i,j] .+ 0.5Δt .* (F₂[i,j] .+ M₂[i,j] - a^β * λ * C₂[i,j] *(A₁[i,j]^2 .+ A₂[i,j]^2 .- η^2) - λ*T^2/3 * C₂[i,j] - α * © * Ȧ₂[i,j] /time
        - mass * η^2 * sin(angle[i,j]) * C₁[i,j]/ (A₁[i,j]^2 .+ A₂[i,j]^2) );
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

function update_mass2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,ω,η,Δx,Δt,time,mass,angle,fₐ)
    a = 4.36000000000006e-18*time - 6.78288102293483e-23
    a₁ = 4.36000000000006e-18*(time + Δt) - 6.78288102293483e-23

    β = 0
    α = 3
    © = 1
    λ = 2π^2/ω^2

    fupdate_mass2D!(F₁,F₂,M₁,M₂,a,©,A₁,A₂,A₁,A₂,Ȧ₁,Ȧ₂,η,time,β,α,λ,mass,angle,fₐ)


    Aupdate_2D!(A₁,A₂,Δt,Ȧ₁,Ȧ₂,F₁,F₂)

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    angler!(angle,A₁,A₂)

    velupdate_mass2D!(Ȧ₁, Ȧ₂, Δt, a₁, ©, F₁, F₂, M₁, M₂, A₁, A₂, A₁, A₂, η, time + Δt, β, α, λ,mass,angle,fₐ)

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

function thresholding!(moo,t)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(moo)
        (i,j) = Tuple(😄)
        moo[i,j] = moo[i,j] > t
    end
    return moo
end

function angler!(angle,A₁,A₂)
    @inbounds Threads.@threads for 😄 ∈ CartesianIndices(angle)
        (i,j) = Tuple(😄)
        angle[i,j] = @fastmath atan(A₂[i,j],A₁[i,j]);
    end
    return nothing
end



function run_2D!(N,t₀,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,i)

    time = t₀

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);

    moo = zeros(N,N);
    angle = zeros(N,N);

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    # k_freq = fftfreq(N)*N
    # kx,ky = meshgrid(k_freq,k_freq)

    # knrm = sqrt.( kx.^2 + ky.^2)
    # knrm = collect(Iterators.flatten(knrm))

    # kbins = range(0.5, N/2+1, step = 1)
    # kvals = 0.5 * (kbins[2:end] + kbins[1:end-1])

    

    THRESHOLD = 190
    thr = 1
    B = zeros(0)

    font                   = cv2.FONT_HERSHEY_SIMPLEX
    topLeftCornerOfText = (30,50)
    fontScale              = 1
    fontColor              = (0,0,0)
    lineType               = 2

    bottomLeftCornerOfText = (30,500)



    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if time % 1 == 0
            mooing!(moo,A₁,A₂);
            setting!(moo);
            angler!(angle,A₁,A₂);

            
            strng = Cores_2D!(N,angle,thr)

            # moo .= moo .< 0.45
            
            # mood = find_contours(moo)

            # if length(mood) == 1
            #     append!(x,length(mood[1]))
            # else
            #     append!(x,0)
            # end

            # f_image = FFTW.fft(moo)
            # f_images = (abs.(f_image)).^2
            # f_images = collect(Iterators.flatten(f_images))
   
            # Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            # Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            # plotd = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false, ylims = (1e2,1e8))
            # Plots.savefig(plotd,"plottting_angle/Fourier/"*lpad( string(trunc(Int,(time-t₀))) ,3,"0")*".png")
            

            # f_image = FFTW.fft(angle)
            # f_images = (abs.(f_image)).^2
            # f_images = collect(Iterators.flatten(f_images))
   
            # Abins,_,_ = stats.binned_statistic(knrm,f_images,statistic = "mean",bins = kbins)
            # Abins = π* Abins.* (kbins[2:end].^2 - kbins[1:end-1].^2)
         
            # plotc = Plots.plot(kvals,Abins,xaxis= :log,yaxis =:log,legend = false)
            # Plots.savefig(plotc,"plottting_angle/Angle/"*lpad( string(trunc(Int,(time-t₀))) ,3,"0")*".png")
            
            #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));

            PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
            PyPlot.imsave("plottting_angle/"*string(i)*"/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")


            im = cv2.imread("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png")

            imgray = cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)

            ret,thresh = cv2.threshold(imgray,THRESHOLD, 255,0)
            contours,hierachy = cv2.findContours(cv2.bitwise_not(thresh),cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)

            

            append!(B,length(contours))

            

            C = cv2.drawContours(cv2.UMat(im),contours,-1,(1,255,1),3)
     
            C = cv2.putText(C,string(length(contours)), 
            topLeftCornerOfText, 
            font, 
            fontScale,
            fontColor,
            lineType)

            C = cv2.putText(C,string(strng),bottomLeftCornerOfText, 
            font, 
            fontScale,
            fontColor,
            lineType)

            cv2.imwrite("plottting_m/1/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",C)


        end
        update_2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,ω,η,Δx,Δt,time)
        time = time + Δt

    end

    return time
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

function plotting_2D!(N,t₀,t₁,t,A₁,A₂,Ȧ₁,Ȧ₂,ω,η,Δx,Δt,fₐ,i)

    time = t₁

    M₁ = zeros(N,N);
    M₂ = zeros(N,N);

    F₁ = zeros(N,N);
    F₂ = zeros(N,N);

    Laplacian_2D!(M₁,M₂,A₁,A₂,Δx)

    # moo = zeros(N,N);
    mas = mass(fₐ,100);
    #mas = 1e-6;
    angle = zeros(N,N);
    angler!(angle,A₁,A₂);

    for _ ∈ 1:round(t/Δt,digits = 0)
        time = round(time,digits = 1);
        if time % 1 == 0
        #     mooing!(moo,A₁,A₂);
        #     setting!(moo);
        #     #save("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png", colorview(Gray,moo));
        #     PyPlot.imsave("plottting_m/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",moo,vmin=0,vmax = 1,cmap = "gray")
             PyPlot.imsave("plottting_angle/"*string(i)*"/"*lpad( string(trunc(Int,time-t₀)) ,3,"0")*".png",angle,vmin=-π,vmax = π,cmap = "twilight")
        end
        update_mass2D!(A₁,A₂,Ȧ₁,Ȧ₂,M₁,M₂,F₁,F₂,ω,η,Δx,Δt,time,mas,angle,fₐ)
        time = time + Δt

    end

    return nothing
end


function Cores_3D!(N,angle,thr)
    s = []
    count = 0
    accept = 0.5 - 0.5*thr/100
    for 😄 ∈ 1:(N-1)
        for 🥪 ∈ 1:(N-1)
            for 🎅 ∈ 1:(N-1)
                anorm₁ = (angle[😄,🥪,🎅] + π)/(2π)
                anorm₂ = (angle[😄+1,🥪,🎅] + π)/(2π)
                anorm₃ = (angle[😄+1,🥪+1,🎅] + π)/(2π)
                anorm₄ = (angle[😄,🥪+1,🎅] + π)/(2π)

                bnorm₁ = (angle[😄,🥪,🎅] + π)/(2π)
                bnorm₂ = (angle[😄+1,🥪,🎅] + π)/(2π)
                bnorm₃ = (angle[😄+1,🥪,🎅+1] + π)/(2π)
                bnorm₄ = (angle[😄,🥪,🎅+1] + π)/(2π)

                cnorm₁ = (angle[😄,🥪,🎅] + π)/(2π)
                cnorm₂ = (angle[😄,🥪+1,🎅] + π)/(2π)
                cnorm₃ = (angle[😄,🥪+1,🎅+1] + π)/(2π)
                cnorm₄ = (angle[😄,🥪,🎅+1] + π)/(2π)

                aθ₁ = min(abs(anorm₂ - anorm₁),1-abs(anorm₂ - anorm₁))
                aθ₂ = min(abs(anorm₃ - anorm₂),1-abs(anorm₃ - anorm₂))
                aθ₃ = min(abs(anorm₄ - anorm₃),1-abs(anorm₄ - anorm₃))

                bθ₁ = min(abs(bnorm₂ - bnorm₁),1-abs(bnorm₂ - bnorm₁))
                bθ₂ = min(abs(bnorm₃ - bnorm₂),1-abs(bnorm₃ - bnorm₂))
                bθ₃ = min(abs(bnorm₄ - bnorm₃),1-abs(bnorm₄ - bnorm₃))

                cθ₁ = min(abs(cnorm₂ - cnorm₁),1-abs(cnorm₂ - cnorm₁))
                cθ₂ = min(abs(cnorm₃ - cnorm₂),1-abs(cnorm₃ - cnorm₂))
                cθ₃ = min(abs(cnorm₄ - cnorm₃),1-abs(cnorm₄ - cnorm₃))

                aθ = aθ₁ + aθ₂ + aθ₃
                bθ = bθ₁ + bθ₂ + bθ₃
                cθ = cθ₁ + cθ₂ + cθ₃
                θₛ = aθ + bθ + cθ

                if θₛ > accept 
                    append!(s,[[😄,🥪,🎅]])
                end
            end
        end
    end

    if length(s) > 0
        for 🇸🇦 ∈ 1:(length(s)-1)
  
            diffᵣ = s[🇸🇦 + 1][1] - s[🇸🇦][1]
            diffₛ = s[🇸🇦 + 1][2] - s[🇸🇦][2]
            diffₜ = s[🇸🇦 + 1][3] - s[🇸🇦][3]

            if diffᵣ == 0 && diffₛ == 1 && diffₜ == 0
                count += 1
            end
            if diffᵣ == 1 && diffₛ == 0 && diffₜ == 0
                count += 1
            end
            if diffᵣ == 0 && diffₛ == 0 && diffₜ == 1
                count += 1
            end
        end
    end

    num = length(s) - count
    num = 1.5num

    return num
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


  
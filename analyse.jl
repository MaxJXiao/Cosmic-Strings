using Base: has_tight_type

saxion = load("Saving/PQ/PQStatssecond9.jld")["saxion"]
axion = load("Saving/PQ/PQStatssecond9.jld")["axion"]
tracker = load("Saving/PQ/PQStringssecond9.jld")["time"]
axenergy = load("Saving/PQ/PQStatssecond9.jld")["axenergy"]
cores = load("Saving/PQ/PQStringssecond9.jld")["number"]

n = 9

N = 2^n

k_freq = fftfreq(N)*N
kx,ky = meshgrid(k_freq,k_freq)

knrm = sqrt.( kx.^2 + ky.^2)
knrm = collect(Iterators.flatten(knrm))

kbins = range(0.5, N/2+1, step = 1)
kvals = 0.5 * (kbins[2:end] + kbins[1:end-1])

lowₛ = minimum(minimum(saxion))
highₛ = maximum(maximum(saxion))
lowₐ = minimum(minimum(axion))
highₐ = maximum(maximum(axion))

for i ∈ 1:length(tracker)

plotc = Plots.plot(kvals,saxion[i],yaxis =:log,legend = false, ylims = (lowₛ,highₛ),title = "Saxion η = "*string(tracker[i]))
Plots.savefig(plotc,"Saving/PQ/Saxion/"*lpad( string(i-1) ,3,"0")*".png")


plotd = Plots.plot(kvals,axion[i],yaxis =:log,legend = false, ylims = (lowₐ,highₐ),title = "Axion (a/fₐ) η = "*string(tracker[i]))
Plots.savefig(plotd,"Saving/PQ/Axion/"*lpad( string(i-1) ,3,"0")*".png")

end



plote = Plots.plot(tracker,axenergy/1e30,yaxis = :log,title = "Axion Energy Density ρₐ = ⟨ ȧ^2 / fₐ^2 ⟩",legend = false)
Plots.savefig(plote,"Saving/PQ/Energy_Density.png")

plotf = Plots.plot(tracker,cores,title = "String Count",legend = false,yaxis = :log)
Plots.savefig(plotf,"Saving/PQ/Strings.png")









saxion = load("Saving/QCD/QCDStatsfourth911.jld")["saxion"]
axion = load("Saving/QCD/QCDStatsfourth911.jld")["axion"]
tracker = load("Saving/QCD/QCDStringsfourth911.jld")["time"]
axenergy = load("Saving/QCD/QCDStatsfourth911.jld")["axenergy"]
cores = load("Saving/QCD/QCDStringsfourth911.jld")["number"]

lowₛ = minimum(minimum(saxion))
highₛ = maximum(maximum(saxion))
lowₐ = minimum(minimum(axion))
highₐ = maximum(maximum(axion))

for i ∈ 1:length(saxion)

    plotc = Plots.plot(kvals,saxion[i],yaxis =:log,legend = false, ylims = (lowₛ,highₛ),title = "λ = 1024 ,ηᵪ = 2.8 ,Saxion η = "*string(tracker[i]))
    Plots.savefig(plotc,"Saving/QCD/Saxion/911/"*lpad( string(i-1) ,3,"0")*".png")
    
end

for i ∈ 1:length(axion)

    plotd = Plots.plot(kvals,axion[i],yaxis =:log,legend = false, ylims = (lowₐ,highₐ),title = "λ = 1024 ,ηᵪ = 2.8 ,Axion (a/fₐ) η = "*string(tracker[i]))
    Plots.savefig(plotd,"Saving/QCD/Axion/911/"*lpad( string(i-1) ,3,"0")*".png")

end



plote = Plots.plot(tracker,axenergy/1e30,yaxis = :log,title = "λ = 1024 ,ηᵪ = 2.8 ,Axion Energy Density ρₐ = ⟨ ȧ^2 / fₐ^2 ⟩",legend = false)
Plots.savefig(plote,"Saving/QCD/Energy_Density911.png")

plotf = Plots.plot(tracker,cores,title = "λ = 1024 ,ηᵪ = 2.8 ,String Count",legend = false)
Plots.savefig(plotf,"Saving/QCD/Strings911.png")







saxion = load("Saving/QCD/QCDStatsfourth915.jld")["saxion"]
axion = load("Saving/QCD/QCDStatsfourth915.jld")["axion"]
tracker = load("Saving/QCD/PQStringsfourth915.jld")["time"]
axenergy = load("Saving/QCD/QCDStatsfourth915.jld")["axenergy"]
cores = load("Saving/QCD/QCDStringsfourth915.jld")["number"]

lowₛ = minimum(minimum(saxion))
highₛ = maximum(maximum(saxion))
lowₐ = minimum(minimum(axion))
highₐ = maximum(maximum(axion))

for i ∈ 1:length(saxion)

    plotc = Plots.plot(kvals,saxion[i],yaxis =:log,legend = false, ylims = (lowₛ,highₛ),title = "λ = 1024 ,ηᵪ = 3.6 ,Saxion η = "*string(tracker[i]))
    Plots.savefig(plotc,"Saving/QCD/Saxion/915/"*lpad( string(i-1) ,3,"0")*".png")
    
end

for i ∈ 1:length(axion)

    plotd = Plots.plot(kvals,axion[i],yaxis =:log,legend = false, ylims = (lowₐ,highₐ),title = "λ = 1024 ,ηᵪ = 3.6 ,Axion (a/fₐ) η = "*string(tracker[i]))
    Plots.savefig(plotd,"Saving/QCD/Axion/915/"*lpad( string(i-1) ,3,"0")*".png")

end



plote = Plots.plot(tracker,axenergy/1e30,yaxis = :log,title = "λ = 1024 ,ηᵪ = 3.6 ,Axion Energy Density ρₐ = ⟨ ȧ^2 / fₐ^2 ⟩",legend = false)
Plots.savefig(plote,"Saving/QCD/Energy_Density915.png")

plotf = Plots.plot(tracker,cores,title = "λ = 1024 ,ηᵪ = 3.6 ,String Count",legend = false)
Plots.savefig(plotf,"Saving/QCD/Strings915.png")








saxion = load("Saving/QCD/QCDStatsfourth951.jld")["saxion"]
axion = load("Saving/QCD/QCDStatsfourth951.jld")["axion"]
tracker = load("Saving/QCD/QCDStringsfourth951.jld")["time"]
axenergy = load("Saving/QCD/QCDStatsfourth951.jld")["axenergy"]
cores = load("Saving/QCD/QCDStringsfourth951.jld")["number"]

lowₛ = minimum(minimum(saxion))
highₛ = maximum(maximum(saxion))
lowₐ = minimum(minimum(axion))
highₐ = maximum(maximum(axion))

for i ∈ 1:length(saxion)

    plotc = Plots.plot(kvals,saxion[i],yaxis =:log,legend = false, ylims = (lowₛ,highₛ),title = "λ = 5504 ,ηᵪ = 2.8 ,Saxion η = "*string(tracker[i]))
    Plots.savefig(plotc,"Saving/QCD/Saxion/951/"*lpad( string(i-1) ,3,"0")*".png")
    
end

for i ∈ 1:length(axion)

    plotd = Plots.plot(kvals,axion[i],yaxis =:log,legend = false, ylims = (lowₐ,highₐ),title = "λ = 5504 ,ηᵪ = 2.8 ,Axion (a/fₐ) η = "*string(tracker[i]))
    Plots.savefig(plotd,"Saving/QCD/Axion/951/"*lpad( string(i-1) ,3,"0")*".png")

end



plote = Plots.plot(tracker,axenergy/1e30,yaxis = :log,title = "λ = 5504 ,ηᵪ = 2.8 ,Axion Energy Density ρₐ = ⟨ ȧ^2 / fₐ^2 ⟩",legend = false)
Plots.savefig(plote,"Saving/QCD/Energy_Density951.png")

plotf = Plots.plot(tracker,cores,title = "λ = 5504 ,ηᵪ = 2.8 ,tring Count",legend = false)
Plots.savefig(plotf,"Saving/QCD/Strings951.png")

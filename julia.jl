x=3

print(x)

print(x^2.5)

mdr = 7

δ = mdr

typeof(δ)

😄 = 4

print("$(😄) a")

nt = (x =5, y = "str", z = 5/3)

nt[1]

nt.x

nt[:x]

keys = Dict("key1" => "😄🤣", "key2" => "ghostbusters")
keys["key1"]

keys


rand(4,3)

fib = [1,2,3,4,5,"hi"]

push!(fib,3)

pop!(fib)

fib

#range(start,end; length = )
#range(start, end;step = )

range(0, 10, step = 1)
range(0,10 , length = 5)



for (key,val) in keys
    print("$key $val abc")

end

N = 1:100
for n in N
    isodd(n) && continue
    println(n)
    n > 10 && break
end

g(x) = x^2

g(2)

f(x, y =5,z =2) = x*y*z

f(5)
f(5,3)
f(5, 5,3)

A = rand(3,3)

g(A)

function add3!(x)
    x[1] += 3
    x
end

x = [5,5,5]
add3!(x)
sort(x)

print(x)

f.([3,4,5])

e = 10.0 .^ (-3:3)

#using BenchmarkTools
#@btime

#put all performance critical parts of your code inside a function to avoid global scope
#ensure that your functions are type stable


#try to create as little new large mutable entries in your functions as possible
#try to make the containers first

using Pkg
Pkg.add("OrdinaryDiffEq")
Pkg.add("StaticArrays")
Pkg.add("PyPlot")


#using PyPlot as plt 
using OrdinaryDiffEq, StaticArrays

# function lorenz_rule(u,p,t)
#     σ,ρ,β = p
#     x,y,z = u 
#     dx = σ*(y-x) 
#     dy = x*(ρ-z) - y
#     dz = x*y - β*z
#     return SVector(dx,dy,dz)
# end

u_θ = SVector(10.0,10.0,10.0)
p_θ = [10,28,8/3]

prob = ODEProblem(lorenz_rule,u_θ,(0.0,100.0),p_θ)

alg = Vern9()

sol = solve(prob; alg = alg)

(sol.t[3],sol.u[3])

using PyPlot

Pkg.add("Plots")
Pkg.add("CSV")

using Plots




plot3D(sol[1,:],sol[2,:],sol[3,:])

t = 0:0.01:10
X,Y,Z = zero(t) , zero(t) , zero(t)
for (i,τ) in enumerate(t)
    X[i],Y[i],Z[i] = sol(τ)
end
plot3D(X,Y,Z)
#display(gcf())
#ion()

#PyPlot.close("all")

#typeof

Pkg.add("Measurements")

using Measurements

uunc = SVector(10.0 ± 0.1, 10.0 ± 0.1, 10.0 ±0.1)
probunc = ODEProblem(lorenz_rule,uunc,(0.0,10.0),p_θ)

solunc = solve(probunc ; alg = alg)

t = 0:0.05:5
Xunc = [solunc(τ)[1] for τ in t]

Xv = Measurements.value.(Xunc) #get values
Xe = Measurements.uncertainty.(Xunc)#get errors

plot(t,Xv)
fill_between(t,Xv .-Xe,Xv.+ Xe, color = "C3",alpha = 0.25)

condition(u,t,integrator) = u[2]

function affect!(integrator)
    u = integrator.u
    integrator.u = (u[1],u[2],0.9u[3])
end

cb = ContinuousCallback(condition,affect!)

sol = solve(prob;alg = Vern9(),callback = cb)

p = plot3D(sol[1,:],sol[2,:],sol[3,:])

Pkg.add("DynamicalSystems")

using DynamicalSystems

lorenz = ContinuousDynamicalSystem(lorenz_rule,u_θ,p_θ)

lyapunovs(lorenz,10000)

using Plots

# tagsave("test.bson", copy(pars);gitpath = pwd())

#safesave("test.bson",pars)



using Catch22
using Clustering
using CSV
using Plots
using StatsBase
using UrlDownload

pyplot()
X = urldownload("https://ndownloader.figshare.com/files/24950795", true, format=:CSV, header=false, delim=",", type=Float64, silencewarnings=true) |> Array; # Ignore the warnings

nomissing = F -> F[.!ismissing.(Array(F))]
X = [nomissing(collect(x)) for x in X];

F = hcat([catch22(Float64.(x)) for x in X]...); # May take a minute
Df = 1.0.-abs.(StatsBase.corspearman(F'))
idxs = Clustering.hclust(Df; linkage=:average, branchorder=:optimal).order

p2 = plot(Df[idxs, idxs], seriestype = :heatmap, aspect_ratio=:equal, xaxis=nothing)

plot!(yticks=(1:size(Df, 1), replace.(string.(Catch22.featureNames[idxs]), '_'=>"\\_")), size=(800, 400), xlims=[0.5, size(Df, 1)+0.5], ylims=[0.5, size(Df, 1)+0.5], box=:on, colorbar_title="1-|ρ|", clims=(0.0, 1.0))

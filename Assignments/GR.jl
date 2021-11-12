# using Pkg
# Pkg.add("DifferentialEquations")
# Pkg.add("Plots")
using DifferentialEquations
using Plots


function derivative(y,p,t)
    t,r,phi,ut,ur,uphi = y 
    m = 1
    dydtau = [ut,ur,uphi, -2m/r^2 *(1 - 2m/r)^(-1) * ut * ur , (1 - 2m/r)*( -m/r^2 * ut^2 + m/r^2* (1 - 2m/r)^(-2)* ur^2 + r* uphi^2), -2/r * uphi * ur ]
    return dydtau
end

#for massive object
R,C,omega,m = 10, (10sqrt(7))^(-1), -1, 1
tau = (0.0, 1000.0)
y = [0,R,0,sqrt( ( ( R^2* C^2 - omega)*R) /(R - 2m)), 0, C ]

prob = ODEProblem(derivative,y,tau)
sol = solve(prob,Tsit5(),reltol= 1e-10,abstol = 1e-10,saveat = (0:0.1:1000))


t = []
r = []
phi = []
ut = []
ur = []
uphi = []

for i in 1:length(sol)
append!(t,sol[i][1])
append!(r,sol[i][2])
append!(phi,sol[i][3])
append!(ut,sol[i][4])
append!(ur,sol[i][5])
append!(uphi,sol[i][6])
end

taur = range(0,1000,length = length(t))


plot(taur,t,legend = true,label = 
"C = (10sqrt(7))^(-1)",title ="Proper time vs time"
,xlabel = "tau",ylabel = 't')

plot(taur,r,legend = true,label = "C = (10sqrt(7))^(-1)"
,title ="Proper time vs radius "
,xlabel = "tau",ylabel = 'r')

plot(taur,phi,legend = true,label = "C = (10sqrt(7))^(-1)",
title ="Proper time vs Azimuthal Angle "
,xlabel = "tau",ylabel = "phi")

x = r .* cos.(phi)
y = r .* sin.(phi)


plot(x,y,legend = true,label = "C = (10sqrt(7))^(-1)",title = "Circular Orbit",xlabel = 'x',ylabel = 'y')


#for massive object
R,C,omega,m = 10, 1.1(10sqrt(7))^(-1), -1, 1
tau = (0.0, 1000.0)
y = [0,R,0,sqrt( ( ( R^2* C^2 - omega)*R) /(R - 2m)), 0, C ]

prob = ODEProblem(derivative,y,tau)
sol = solve(prob,Tsit5(),reltol= 1e-10,abstol = 1e-10,saveat = (0:0.1:1000))


t = []
r = []
phi = []
ut = []
ur = []
uphi = []

for i in 1:length(sol)
    append!(t,sol[i][1])
    append!(r,sol[i][2])
    append!(phi,sol[i][3])
    append!(ut,sol[i][4])
    append!(ur,sol[i][5])
    append!(uphi,sol[i][6])
end

taur = range(0,1000,length = length(t))


x = r .* cos.(phi)
y = r .* sin.(phi)


plot(x,y,legend = true,label = "C = 1.1(10sqrt(7))^(-1)",title = "Path of massive object",xlabel = 'x',ylabel = 'y')


m = 1
nu = -(1 .-2m./r).*ut.^2 .+ (1 .- 2m./r).^(-1) .* ur.^2 .+ r.^2 .*uphi.^2
plot(taur ,nu, ylims =(-1-1e-8,-1 + 1e-8),title = "Normalisation of 4 velocity",xlabel = "tau", ylabel = "uu",label = "uu")



    


function derivative_p(y,p,t)
    t,r,phi,ut,ur,uphi = y 
    m = 1
    dydtau = [ut,ur,uphi, -2m/r^2 *(1 - 2m/r)^(-1) * ut * ur , (1 - 2m/r)*( -m/r^2 * ut^2 + m/r^2* (1 - 2m/r)^(-2)* ur^2 + r* uphi^2), -2/r * uphi * ur ]
    return dydtau
end


#for a photon
R,C,omega,m = 10, 0.5, 0, 1
lambda = (0.0, 10.0)
y = [0,R,0,sqrt( ( ( R^2* C^2 - omega)*R) /(R - 2m)), 0, C ]


prob = ODEProblem(derivative_p,y,lambda)
sol = solve(prob,Tsit5(),reltol= 1e-10,abstol = 1e-10,saveat = (0:0.01:10))


xs = []
ys = []


t = []
r = []
phi = []
ut = []
ur = []
uphi = []

for i in 1:length(sol)
    append!(t,sol[i][1])
    append!(r,sol[i][2])
    append!(phi,sol[i][3])
    append!(ut,sol[i][4])
    append!(ur,sol[i][5])
    append!(uphi,sol[i][6])
end

x = r .* cos.(phi)
y = r .* sin.(phi)

append!(xs,[x])
append!(ys,[y])


plot(x,y,legend = true,label = "C = 0.5",title = "Path taken by light",xlabel = 'x',ylabel = 'y')



#for a photon
R,C,omega,m = 10, 1, 0, 1
lambda = (0.0, 10.0)
y = [0,R,0,sqrt( ( ( R^2* C^2 - omega)*R) /(R - 2m)), 0, C ]


prob = ODEProblem(derivative_p,y,lambda)
sol = solve(prob,Tsit5(),reltol= 1e-10,abstol = 1e-10,saveat = (0:0.01:10))


t = []
r = []
phi = []
ut = []
ur = []
uphi = []

for i in 1:length(sol)
    append!(t,sol[i][1])
    append!(r,sol[i][2])
    append!(phi,sol[i][3])
    append!(ut,sol[i][4])
    append!(ur,sol[i][5])
    append!(uphi,sol[i][6])
end

x = r .* cos.(phi)
y = r .* sin.(phi)

append!(xs,[x])
append!(ys,[y])


plot(x,y,legend = true,label = "C = 1",title = "Path taken by Light",xlabel = 'x',ylabel = 'y')


#for a photon
R,C,omega,m = 10, 2, 0, 1
lambda = (0.0, 10.0)
y = [0,R,0,sqrt( ( ( R^2* C^2 - omega)*R) /(R - 2m)), 0, C ]


prob = ODEProblem(derivative_p,y,lambda)
sol = solve(prob,Tsit5(),reltol= 1e-10,abstol = 1e-10,saveat = (0:0.01:10))



t = []
r = []
phi = []
ut = []
ur = []
uphi = []

for i in 1:length(sol)
    append!(t,sol[i][1])
    append!(r,sol[i][2])
    append!(phi,sol[i][3])
    append!(ut,sol[i][4])
    append!(ur,sol[i][5])
    append!(uphi,sol[i][6])
end

x = r .* cos.(phi)
y = r .* sin.(phi)

append!(xs,[x])
append!(ys,[y])


plot(x,y,legend = true,label = "C = 2",title = "Path taken by light",xlabel = 'x',ylabel = 'y')


plot()
plot!(xs[1],ys[1],legend = true,label = "C = 0.5",title = "Path taken by light",xlabel = 'x',ylabel = 'y')
plot!(xs[2],ys[2],label = "C = 1")
plot!(xs[3],ys[3],label = "C = 2")



m = 1
taur = range(0,10,length = length(t))
nu = -(1 .-2m./r).*ut.^2 .+ (1 .- 2m./r).^(-1) .* ur.^2 .+ r.^2 .*uphi.^2
plot(taur ,nu,title = "Normalisation of 4 velocity",xlabel = "lambda", ylabel = "uu",label = "uu")
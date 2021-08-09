using Plots
using LinearAlgebra
using StaticArrays

function Wpoly6(xi,xj,h)
    sigma = 315/(64*pi*h^9)
    q = norm(Differenz(xi,xj))
    if 0<=q<=h
        W = sigma*(h^2-q^2)^3
    else
        W=0.0
    end
    return W
end

function dWspiky(xi,xj,h)
    r = Differenz(xi,xj)
    q = norm(r)
    sigma = -45/(q*pi*h^6)
    if 0<=q<=h
        W= r*(sigma * (h-q)^2 )
    else
        W=SVector(0.0, 0.0)
    end
    return W
end

function d2Wviscosity(xi,xj,h)
    sigma = 45/(pi*h^6)
    q = norm(Differenz(xi,xj))
    if 0<=q<=h
        W = sigma *(h-q)
    else
        W=0.0
    end
    return W
end

function Gitter(partikelrand_links, partikelrand_rechts, partikelrand_unten, partikelrand_oben,
                gittergroeße)
    ax = partikelrand_links+((abs(partikelrand_rechts)+(abs(partikelrand_links)))/(2*gittergroeße))
    bx = partikelrand_rechts-((abs(partikelrand_rechts)+(abs(partikelrand_links)))/(2*gittergroeße))
    abx = bx-ax
    xsmall=collect(ax:abx/(gittergroeße-1):bx)
    ay = partikelrand_unten+((abs(partikelrand_oben)+(abs(partikelrand_unten)))/(2*gittergroeße))
    by = partikelrand_oben - ((abs(partikelrand_oben)+(abs(partikelrand_unten)))/(2*gittergroeße))
    aby = by-ay
    ysmall=collect(ay:aby/(gittergroeße-1):by)
    x = zeros(length(xsmall)*length(ysmall))
    y = zeros(length(xsmall)*length(ysmall))
    for i in 1:length(xsmall)
        for j in 1:length(ysmall)
            x[(i-1)*length(xsmall)+j]=xsmall[i]
            y[(i-1)*length(xsmall)+j]=ysmall[j]
        end
    end
    return x,y,[x y],xsmall,ysmall
end

function Wandpartikel(wandwert_links, wandwert_unten, wandwert_rechts)
    wy_unten = zeros(length(ysmall))
    wy_unten = wy_unten .+wandwert_unten
    wand_unten = [xsmall wy_unten]
    wx_links = zeros(length(xsmall))
    wx_links = wx_links.+wandwert_links
    wand_links = [wx_links ysmall]
    wx_rechts = zeros(length(xsmall))
    wx_rechts = wx_rechts.+wandwert_rechts
    wand_rechts = [wx_rechts ysmall]
    return wand_links,wand_unten,wand_rechts
end

function Wandpartikel2()
    wand_x = collect(wandwert_links:0.05:wandwert_rechts)
    wand_y = collect(wandwert_unten:0.05:wandwert_oben)
    hilfswand_unten = zeros(length(wand_x)) .+ wandwert_unten
    hilfswand_links = zeros(length(wand_y)) .+ wandwert_links
    hilfswand_rechts = zeros(length(wand_y)) .+ wandwert_rechts
    wand_unten = [wand_x hilfswand_unten]
    wand_links = [hilfswand_links wand_y]
    wand_rechts = [hilfswand_rechts wand_y]
    return wand_links,wand_unten,wand_rechts
end

function Differenz(xi,xj)
    xi_x = xi[1,1]
    xi_y = xi[2,1]
    xj_x = xj[1,1]
    xj_y = xj[2,1]
    minimum_x = xi_x - xj_x
    minimum_y = xi_y - xj_y
    return SVector(minimum_x, minimum_y)
end

function Dichte(pos, m, h)
    rho = zeros(size(pos,1))
    for i in 1:length(rho)
        xi = SVector(pos[i, 1], pos[i, 2])
        rho[i] = rho[i] + m*Wpoly6(xi,xi,h)
        for j in (i+1):length(rho)
            xj = SVector(pos[j, 1], pos[j, 2])
            rho_ij = m*Wpoly6(xi,xj,h)
            rho[i] = rho[i] + rho_ij
            rho[j] = rho[j] + rho_ij
        end
    end
    return rho
end

function Druck(rho,rest_density,k)
    "P = k*(rho.*rho)"
    P = k*(rho.-rest_density)
    return P
end

function F_pressure(rho,P,pos,m,h)
    partikelanzahl = length(rho)
    f_pres = zeros(partikelanzahl,2)
    for i in 1:partikelanzahl
        sum = SVector(0.0, 0.0)
        xi = SVector(pos[i, 1], pos[i, 2])
        for j in 1:partikelanzahl
            if i!=j
                xj = SVector(pos[j, 1], pos[j, 2])
                sum = sum + (m./rho[j])*(((P[i]+P[j])./2)*dWspiky(xi,xj,h))
                "sum = sum - rho[i]*m*((P[i]/(rho[i].*rho[i]))+(P[j]/(rho[j].*rho[j]))) .* dWspiky(xi,xj)"
            end
        end
        f_pres[i,1] = sum[1]
        f_pres[i,2] = sum[2]
    end
    return f_pres
end

function F_viscosity(rho,pos,v,m,mu,h)
     partikelanzahl = length(rho)
     f_vis = zeros(partikelanzahl,2)
     for i in 1:partikelanzahl
         sum = SVector(0.0, 0.0)
         xi = SVector(pos[i, 1], pos[i, 2])
         vi = SVector(v[i, 1], v[i, 2])
         for j in 1:partikelanzahl
             if i !=j
                 xj = SVector(pos[j, 1], pos[j, 2])
                 vj = SVector(v[j, 1], v[j, 2])
                 "sum = sum + (m/rho[j])*d2Wviscosity(xi,xj).*v[j,:]"
                 sum = sum + (m/rho[j])*d2Wviscosity(xi,xj,h).*(vj.-vi)
             end
         end
         f_vis[i,1] = mu.*sum[1,1]
         f_vis[i,2] = mu.*sum[2,1]
     end
     return f_vis
end

function F_gravity(partikelanzahl)
    f_grav = zeros(partikelanzahl,2)
    f_grav[:,2] .= -9.81
    "return f_grav.*rest_density"
    return f_grav
end

"Initialisierung ds Gitters"

function main()
partikelrand_links = 0
partikelrand_rechts = 1
partikelrand_unten = 0
partikelrand_oben = 2


gittergroeße = 25


wandwert_rechts = 4
wandwert_links = 0
wandwert_unten = 0
wandwert_oben = 20
x,y,anfangspos,xsmall,ysmall = Gitter(partikelrand_links, partikelrand_rechts, partikelrand_unten,
                                      partikelrand_oben, gittergroeße)
partikelanzahl = length(x)

"Initialisierung der Konstanten"

star_mass = 2
m = star_mass/partikelanzahl
"m = 0.02"
h = 0.0457
mu = 3.5
n = 1
"k = 3"
k = 0.5
"rest_density = 998.29"
rest_density = 2.816
surface_tension = 0.0728
restitution = 0
#plotanzahl = 16
#plotanzeige = 4
#plotindex = 1
#ausgabe_index = 2
motion_damping = 1

"Initialisierung der Variablen"

t_end = 800
dt = 0.01
v = ones(partikelanzahl,2)
v[:,1] = v[:,1].*0
v[:,2] = v[:,2].*0
pos = anfangspos
#ausgabe_x = zeros(partikelanzahl,plotanzahl)
#ausgabe_y = zeros(partikelanzahl,plotanzahl)
#ausgabe_x[:,1] = pos[:,1]
#ausgabe_y[:,1] = pos[:,2]

"Startbedingungen"

rho = Dichte(pos, m, h)
P = Druck(rho,rest_density,k)
f_pres = F_pressure(rho,P,pos,m,h)
f_vis = F_viscosity(rho,pos,v,m,mu,h)
f_grav = F_gravity(partikelanzahl)
#f_gesamt = f_pres + f_vis + f_grav
f_gesamt = -f_pres + f_grav
a = f_gesamt ./ rho
#@show sum(a) #Test der Werte, damit sie gleich bleiben
time_plot = 1

"Leap-frog-time-integration"
v_minus = v - (dt/2)*a

for i in 1:t_end

    v_plus = v_minus + a*dt
    pos = pos + v_plus*dt
    v_minus = v_plus

    "Statische Ränder"
    for j in 1:partikelanzahl
        if pos[j,1] .> wandwert_rechts
            pos[j,1] = 2*wandwert_rechts - pos[j,1]
            v[j,1] = v[j,1] -(1+restitution)*v[j,1]
        end
        if pos[j,1] .< wandwert_links
            pos[j,1] = 2*wandwert_links - pos[j,1]
            v[j,1] = v[j,1] -(1+restitution)*v[j,1]
        end
        if pos[j,2] .< wandwert_unten
            pos[j,2] = 2*wandwert_unten - pos[j,2]
            v[j,2] = v[j,2] -(1+restitution)*v[j,2]
        end
        if pos[j,2] .> wandwert_oben
            pos[j,2] = 2*wandwert_oben - pos[j,2]
            v[j,2] = v[j,2] -(1+restitution)*v[j,2]
        end
    end
    #if i == 1 || i == round((t_end/(plotanzahl-1)))*plotindex || i ==t_end
    #    ausgabe_x[:,ausgabe_index] = pos[:,1]
    #    ausgabe_y[:,ausgabe_index] = pos[:,2]
    #    ausgabe_index = ausgabe_index + 1
    #    plotindex = plotindex + 1
    #end
    v = (v_minus + v_plus).*0.5
    rho = Dichte(pos, m, h)
    P = Druck(rho,rest_density,k)
    f_pres = F_pressure(rho,P,pos,m,h)
    f_vis = F_viscosity(rho,pos,v,m,mu,h)
    #f_gesamt = f_pres + f_vis + f_grav
    f_gesamt = -f_pres + f_grav
    a = f_gesamt ./ rho
    #@show sum(a) #Test der Werte, damit sie gleich bleiben
    if i == time_plot
        x = pos[:,1]
        y= pos[:,2]
        display(scatter(x,y))
        time_plot = time_plot + 20
    end
end

#scatter(ausgabe_x,ausgabe_y,layout=(plotanzeige,plotanzeige))

end # function main()

main()

# hydrogeophysical_process_simulation

In this project a subsurface flow is simulated using the PyGimli Python library. This is a library developed for the geophysical application such as this one (read [here](https://pdf.sciencedirectassets.com/271720/1-s2.0-S0098300417X00101/1-s2.0-S0098300417300584/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEJv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIC%2B6YiOqugQrWmkhVuWyCvglTC9fjTbmm2GSDuBGsVLVAiEAxPmlPeMOBY6uKZhZniLanPb8X201E%2B4ynHnYx9QbScQqvQMI1P%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARADGgwwNTkwMDM1NDY4NjUiDPsLqWI49zcOb0YOrSqRA4MMOF4t98%2Bv0CAro06M9VAaB0rhoabloH13JdSDH0Q4ST%2FGJ1J2aotb1RbmiYNSPhS1QHiM06fqJAvfuw6zpd%2Bpn%2BtAJWB6NYyGa6kyJHHuMbUpVq2%2BrEQFQEPQCwbJ31Xt4k6QK39m0%2Btynoh7rFIxRnxlIfKiMQlNeaGtC5D7Ouh%2Fh59MsSzYQ7a9UOlx98pO5qJB2uFVJnbra3pRn5e08x9%2B8z22Y6Ee1KESvrIbdopUMasX85t7YApganTCFUTUQbZ1p3ZLre%2FeUD3Nq7r9abAuifa2qrZQNNlPSe%2FmkXwtie57FjcnUtB2m0VV3imSAR4XzYtbhGme4HyvtSls7VG0LmQysdOXLTtj4fDXYAlgQvH%2BMhBOwvUebMuYWh4S9ZeL9csjq01u7EM%2FdYs0FH0RML13BlB01OJAKpUoCdN4jruWLsF8ZtCqX%2FulTJBNg03EK7aVwkUrRKNR1MFpkGVkocjMpqJ8hvp9alMUVwcztNpzfTCZsPWFTBUspZctmO84y62oYoVUj81gXmiQML3WoPwFOusBSFLDYvaoJ4eA13YSxggfWsYn1oBk7uBpHA89ZqBII601NQXXeHUE5eSAHZxZuP%2BAQwg%2BW7icuXeuUn2vkAa6ILDr9W1tr58QM8uaiwgYLO7OrJuCPo3slLTivURt5%2Bb1mwkUOxI8S9RKFlUgD5tEqkoKP80wpzY425xgpk09J9ItNtCSLkY3%2FYVbjCnEDZDlGw1cGiUj%2FEmWlYrX9e0%2BWY4uvKb1Xc%2BEgSKWciqwwW7I2Izbo%2Fum4p7i1C11%2Bw9%2B%2BILh9aTYGJwt4iufUkQ%2Fr%2FtG3aG55mlLZtv3C1gnSjalvOWtUgLusFvSKw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20201015T112554Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYR3X7JXNJ%2F20201015%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=20ecadf43a5ebdb8bdae851877925d733aa4da1f0f206fa08cec98b9dafb4705&hash=e1a1a8d09e8cb6acb4a5bc410feb2f089651953f8fe630bc83df320d56602b9b&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0098300417300584&tid=spdf-45c5ca2d-5fd9-48ad-8d79-20680ec048be&sid=a31e3f37415365457c6b5e4088a7ab5e041agxrqb&type=client)). 
Whenever a fluid runs through a porous medium, Darcy's law states the flow rate through the medium as 
<img src="https://render.githubusercontent.com/render/math?math=Q = - K*A* dh/dl> , with
Q: flow rate [m^3/s]
K: hydraulic conductivity [m/s]
A: cross-sectional area of medium [m^2]
dh: hydraulic head [m]
dl: distance of flow [m]

\\\insert pic

Due to this subsurface flow and its contact with the electrically charged rock interfaces, electrical charges in the pore water are draged along with the flow, creating an electric current (also streaming current).
This leads to the generation of an electric field, called the streaming potential $\varphi$.:
$$\nabla \cdot (\sigma \nabla \varphi) = \nabla \cdot (L \textbf{v})$$.
$\varphi$: streaming potential
$\sigma$: electrical condutivity [S/m]
L: dimensionless parameter
h: hydraulic head [m]

\\\insert pic with rock matrix

With this, we have a subsurface flow that generates an electric field. Now we can go forward an create a world with any geometry and simulate a subsurface flow. We can solve the above equations using a Finite Element scheme implemented in the PyGimli library. We could then insert a tracer that acts as a fluid to run through the porous medium.
The transport process of the fluid (advection, diffusion) is governed by the advection-diffusion equation:
$$\frac{\partial c}{\partial t} = \underbrace{\nabla\cdot(D \nabla c)}_{\text{Diffusion / Dispersion}} - \underbrace{\mathbf{v}\nabla c}_{\text{Advection}} + S$$
c: tracer concentration [g/l]
D: dispersion rate [m/s], D = $\alpha \cdot \textbf{v}_{abs}$
v: velocity [m/s]
S: tracer amount [g/ls]

Now, it's possible to translate the simulated fluid concentration into geophysical properties if we assocaite the tracer concentration with salt content and assume dominance of electrolytic conduction. Then, the fluid resistivity $\rho_f$ can be obtained by:
$$1/(0.1*c + \sigma_0)$$, with
c: (salt-)tracer concentration [g/l]
$\sigma_0$: conductivity of groundwater [S/m]

The Archie equation relates the bulk electrical resistivity of a medium $\rho$ to the fluid resistivity $\rho_f$ and saturation $S$ depending on the porosity $\phi$:
$$\rho_b = a \cdot (1/\sigma_W)\cdot\phi^{-m}\cdot S^{-n}$$
a: tortuosity factor (usually a =1)
n: saturation exponent (usually n = 2)
m: cementation exponent
$\phi$: porosity

Now we can simulate the electrical resistivity of the subsurface with electrical resistivity tomography (ERT). In general, for geoelectric measurements, a current is induced in the ground by electrodes. The way to arrange the electrodes can vary. Here we use the dipole-dipole array, which means one pair of electrodes induces the electric current while another pair of electrodes gauges it. The distance between each pair of electrodes is small compared to the distance between the two pairs. This gives us an image of the apparent resitivity values. This is called pseudosection. From here, usually a geophysical inversion is done to obtain the real resistivity values. This part is skipped here.


## Chosen values for this geometry

In my example, an aquifer lies... 

### world dimensions: 
x-direction [m]: 0-34  
y-direction [m]: 20-(-3)  

### hydraulic potential:
h1 [m] = 5  
h2 [m] = 0.2  
hydraulic conductivity K [m/s]:  
world: $10^{-4}$  
circle: $10^{-7}$  
polygon: $10^{-8}$  
boundary conditions: Dirichlet (left: h1, right: h2)  
velocity v [m/s]: $\textbf{v} = \frac{-K \nabla h} {\phi}$  
$\textbf{v}_{abs} = \sqrt{\textbf{v}_x^2 + \textbf{v}_y^2}$  
porosity $\phi$ = 0.3  
### streaming potential:
solve streaming potential equation with Finite Elements for streaming potential $\varphi$ [V]:  
$$\nabla \cdot (\sigma \nabla \varphi) = \nabla \cdot (L \textbf{v})$$  
L = 10  
$\varphi1$ [V] = 0  
$\varphi2$ [V] = 0  
boundary conditions: Dirichlet (left: $\varphi1$, right: $\varphi2$)    
electrical conductivity $\sigma$ [S/m]:  
world: $9 \cdot 10^{-1}$  
circle: $3 \cdot 10^{-4}$  
polygon: $10^{-4}$
### tracer transport:
solve transport equation with Finite Volume for concentration c [g/l]:  
$$\frac{\partial c}{\partial t} = \underbrace{\nabla\cdot(D \nabla c)}_{\text{Diffusion / Dispersion}} - \underbrace{\mathbf{v}\nabla c}_{\text{Advection}} + S$$
3 tracer injection positions S[x,y]:  
position1 = [5,2.5]  
position2 = [2.5,7]  
position3 = [8,6.5]  
injection load [g/l]: 0.001  
dispersion rate D [m/s]:   
$\alpha = 10^{-3}$  
###### timestepping:
steps: 100  
timeperiod [days] = 5  
intervall $\Delta t$ = timeperiod[sec]/(steps-1) = 4363.6 sec  
injections stopped after: steps/2  
##### CFL criterion:  
$\Delta x$ = 0.15 m  
$v_{max}$ = 0.000166 m/s  
c = $\frac{v_{max} \Delta t} {\Delta x}$ = 4.78  
### electrical conductivity:
$\sigma_0$ [S/m] = 0.01  
$\sigma_W =$ 0.1$ c + \sigma_0$  
$\sigma_W$  = 0.01-1000 S/m  
### electrical resitivity:
Archie's equation for $\rho_b$ [$\Omega m$]:  
$$\rho_b = a \cdot (1/\sigma_W)\cdot\phi^{-m}\cdot S^{-n}$$
$\rho_b = 10^{-3}-25  \Omega m$  
a=1  
m=2  
S=1  
n=2  

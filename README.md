# Hydrogeophysical Process Simulation

In this project I simulated a flow in the subsurface using the PyGimli Python library. This is a library developed for the geophysical application such as this one (read [here](https://pdf.sciencedirectassets.com/271720/1-s2.0-S0098300417X00101/1-s2.0-S0098300417300584/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEJv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIC%2B6YiOqugQrWmkhVuWyCvglTC9fjTbmm2GSDuBGsVLVAiEAxPmlPeMOBY6uKZhZniLanPb8X201E%2B4ynHnYx9QbScQqvQMI1P%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARADGgwwNTkwMDM1NDY4NjUiDPsLqWI49zcOb0YOrSqRA4MMOF4t98%2Bv0CAro06M9VAaB0rhoabloH13JdSDH0Q4ST%2FGJ1J2aotb1RbmiYNSPhS1QHiM06fqJAvfuw6zpd%2Bpn%2BtAJWB6NYyGa6kyJHHuMbUpVq2%2BrEQFQEPQCwbJ31Xt4k6QK39m0%2Btynoh7rFIxRnxlIfKiMQlNeaGtC5D7Ouh%2Fh59MsSzYQ7a9UOlx98pO5qJB2uFVJnbra3pRn5e08x9%2B8z22Y6Ee1KESvrIbdopUMasX85t7YApganTCFUTUQbZ1p3ZLre%2FeUD3Nq7r9abAuifa2qrZQNNlPSe%2FmkXwtie57FjcnUtB2m0VV3imSAR4XzYtbhGme4HyvtSls7VG0LmQysdOXLTtj4fDXYAlgQvH%2BMhBOwvUebMuYWh4S9ZeL9csjq01u7EM%2FdYs0FH0RML13BlB01OJAKpUoCdN4jruWLsF8ZtCqX%2FulTJBNg03EK7aVwkUrRKNR1MFpkGVkocjMpqJ8hvp9alMUVwcztNpzfTCZsPWFTBUspZctmO84y62oYoVUj81gXmiQML3WoPwFOusBSFLDYvaoJ4eA13YSxggfWsYn1oBk7uBpHA89ZqBII601NQXXeHUE5eSAHZxZuP%2BAQwg%2BW7icuXeuUn2vkAa6ILDr9W1tr58QM8uaiwgYLO7OrJuCPo3slLTivURt5%2Bb1mwkUOxI8S9RKFlUgD5tEqkoKP80wpzY425xgpk09J9ItNtCSLkY3%2FYVbjCnEDZDlGw1cGiUj%2FEmWlYrX9e0%2BWY4uvKb1Xc%2BEgSKWciqwwW7I2Izbo%2Fum4p7i1C11%2Bw9%2B%2BILh9aTYGJwt4iufUkQ%2Fr%2FtG3aG55mlLZtv3C1gnSjalvOWtUgLusFvSKw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20201015T112554Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYR3X7JXNJ%2F20201015%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=20ecadf43a5ebdb8bdae851877925d733aa4da1f0f206fa08cec98b9dafb4705&hash=e1a1a8d09e8cb6acb4a5bc410feb2f089651953f8fe630bc83df320d56602b9b&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0098300417300584&tid=spdf-45c5ca2d-5fd9-48ad-8d79-20680ec048be&sid=a31e3f37415365457c6b5e4088a7ab5e041agxrqb&type=client)). 
Whenever a fluid runs through a porous medium, Darcy's law states the flow rate through the medium as <br/>
<img src="https://render.githubusercontent.com/render/math?math={Q = - K \cdot A \cdot \frac{dh}{dl}}"> , with<br/>
Q: flow rate [m^3/s]<br/>
K: hydraulic conductivity [m/s]<br/>
A: cross-sectional area of medium [m^2]<br/>
dh: hydraulic head [m]<br/>
dl: distance of flow [m]<br/>
<br/>
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/Darcy.png" width="674" height="383">
<br/>
(image credit: F.Wagner)<br/>

Due to this subsurface flow and its contact with the electrically charged rock interfaces, electrical charges in the pore water are draged along with the flow, creating an electric current (also streaming current).<br/>
This leads to the generation of an electric field, called the streaming potential <img src="https://render.githubusercontent.com/render/math?math={\varphi}">.:<br/>
<img src="https://render.githubusercontent.com/render/math?math={\nabla \cdot (\sigma \nabla \varphi) = \nabla \cdot (L \textbf{v})}">.<br/>
<img src="https://render.githubusercontent.com/render/math?math={\varphi}">: streaming potential<br/>
<img src="https://render.githubusercontent.com/render/math?math={\sigma}">: electrical condutivity [S/m]<br/>
L: dimensionless parameter<br/>
h: hydraulic head [m]<br/>
<br/>
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/SP.png" width="400" height="274">
<br/>
(image credit: Revil et al.(2017))<br/>


With this, we have a subsurface flow that generates an electric field. Now we can go forward an create a world with any geometry and simulate a subsurface flow. We can solve the above equations using a Finite Element scheme implemented in the PyGimli library. We could then insert a tracer that acts as a fluid to run through the porous medium.
The transport process of the fluid (advection, diffusion) is governed by the advection-diffusion equation:<br/>
<img src="https://render.githubusercontent.com/render/math?math={\frac{\partial c}{\partial t} = \underbrace{\nabla\cdot(D \nabla c)}_{\text{Diffusion / Dispersion}} - \underbrace{\mathbf{v}\nabla c}_{\text{Advection}}}"> + S <br/>
c: tracer concentration [g/l]<br/>
D: dispersion rate [m/s], D = <img src="https://render.githubusercontent.com/render/math?math={\alpha \cdot \textbf{v}_{abs}}"><br/>
v: velocity [m/s]<br/>
S: tracer amount [g/ls]<br/>
<br/>
Now, it's possible to translate the simulated fluid concentration into geophysical properties if we assocaite the tracer concentration with salt content and assume dominance of electrolytic conduction. Then, the fluid resistivity <img src="https://render.githubusercontent.com/render/math?math={\rho_f}"> can be obtained by:<br/>
<img src="https://render.githubusercontent.com/render/math?math={\frac{1}{0.1 \cdot c + \sigma_0}}">, with<br/>
c: (salt-)tracer concentration [g/l]<br/>
<img src="https://render.githubusercontent.com/render/math?math={\sigma_0}">: conductivity of groundwater [S/m]<br/>
<br/>
The Archie equation relates the bulk electrical resistivity of a medium <img src="https://render.githubusercontent.com/render/math?math={\rho}"> to the fluid resistivity <img src="https://render.githubusercontent.com/render/math?math={\rho_f}"> and saturation S depending on the porosity <img src="https://render.githubusercontent.com/render/math?math={\phi}">:<br/>
<img src="https://render.githubusercontent.com/render/math?math={\rho_b = a \cdot (1/\sigma_W)\cdot\phi^{-m}\cdot S^{-n}}"><br/>
a: tortuosity factor (usually a =1)<br/>
n: saturation exponent (usually n = 2)<br/>
m: cementation exponent<br/>
<img src="https://render.githubusercontent.com/render/math?math={\phi}">: porosity<br/>
<br/>
Now we can simulate the electrical resistivity of the subsurface with electrical resistivity tomography (ERT). In general, for geoelectric measurements, a current is induced in the ground by electrodes. The way to arrange the electrodes can vary. Here we use the dipole-dipole array, which means one pair of electrodes induces the electric current while another pair of electrodes gauges it. The distance between each pair of electrodes is small compared to the distance between the two pairs. This gives us an image of the apparent resitivity values. This is called pseudosection. From here, usually a geophysical inversion is done to obtain the real resistivity values. This part is skipped here.
<br/>




## Setup and Results

In my example (gravel_aquifer_between_sandstone.ipynb), a gravel layer that lies inbetween two sandstone layers acts as aquifer. This means, its hydraulic conductivity is high and it can bear water well. The sandstone layers have a low hydraulic conductivity and do not carry water as well. Below these three layers follows another layer of gravel. Within the upper gravel layer is a limestone dolomite that again is a low quality aquifer. In the lower gravel layer is a fine sand body that  bears water very well, hence has a high hydraulic conductivity.
In the following, the chosen values for each step of the process are shown, as well as the result of the above mentioned equations.

### world dimensions: 
x-direction [m]: 0-34  <br/>
y-direction [m]: 20-(-3)  <br/>

This image below shows the geometry of the setup 
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/geometry.png" width="674" height="383"><br/>
The mesh quality is shown by red triangles<br/>
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/meshquality.png" width="674" height="383">

### hydraulic potential:
h1 [m] = 5  <br/>
h2 [m] = 0.2  <br/>
world: <img src="https://render.githubusercontent.com/render/math?math={10^{-4}}">  <br/>
rock: <img src="https://render.githubusercontent.com/render/math?math={10^{-7}}">  <br/>
rock2: <img src="https://render.githubusercontent.com/render/math?math={10^{-8}}">  <br/>
boundary conditions: Dirichlet (left: h1, right: h2)  <br/>
velocity v [m/s]: <img src="https://render.githubusercontent.com/render/math?math={\textbf{v} = \frac{-K \nabla h} {\phi}}">  <br/>
<img src="https://render.githubusercontent.com/render/math?math={\textbf{v}_{abs} = \sqrt{\textbf{v}_x^2 + \textbf{v}_y^2}}">  <br/>
porosity <img src="https://render.githubusercontent.com/render/math?math={\phi}"> = 0.3  <br/>
Here is the hydraulic cinductivity for each layer and feature
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/hydr_cond.png" width="674" height="383"> <br/>
and the hydraulic head distribution <br/>
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/hydraulic_grad.png" width="674" height="383">


### streaming potential:
L = 10  <br/>
<img src="https://render.githubusercontent.com/render/math?math={\varphi1}"> [V] = 0  <br/>
<img src="https://render.githubusercontent.com/render/math?math={\varphi2}"> [V] = 0  <br/>
boundary conditions: Dirichlet (left: <img src="https://render.githubusercontent.com/render/math?math={\varphi1, right: \varphi2}">)    <br/>
electrical conductivity <img src="https://render.githubusercontent.com/render/math?math={\sigma}"> [S/m] for each part of the geometry:  <br/>
world: <img src="https://render.githubusercontent.com/render/math?math={9 \cdot 10^{-1}}">  <br/>
circle: <img src="https://render.githubusercontent.com/render/math?math={3 \cdot 10^{-4}}">  <br/>
polygon: <img src="https://render.githubusercontent.com/render/math?math={10^{-4}}"><br/>

The streaming potential of the given geometry can be seen here
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/stream_pot.png" width="674" height="383">

### tracer transport:
3 tracer injection positions S[x,y]:  <br/>
position1 = [5,2.5]  <br/>
position2 = [2.5,7]  <br/>
position3 = [8,6.5]  <br/>
injection load [g/l]: 0.001  <br/>
dispersion rate D [m/s]:   <br/>
<img src="https://render.githubusercontent.com/render/math?math={\alpha = 10^{-3}}">  <br/>

###### timestepping:
steps: 100  <br/>
timeperiod [days] = 5  <br/>
intervall <img src="https://render.githubusercontent.com/render/math?math={\Delta t}"> = timeperiod[sec]/(steps-1) = 4363.6 sec  <br/>
injections stopped after: steps/2  <br/>

##### CFL criterion:  
<img src="https://render.githubusercontent.com/render/math?math={\Delta x}"> = 0.15 m  <br/>
<img src="https://render.githubusercontent.com/render/math?math={v_{max}}"> = 0.000166 m/s  <br/>
c = <img src="https://render.githubusercontent.com/render/math?math={\frac{v_{max} \Delta t} {\Delta x}}"> = 4.78  <br/>

The final simulation of the injected salt tracer concentration after 12, 21.6, 27.6, 36, 60 and 84 hours is shown below from left to right and top down.<br/>
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/conc0.png" width="500" height="300"><img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/conc1.png" width="500" height="300">
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/conc2.png" width="500" height="300"><img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/conc3.png" width="500" height="300">
<img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/conc4.png" width="500" height="300"><img src="https://github.com/MichaelSchffl/hydrogeophysical_process_simulation/blob/master/images/conc5.png" width="500" height="300">

### electrical conductivity:
<img src="https://render.githubusercontent.com/render/math?math={\sigma_0}"> [S/m] = 0.01  <br/>
<img src="https://render.githubusercontent.com/render/math?math={\sigma_W = 0.1 c + \sigma_0}">  <br/>
<img src="https://render.githubusercontent.com/render/math?math={\sigma_W}">  = 0.01-1000 S/m  <br/>

### electrical resitivity:
Archie's equation for <img src="https://render.githubusercontent.com/render/math?math={\rho_b [\Omega m]}">:  <br/>
<img src="https://render.githubusercontent.com/render/math?math={\rho_b = a \cdot (1/\sigma_W)\cdot\phi^{-m}\cdot S^{-n}}"> <br/>
<img src="https://render.githubusercontent.com/render/math?math={\rho_b = 10^{-3}-25  \Omega m}">  <br/>
a=1  <br/>
m=2  <br/>
S=1  <br/>
n=2  <br/>

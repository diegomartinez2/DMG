# 'l' es la distancia del sistema (¿o de cada bloque?)
l =135.775
k(n) = (pi * n) / l
k1(n) = (pi * n)
K2=0.00001
# aquí se ponen las temperaturas de los bloques caliente y frio
Tmax = 400.0
Tmin = 200.0

A(n) = ((1/k1(n))*(-sin(k1(n))*(Tmax+Tmin)))
B(n) = ((1/k1(n))*(-cos(k1(n))*(Tmin-Tmax)-Tmax+Tmin))
#tmax-tmin=delta_t, c(n):delta_t coeff.
C(n)=(-2*B(n)+2*B(n)*(cos(k1(n))))*(1/k1(n))
#delta_t
f(n,x) =exp(-K2*x*(k(n)**2))*(-2*B(n)+2*B(n)*(cos(k1(n))))*(1/k1(n))

#*
#*exp(-K2*k(n)**2*x)

F0(x)= f(1,x) +        f(2,x) +        f(3,x) +        f(4,x) +        f(5,x) +        f(6,x) +        f(7,x) +        f(8,x) +        f(9,x) +        f(10,x) +       f(11,x) +        f(12,x) +        f(13,x) +        f(14,x) +        f(15,x) +        f(16,x) +        f(17,x) +        f(18,x) +        f(19,x) +        f(20,x)

##set xrange [0:2634200]

set fit errorvariables

fit  F0(x) "out1"  u ($1-20000):4 via K2
p 'out1' u ($1-20000):4, F0(x) w l lw 2

pause -1

set print "kbarr.dat"
# Modificar los valores ya que he pasado de 'metal' a 'real' en LAMMPS
print "kbarr ",(K2*1e-5*3*1600*1.38*1e-23)/((5.431*50)*(2*5.431)*(2*5.431)*1e-30)
# No hay que cambiar el 1e-5 para pasar de 'metal' a 'real' en LAMMPS, ya que A˚2/fs→m2/s: $C = \frac{(1 \r{A})^2}{1 \text{ fs}} = \frac{(10^{-10} \text{ m})^2}{10^{-15} \text{ s}} = \frac{10^{-20} \text{ m}^2}{10^{-15} \text{ s}} = \textbf{10^{-5} m^2/s}$
# 1.38*1e-23 es la constante de Boltzmann; 1600 es el número de partículas (y 3 porque es 3D).
#(5.431*50)*(2*5.431)*(2*5.431)*1e-30 es el volumen (con el factor 1e-30 para pasarlo a metros)

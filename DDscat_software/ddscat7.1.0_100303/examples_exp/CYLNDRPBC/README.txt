x=10 calculation for cylinder with |m|=1.33
shpar1=1      : cylinder length along axis
shpar3=1      : cylinder axis || y_TF

shpar2 = (cylinder diameter 2R)/d

Take wavelength = 2*pi units
Cylinder radius = r
Want 2*pi*r/lambda = 10 units
Hence set r=10

take shpar2 = 64.49
The estimated number of dipoles in this slice would be N=pi*(64.49/2)^2=3266.
The actual number turns out to be N=3260.
Thus the "radius" of this slice r=sqrt(3260/pi)d=32.213d
We want the radius to be 10
Thus d=10/32.213

Now need to determine aeff
V=3260*d^3

Thus

aeff=(3*V/(4*pi))^{1/3}=(3*3260/(4*pi))^{1/3}=9.1983*d
    = 9.1983*(10/32.213)=2.8555 units

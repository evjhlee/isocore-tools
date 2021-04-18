# isocore-tools

To get planetary radius use
```
rratio = gcr_to_radius.get_radius(gcr, Mcore, a, tKH, rhocfac)
Rp = Rcore / rratio
```

To get the photospheric radius use
```
Rphot = gcr_to_radius.phot_corr(Rcore/Rp, gcr, Mcore, a, rhocfac)
```

where
```
gcr = gas-to-core mass ratio
Mcore = core mass in Earth masses
Rcore = core radius in Earth radii
a = orbital distance in au
tKH = Kelvin-Helmholtz timescale in Myrs
rhocfac = core density / Earth density
```

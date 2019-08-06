# Testing the combined monochromatic and erf profile.
#
# Now going to compare against full ERA-Interim profile calculated at 60N, 340E
# The full profile is calculated by Pro/MyWave/Stokes_profile/decode_point_spectra_to_stokes_profiles.F
# This code reads GRIB spectra. Python/read_stokes_profile.py plots it.
# Must recode this bit in Julia. Then compare the combined parametric profile to the full profile.
#
# run_stokes.jl generates the profile. This code is only used for plotting

using Stokes, Sphere
using PyPlot
#using Plots
using PyCall
using Printf
using NetCDF
using Dates
basemap = pyimport("mpl_toolkits.basemap")
pyplt = pyimport("matplotlib.pyplot")
mplot3d = pyimport("mpl_toolkits.mplot3d")

# Juno: gcf() to display (in REPL).

function testing(x)
    return x^2
end # function

lon = 340.0
lat = 60.0
infiles = ["../Pro/Stokes_profile/stokes_shear_ei.20100101.nc","../Pro/Stokes_profile/stokes_shear_ei.20100102.nc"]
infiles = ["../Stokesdrift/Stokes_shear_ei/stokes_shear_ei.20100701.nc"] # nice profile
infiles = ["../Pro/Stokes_profile/stokes_shear_ei.20100709.nc"] # weird profile
infiles = ["../Pro/Stokes_profile/stokes_shear_ei.20100708.nc"] # good but weird
infile = infiles[1]
outfile = "o1"
plt = 1
l = 0

b = 1.0
zvec = 0.0:-0.1:-29.9
BIG = 1000.0
miss = 0
TOL = 1.0/BIG
varnames = ["mp1", "ust", "vst", "swh", "shww", "p1ww", "p1ps", "shts", "mdts"]
#v0spdws = v0spd = tm01 = fm01 = hm0 = Vspd = ust = vst = []

# Read lons and lats from first file
lons = ncread(infile, "longitude")
lats = ncread(infile, "latitude")
times = ncread(infile, "time")
timeunits = ncgetatt(infile, "time", "units")
m = match(r"(\d{4})-(\d\d)-(\d\d).(\d\d):(\d\d):(\d\.\d)", timeunits)
t0 = DateTime(m.match, "Y-m-d H:M:S.s")
# CCC fix dates

dlon = Sphere.ang180(lons[2]-lons[1])
dlat = lats[2]-lats[1]
i0 = Int(ceil(Sphere.ang360(lon-lons[1])/dlon))+1
j0 = Int(ceil(lat-lats[1])/dlat)+1
k = 1 # Selecting first time step for now CCC
#println("CCC lon=$lon, lat=$lat, i0=$i0, j0=$j0")

vars = Dict()
dry = []

fout = open(outfile, "w")
#nc = NetCDF.open(infile) # Not needed, and causes erroneous listing of
# contents of file in Juno IDE
for (i, varname) in enumerate(varnames)
    offset = ncgetatt(infile, varname, "add_offset")
    scaling = ncgetatt(infile, varname, "scale_factor")
    ifield = ncread(infile, varname)
    miss = ncgetatt(infile, varname, "missing_value")
    dry = ifield.==miss
    v = ifield*scaling.+offset
    # Reset masked (land) to undef value
    v[dry] .= miss
    vars[varname] = v
end # for
l+=1

### Total sea
tm01 = vars["mp1"]
tm01[tm01.<TOL] .= BIG
fm01 = 1.0./tm01
fm01[fm01.==miss] .= 0.0
hm0 = vars["swh"]
hm0[hm0.==miss] .= 0.0
Vspd = 2π*fm01.*hm0.^2

# Total sea Stokes parameters
Vspd = Stokes.transport(hm0, fm01)
ust = vars["ust"]
vst = vars["vst"]
v0spd = hypot.(ust, vst)

### Swell
# Significant wave height of total swell
shts = vars["shts"]
shts[shts.==miss] .= 0.0
# Swell first moment period
p1ps = vars["p1ps"]
p1ps[p1ps.<TOL] .= BIG
fm01sw = 1.0./p1ps
fm01sw[fm01sw.==miss] .= 0.0
# Wind and wave directions coming from convention, so add 180 deg
mdts = Sphere.ang180(vars["mdts"].+180.0)
mdts[mdts.==miss] .= 361.0

### Swell Stokes parameters
# Swell transport
Vspdsw = Stokes.transport(shts, fm01sw)
Vspdsw[dry] .= 0.0
# Swell wavenumber
ksw = 4π^2*fm01sw.^2/GEARTH
# Swell surface Stokes drift
v0spdsw = 2ksw.*Vspdsw
# ℯ
# Swell Stokes drift direction
sdirsw = mdts
v0eastsw = v0spdsw.*sind.(mdts)
v0northsw = v0spdsw.*cosd.(mdts)
v0eastsw[dry] .= 0.0
v0northsw[dry] .= 0.0
Veastsw = Vspdsw.*sind.(mdts)
Vnorthsw = Vspdsw.*cosd.(mdts)

### Wind sea
# Significant height of wind waves
shww = vars["shww"]
shww[shww.==miss] .= 0.0
# Wind sea first moment period
p1ww = vars["p1ww"]
p1ww[p1ww.<TOL] .= BIG
fm01ws = 1.0./p1ww
fm01ws[fm01ws.==miss] .= 0.0

### Wind sea Stokes parameters
# Wind sea surface Stokes drift
v0eastws = ust-v0eastsw
v0northws = vst-v0northsw
sdirws = rad2deg.(atan.(v0eastws, v0northws))
v0spdws = hypot.(v0eastws,v0northws)
v0spdws[v0eastws.==miss] .= 0.0
# Wind sea Stokes transport
Vspdws = Stokes.transport(shww, fm01ws)
Veastws = Vspdws.*sind.(sdirws)
Vnorthws = Vspdws.*cosd.(sdirws)
# Wind sea wave number
kws = Stokes.phillips_wavenumber(v0spdws, Vspdws, beta=b)

if plt<2
    # print profile in location (lon,lat) to file
    vspdsw = Stokes.mono_profile(v0spdsw[i0,j0,k],ksw[i0,j0,k],zvec)
    veastsw = vspdsw*sind(sdirsw[i0,j0,k])
    vnorthsw = vspdsw*cosd(sdirsw[i0,j0,k])
    vspdws = Stokes.phillips_profile(v0spdws[i0,j0,k],kws[i0,j0,k],zvec)
    veastws = vspdws*sind(sdirws[i0,j0,k])
    vnorthws = vspdws*cosd(sdirws[i0,j0,k])

    # Plot profile in one location
    if plt==1
        #fig=figure()
        fig=matplotlib.pyplot.figure()
        ax = fig.gca(projection="3d")
        plot(veastws,vnorthws,zvec)
        plot(veastws+veastsw,vnorthws+vnorthsw,zvec)
        plot(veastsw,vnorthsw,zvec)
        legendtexts = ("Phillips (wind sea)","Combined", "Monochromatic (swell)")
        #legend(legendtexts,loc="upper left")
        legend(legendtexts,loc="upper right")
        xlabel(L"$u_{east}$ [m/s]")
        ylabel(L"$u_{north}$ [m/s]")
        zlabel(L"Depth [m]")
        profile3dfig = "stokes_combined3d"
        savefig("Fig/$profile3dfig.pdf")
        savefig("Fig/$profile3dfig.png")

        fig2=matplotlib.pyplot.figure()
        plot(veastws,vnorthws)
        plot(veastws+veastsw,vnorthws+vnorthsw)
        plot(veastsw,vnorthsw)
        #legend(legendtexts,loc="upper left")
        legend(legendtexts,loc="upper right")
        xlabel(L"$u_{east}$ [m/s]")
        ylabel(L"$u_{north}$ [m/s]")
        profile2dfig = "stokes_combined2d"
        savefig("Fig/$profile2dfig.pdf")
        savefig("Fig/$profile2dfig.png")
    end # plt
end # if plt<2

# Plot maps
if plt>1
    fig=figure()
    #mp = basemap.Basemap(projection="mill", lon_0=-30)
    mp = basemap.Basemap(projection="ortho", lat_0=30, lon_0=-30, resolution="l")
    mp.drawcoastlines(linewidth=0.5)
    mp.drawcountries(linewidth=0.25)
    mp.fillcontinents(color="coral", lake_color="aqua")
    # Draw the edge of the map projection region (the projection limb)
    mp.drawmapboundary(fill_color="aqua")

    # Draw lat/lon grid lines every 10 degrees.
    mp.drawmeridians(collect(0:10:350))
    mp.drawparallels(collect(-80:10:80))

    # Compute native map projection coordinates of lat/lon grid.
    # First meshgrid the Julian way
    Lons = [lo for lo in lons, la in lats]
    Lats = [la for lo in lons, la in lats]
    # Nasty little bug in Basemap forces me to stretch the vars
    x, y = mp(Lons[:], Lats[:])
    xx = reshape(x, size(Lons))
    yy = reshape(y, size(Lats))
    v = v0spdws[:,:,1]
    v[v.>2] .= 0.0
    mp.pcolor(xx, yy, v[:,:])
    mp.colorbar()
    #gcf()

end # if plt

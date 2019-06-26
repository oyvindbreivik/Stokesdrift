# Comparing the full 2D Stokes drift profile calculated from spectra to the combined profile.
#
# Now going to compare against full ERA-Interim profile calculated at 60N, 340E
# The full profile is calculated by Pro/MyWave/Stokes_profile/decode_point_spectra_to_stokes_profiles.F
# That code reads GRIB spectra. Python/read_stokes_profile.py plots it.
# Must recode this bit in Julia. Then compare the combined parametric profile to the full profile.
#
# run_stokes.jl generates the profile. This code is only used for plotting.

using Stokes, Sphere
using PyPlot
using PyCall
using Printf
#using NetCDF
using Dates
using DelimitedFiles
using CSV

lon = 340.0
lat = 60.0
infiles = ["output_stokes_profile_erai_natl_60N20W_2010.asc","run_stokes_12h.asc"]
outfile = "o1"
plt = 1
l = 0

b = 1.0
zvec = 0.0:-0.1:-29.9
BIG = 1000.0
miss = 0
TOL = 1.0/BIG
#v0spdws = v0spd = tm01 = fm01 = hm0 = Vspd = ust = vst = []

# Read lons and lats from first file

#fout = open(outfile, "w")
i=0
#for row in CSV.File(infiles[2], delim=' ', comment="#")
#    i+=1
#    if i > 12
#        break
#    end
#    println("CCC i, row: ", i, row)
#end
#
for line in eachline(infiles[2])
    print(line)
end

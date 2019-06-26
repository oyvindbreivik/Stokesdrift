using Stokes, Sphere
using Printf
using NetCDF
using Dates
using ArgParse

function parse_commandline(args)
    s = ArgParseSettings("""
    Compute Stokes shear and inverse Stokes depth from NetCDF input. Dump ASCII profile to selected location

    Example:
    time julia run_stokes.jl -i /home/oyvindb/Data/Mdata/Stokes_shear/tst.nc -o stokes_combined.asc --lon 340.0 --lat 60.0 --dep 30 --dz 0.1
    """)

    @add_arg_table s begin
        "--infiles", "-i"
            nargs = '+'
            help = "input files (NetCDF)"
            arg_type = String
        "--outfile", "-o"
            help = "output file (ASCII)"
            default = "stokes_combined.asc"
            #arg_type = String
        "--lon"
            help = "longitude of profile"
            arg_type = Float64
            default = 340.0
        "--lat"
            help = "latitude of profile"
            arg_type = Float64
            default = 60.0
        "--dep"
            help = "depth of profile [m]"
            arg_type = Float64
            default = 30.0
        "--dz"
            help = "depth resolution of profile [m]"
            arg_type = Float64
            default = 0.1
    end

    return parse_args(args, s)
end

function read_stokes_write_combined_profile(infiles, outfile, lon, lat, zvec=0.0:-0.1:-30.0)
    b = 1.0
    #zvec = 0.0:-0.1:-30.0
    BIG = 1000.0
    miss = 0
    TOL = 1.0/BIG
    varnames = ["mp1", "ust", "vst", "swh", "mwd", "shww", "mdww", "p1ww", "p1ps", "shts", "mdts", "wind"]

    # Read lons and lats from first file
    lons = ncread(infiles[1],"longitude")
    lats = ncread(infiles[1],"latitude")

    # Read time units and extract offset (hours since 1900-01-01)
    times = ncread(infiles[1], "time")
    timeunits = ncgetatt(infiles[1], "time", "units")
    m=match(r"(\d{4})-(\d\d)-(\d\d).(\d\d):(\d\d):(\d\.\d)", timeunits)
    t0=DateTime(m.match, "Y-m-d H:M:S.s")

    dlon = Sphere.ang180(lons[2]-lons[1])
    dlat = lats[2]-lats[1]
    i0 = Int(ceil(Sphere.ang360(lon-lons[1])/dlon))+1
    j0 = Int(ceil(lat-lats[1])/dlat)+1

    vars = Dict()
    dry = []
    v0spdws = []

    fout = open(outfile, "w")

    # Loop over files
    for (ifile, infile) in enumerate(infiles)
        times = ncread(infile,"time")
        #nc = NetCDF.open(infile)
        for (i, varname) in enumerate(varnames)
            offset = ncgetatt(infile, varname,"add_offset")
            scaling = ncgetatt(infile, varname,"scale_factor")
            ifield = ncread(infile, varname)
            miss=ncgetatt(infile, varname, "missing_value")
            dry = ifield.==miss
            v = ifield*scaling.+offset
            # Reset masked (land) to undef value
            v[dry] .= miss
            vars[varname] = v
        end

        ### Total sea
        tm01 = vars["mp1"]
        tm01[tm01.<TOL] .= BIG
        fm01 = 1.0./tm01
        fm01[fm01.==miss] .= 0.0
        hm0 = vars["swh"]
        hm0[hm0.==miss] .= 0.0
        mwd = Sphere.ang360(vars["mwd"].+180.0)
        mwd[mwd.==miss] .= 361.0
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
        mdts = Sphere.ang360(vars["mdts"].+180.0)
        mdts[mdts.==miss] .= 361.0

        ### Swell Stokes parameters
        # Swell transport
        Vspdsw = Stokes.transport(shts, fm01sw)
        # Swell wavenumber
        ksw = 4π^2*fm01sw.^2/GEARTH
        # Swell surface Stokes drift
        v0spdsw = 2ksw.*Vspdsw
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
        mdww = Sphere.ang360(vars["mdww"].+180.0)
        mdww[mdww.==miss] .= 361.0

        ### Wind sea Stokes parameters
        # Wind sea surface Stokes drift
        v0eastws = ust-v0eastsw
        v0northws = vst-v0northsw
        sdirws = rad2deg.(atan.(v0eastws, v0northws))
        v0spdws = hypot.(v0eastws, v0northws)
        #show(size(v0spdws))

        # Wind sea Stokes transport
        Vspdws = Stokes.transport(shww, fm01ws)
        Veastws = Vspdws.*sind.(sdirws)
        Vnorthws = Vspdws.*cosd.(sdirws)
        # Wind sea wave number
        kws = Stokes.phillips_wavenumber(v0spdws, Vspdws, beta=b)

        # Wind speed
        wspd = vars["wind"]
        wspd[wspd.==miss] .= 0.0

        ### Dump profile at selected locations
        for (k,t) in enumerate(times)
            # Lat, lon, date
            # Convert to real date and time by adding time units offset to number of hours from start
            t1 = t0+Dates.Hour(t)
            @printf(fout, "# %7.3f %8.4f %s\n", lats[j0], lons[i0], "$t1")
            # Transport header
            @printf(fout, "# Vspd Vspdsw Vspdws [m^2/s] ksw kws [rad/m] hm0 [m] tm01 [s] mwd [deg from N going to] shts p1ps mdts shww p1ww mdww wspd [m/s]\n")
            @printf(fout, "# %11.4e %11.4e %11.4e %11.4e %11.4e %6.2f %6.2f %8.2f %6.2f %6.2f %8.2f %6.2f %6.2f %8.2f %6.2f\n", 
                    Vspd[i0,j0,k], Vspdsw[i0,j0,k], Vspdws[i0,j0,k], ksw[i0,j0,k], kws[i0,j0,k], 
                    hm0[i0,j0,k], tm01[i0,j0,k], mwd[i0,j0,k], shts[i0,j0,k], p1ps[i0,j0,k], mdts[i0,j0,k],
                    shww[i0,j0,k], p1ww[i0,j0,k], mdww[i0,j0,k], wspd[i0,j0,k])

            # Swell profile
            vspdsw = Stokes.mono_profile(v0spdsw[i0,j0,k], ksw[i0,j0,k], zvec)
            veastsw = vspdsw*sind(sdirsw[i0,j0,k])
            vnorthsw = vspdsw*cosd(sdirsw[i0,j0,k])

            # Wind sea profile
            vspdws = Stokes.phillips_profile(v0spdws[i0,j0,k], kws[i0,j0,k], zvec)
            veastws = vspdws*sind(sdirws[i0,j0,k])
            vnorthws = vspdws*cosd(sdirws[i0,j0,k])

            # Loop over vertical
            for (i,z) in enumerate(zvec)
                @printf(fout, "%5.2f %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n", abs(z), veastsw[i]+veastws[i], vnorthsw[i]+vnorthws[i], veastsw[i], vnorthsw[i], veastws[i], vnorthws[i])
            end # for z

            @printf(fout, "\n") # A blank line
        end # for t

    end # for infile

    close(fout)

end # function read_stokes_write_combined_profile

function main()
    parsed_args = parse_commandline(ARGS)
    infiles = parsed_args["infiles"]
    outfile = parsed_args["outfile"]#[1]
    lon = parsed_args["lon"]#[1]
    lat = parsed_args["lat"]#[1]
    zmin = -parsed_args["dep"]#[1]
    dz = -parsed_args["dz"]#[1]
    zvec = 0.0:dz:zmin
    #println("CCC testing $infiles $outfile $lon $lat $zmin $dz")
    #println("CCC zvec: ", zvec)

    #infiles=["../Pro/Stokes_profile/stokes_shear_ei.20100101.nc","../Pro/Stokes_profile/stokes_shear_ei.20100102.nc"]
    read_stokes_write_combined_profile(infiles, outfile, lon, lat, zvec)
end # main

main()

using Stokes, Sphere
using Printf
using NetCDF
using Dates
using ArgParse
using PyPlot
using PyCall
using Statistics

basemap = pyimport("mpl_toolkits.basemap")
pyplt = pyimport("matplotlib.pyplot")
mplot3d = pyimport("mpl_toolkits.mplot3d")

function parse_commandline(args)
    s = ArgParseSettings("""
    Compute Stokes shear and inverse Stokes depth from NetCDF input. Dump ASCII profile to selected location

    Example:
    time julia run_stokes.jl -i /lustre/storeB/project/fou/om/STP40/Stokesdrift/Stokes_shear_ei/stokes_shear_ei.201?0[1-3]*.nc -o run_stokes_12h_JFM.asc --lon 340.0 --lat 60.0 --dep 29.9 --dz 0.1 --strd 2
    time julia run_stokes.jl -i /lustre/storeB/project/fou/om/STP40/Stokesdrift/Stokes_shear_ei/stokes_shear_ei.201?0[7-9]*.nc -o run_stokes_12h_JAS.asc --lon 340.0 --lat 60.0 --dep 29.9 --dz 0.1 --strd 2
    time julia run_stokes.jl -i ../../Stokesdrift/Stokes_shear_ei/stokes*2010*.nc -o run_stokes_12h.asc --lon 340.0 --lat 60.0 --dep 29.9 --dz 0.1 --strd 2
    time julia run_stokes.jl -i /lustre/storeB/project/fou/om/STP40/Stokesdrift/Stokes_shear_ei/stokes_shear_ei.2010*.nc -o run_stokes_12h.asc --lon 340.0 --lat 60.0 --dep 29.9 --dz 0.1 --strd 2
    time julia run_stokes.jl -i /lustre/storeB/project/fou/om/STP40/Stokesdrift/Stokes_shear_ei/stokes_shear_ei.201*.nc -o run_stokes_12h_2010-2011.asc --lon 340.0 --lat 60.0 --dep 29.9 --dz 0.1 --strd 2
    time julia run_stokes.jl -i ../Data/stokes*2010*.nc -o run_stokes_short.asc --lon 340.0 --lat 60.0 --dep 29.9 --dz 0.1 --strd 2

    ### Or from the REPL:
    cd("Stokesdrift")
    empty!(ARGS)
    push!(ARGS, "-i", "../Data/stokes_shear_ei.20100718.nc", "../Data/stokes_shear_ei.20100719.nc", "../Data/stokes_shear_ei.20100720.nc", "../Data/stokes_shear_ei.20100721.nc", "../Data/stokes_shear_ei.20100722.nc", "../Data/stokes_shear_ei.20100723.nc")
    push!(ARGS, "-o", "run_stokes_short.asc", "--lon", "340.0", "--lat", "60.0", "--dep", "29.9", "--dz", "0.1", "--strd", "2")
    include("run_stokes.jl")

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
            default = 29.9
        "--dz"
            help = "depth resolution of profile [m]"
            arg_type = Float64
            default = 0.1
        "--strd"
            help = "stride of output"
            arg_type = Int64
            default = 1
    end

    return parse_args(args, s)
end

"""
Read NetCDF files with wave parameters and compute the profile and the transport under the combined
Stokes profile of a general Phillips-type for wind sea spectrum with a spectral shape parameter `beta`
and a monochromatic profile for the swell to depth `z`.

# Arguments:
 * `infiles`: NetCDF
 * `outfile`: ASCII output file from selected location
 * `zvec`: Depth (negative) below surface [``m``]
 * `lon`: longitude [``deg``]
 * `lat`: latitude [``deg``]

# Returns:
 * `nothing`

2019-07-02
Oyvind.Breivik@met.no
"""
function read_stokes_write_combined_profile(infiles, outfile, lon, lat, zvec=0.0:-0.1:-29.9, strd=1)
    #= * Balancing depth swell v wind sea
       * Transport swell v wind sea
       Both globally and in the grid point of interest
    =#

    b = 1.0
    plt = true
    BIG = 1000.0
    miss = 0
    TOL = 1.0/BIG
    TOL2 = 0.01
    varnames = ["mp1", "ust", "vst", "swh", "mwd", "shww", "mdww", "p1ww", "p1ps", "shts", "mdts", "wind"]

    # Read lons and lats from first file
    lons = ncread(infiles[1], "longitude")
    lats = ncread(infiles[1], "latitude")

    # Read time units and extract offset (hours since 1900-01-01)
    times = ncread(infiles[1], "time")
    timeunits = ncgetatt(infiles[1], "time", "units")
    m = match(r"(\d{4})-(\d\d)-(\d\d).(\d\d):(\d\d):(\d\.\d)", timeunits)
    t0 = DateTime(m.match, "Y-m-d H:M:S.s")

    dummy = ncread(infiles[1], "mp1")
    nx, ny = size(dummy)
    nt = length(times)*length(infiles)
    Vratio = zeros(Float64, nx, ny, nt)
    depthratio = similar(Vratio)
    equaldepth = similar(Vratio)
    vcross = similar(Vratio)

    dlon = Sphere.ang180(lons[2]-lons[1])
    dlat = lats[2]-lats[1]
    i0 = Int(ceil(Sphere.ang360(lon-lons[1])/dlon))+1
    j0 = Int(ceil(lat-lats[1])/dlat)+1

    vars = Dict()
    #dry = 0.0

    fout = open(outfile, "w")

    # Loop over files
    for (ifile, infile) in enumerate(infiles)
        times = ncread(infile, "time")
        for (i, varname) in enumerate(varnames)
            offset = ncgetatt(infile, varname, "add_offset")
            scaling = ncgetatt(infile, varname, "scale_factor")
            ifield = ncread(infile, varname)
            miss=ncgetatt(infile, varname, "missing_value")
            global dry = ifield.==miss
            v = ifield*scaling .+ offset
            # Reset masked (land) to undef value
            v[dry] .= miss
            vars[varname] = v
        end
        #ncclose(infile)
        ncclose()

        ### Total sea
        tm01 = vars["mp1"]
        tm01[tm01.<TOL] .= BIG
        fm01 = 1.0./tm01
        fm01[fm01.==miss] .= 0.0
        hm0 = vars["swh"]
        hm0[hm0.==miss] .= 0.0
        mwd = Sphere.ang360(vars["mwd"].+180.0)
        mwd[mwd.==miss] .= 361.0

        # Total sea Stokes parameters
        Vspd = Stokes.transport(hm0, fm01)
        ust = vars["ust"]
        vst = vars["vst"]
        v0spd = hypot.(ust, vst)
        sdir = atand.(ust, vst)
        ktot = Stokes.phillips_wavenumber(v0spd, Vspd)

        ### Swell

        # Significant wave height of total swell
        shts = vars["shts"]
        shts[shts.==miss] .= 0.0
        # Sometimes there are calm days, even in the North Atlantic
        shts[shts.<0.1] .= 0.1
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
        # sdirsw = mdts
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
        # Sometimes there are calm days, even in the North Atlantic
        shww[shww.<0.1] .= 0.1
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
        sdirws = atand.(v0eastws, v0northws)
        v0spdws = hypot.(v0eastws, v0northws)

        # Wind sea Stokes transport
        Vspdws = Stokes.transport(shww, fm01ws)
        Veastws = Vspdws.*sind.(sdirws)
        Vnorthws = Vspdws.*cosd.(sdirws)

        # Wind sea wavenumber
        kws = Stokes.phillips_wavenumber(v0spdws, Vspdws)

        # Wind sea monochromatic Stokes drift (for balancing depth)
        kws_mono = Stokes.mono_wavenumber(v0spdws, Vspdws)

        # Wind sea monochromatic Stokes drift (for balancing depth)
        v0spdws_mono = Vspdws./kws_mono

        # Total Stokes transport magnitude
        Vspd = hypot.(Veastws+Veastsw, Vnorthws+Vnorthsw)

        println("CCC ifile ", ifile)
        nt0 = length(times)
        k0 = (ifile-1)*nt0+1
        k1 = k0+nt0-1
        Vratio[:,:,k0:k1] = Vspdsw./(Vspd.+TOL)

     ### Alternatively, calculate swell and windsea Stokes surface drift speed given swell and windsea direction

        # Normalized vector components of swell and windsea direction
        sweast = sind.(mdts)
        swnorth = cosd.(mdts)
        wseast = sind.(mdww)
        wsnorth = cosd.(mdww)

        # Calculate magnitude (speed) of Stokes surface swell and windsea vectors
        # Must choose one or the other, depending on the magnitude
        v0spdsw = (wsnorth.*ust-wseast.*vst)./(sweast.*wsnorth-swnorth.*wseast)
        #v0spdsw = (wsnorth.*ust-wseast.*vst)./max.(sweast.*wsnorth-swnorth.*wseast, 0.00001)
        v0spdws = (swnorth.*ust-sweast.*vst)./(wseast.*swnorth-wsnorth.*sweast)
        #println("CCC minimum nom denom ", minimum(wsnorth.*ust-wseast.*vst), minimum(sweast.*wsnorth-swnorth.*wseast))
        lbigswell = v0spdsw .> v0spdws

        # Swell and wind sea Stokes drift components
        v0eastsw = v0spdsw.*sweast
        v0northsw = v0spdsw.*swnorth
        v0eastws = v0spdws.*wseast
        v0northws = v0spdws.*wsnorth

        # Let the Stokes swell component be the difference between the total Stokes drift and the wind sea part ...
        v0eastsw[.!lbigswell] = ust[.!lbigswell] - v0eastws[.!lbigswell]
        v0northsw[.!lbigswell] = vst[.!lbigswell] - v0northws[.!lbigswell]
        # ... except where the swell surface Stokes drift is the bigger of the two
        v0eastws[lbigswell] = ust[lbigswell] - v0eastsw[lbigswell]
        v0northws[lbigswell] = vst[lbigswell] - v0northsw[lbigswell]

        # Fix dry points
        v0eastsw[dry] .= 0.0
        v0northsw[dry] .= 0.0
        v0eastws[dry] .= 0.0
        v0northws[dry] .= 0.0

        ### Recompute surface Stokes drift speed

        # Swell Stokes drift
        sdirsw = atand.(v0eastsw, v0northsw)
        v0spdsw = hypot.(v0eastsw, v0northsw)
        v0spdsw_phil2 = hypot.(v0eastsw, v0northsw) # identical to above, just for book-keeping

        # Wind sea Stokes drift
        sdirws = atand.(v0eastws, v0northws)
        sdirws_phil2 = atand.(v0eastws, v0northws)  # identical to above, just book-keeping
        v0spdws = hypot.(v0eastws, v0northws)
        v0spdws_phil2 = hypot.(v0eastws, v0northws) # identical to above, just pleasing the code above

        # Calculate corresponding wavenumbers
        ksw = Stokes.mono_wavenumber(v0spdsw, Vspdsw)
        ksw_phil2 = Stokes.phillips_wavenumber(v0spdsw_phil2, Vspdsw)

        kws = Stokes.phillips_wavenumber(v0spdws, Vspdws)
        kws_phil2 = Stokes.phillips_wavenumber(v0spdws_phil2, Vspdws) # identical to above ...

        # Wind speed
        wspd = vars["wind"]
        wspd[wspd.==miss] .= 0.0

        ### Compute the ratio of the swell Stokes e-folding depth to the total Stokes e-folding depth

        #tmp = Vspdsw.*v0spd./(Vspd.*v0spdsw.+TOL)
        #tmp = Vspdsw.*v0spdws./(Vspdws.*v0spdsw.+TOL)
        tmp = Vspdsw.*v0spdws_mono./(Vspdws.*v0spdsw.+TOL)
        tmp[dry] .= 0.0
        tmp[abs.(tmp).>25.0] .= 0.0
        depthratio[:,:,k0:k1] = tmp

        ### Compute the normalized cross product of surface swell and wind sea Stokes drift
        tmp = v0spdws.*v0spdsw.*sind.(sdirws-sdirsw)./(v0spd.+TOL2).^2
        tmp[abs.(tmp).>5.0] .= 0.0
        vcross[:,:,k0:k1] = tmp

    ### Compute the balancing depth where swell Stokes drift equals the wind sea Stokes drift
        tmp = zeros(Float64, size(dry))
        #tmp2 = zeros(Float64, size(dry))
        if any(v0spdsw.<0.0) | any(v0spdws.<0.0)
            println("CCC negativity abounds")
        end
        #tmp[.!dry] = log.(v0spdsw[.!dry]./v0spdws[.!dry])./(2(kws[.!dry].-ksw[.!dry]))
        #lsmallws = tmp.>0.0
        # Swell greater than wind sea?
        #tmp2[.!dry] = log.(v0spdws[.!dry]./v0spdsw[.!dry])./(2(ksw[.!dry].-kws[.!dry]))
        #tmp[lsmallws] = tmp2[lsmallws]
        #lundef = tmp.>0
        #println("CCC sum(dry), sum(lsmallws), sum(lundef), size(lundef) ", sum(dry), " ", sum(lsmallws), " ", sum(lundef), " ", prod(size(lundef)))
        #tmp[lundef] .= 0.0
        #equaldepth[:,:,k0:k1] = abs.(tmp)
        #tmp[.!dry] = log.(v0spdws[.!dry]./v0spdsw[.!dry])./(2(kws[.!dry].-ksw[.!dry]))
        tmp[.!dry] = log.(v0spdws_mono[.!dry]./v0spdsw[.!dry])./(2(kws_mono[.!dry].-ksw[.!dry]))
        lundef = tmp.<0.0
        tmp[lundef] .= 0.0
        equaldepth[:,:,k0:k1] = tmp


    ### Dump profile at selected locations, use stride (defaults to 1)

        for (k,t) in enumerate(times[1:strd:end])
            println("CCC equaldepth[i0,j0,(ifile-1)*nt0+k] ", equaldepth[i0,j0,(ifile-1)*nt0+k])
            # Lat, lon, date
            # Convert to real date and time by adding time units offset to number of hours from start
            t1 = t0+Dates.Hour(t)
            @printf(fout, "# %7.3f %8.4f %s\n", lats[j0], lons[i0], "$t1")
            # Transport header
            @printf(fout, "# Vspd Vspdsw Vspdws [m^2/s] ksw kws [rad/m] hm0 [m] tm01 [s] mwd [deg from N going to] shts p1ps mdts shww p1ww mdww wspd [m/s]\n")
            @printf(fout, "# %11.4e %11.4e %11.4e %11.4e %11.4e %6.2f %6.2f %8.2f %6.2f %6.2f %8.2f %6.2f %6.2f %8.2f %6.2f\n",
                    Vspd[i0,j0,k], Vspdsw[i0,j0,k], Vspdws[i0,j0,k], ksw[i0,j0,k], kws[i0,j0,k],
                    hm0[i0,j0,k], tm01[i0,j0,k], mwd[i0,j0,k], shts[i0,j0,k], p1ps[i0,j0,k], sdirsw[i0,j0,k],
                    shww[i0,j0,k], p1ww[i0,j0,k], mdww[i0,j0,k], wspd[i0,j0,k])

            # Total sea Phillips Stokes profile
            vspd_phil = Stokes.phillips_profile(v0spd[i0,j0,k], ktot[i0,j0,k], zvec)
            veast_phil = vspd_phil*sind(sdir[i0,j0,k])
            vnorth_phil = vspd_phil*cosd(sdir[i0,j0,k])

            # Swell monochromatic profile
            vspdsw = Stokes.mono_profile(v0spdsw[i0,j0,k], ksw[i0,j0,k], zvec)
            veastsw = vspdsw*sind(sdirsw[i0,j0,k])
            vnorthsw = vspdsw*cosd(sdirsw[i0,j0,k])

            # Wind sea Phillips profile
            vspdws = Stokes.phillips_profile(v0spdws[i0,j0,k], kws[i0,j0,k], zvec)
            veastws = vspdws*sind(sdirws[i0,j0,k])
            vnorthws = vspdws*cosd(sdirws[i0,j0,k])

            # Phillips swell profile
            vspdsw_phil2 = Stokes.phillips_profile(v0spdsw_phil2[i0,j0,k], ksw_phil2[i0,j0,k], zvec)
            veastsw_phil2 = vspdsw_phil2*sind(sdirsw[i0,j0,k])
            vnorthsw_phil2 = vspdsw_phil2*cosd(sdirsw[i0,j0,k])

            # Wind sea Phillips profile adjusted to swell Phillips profile (the sum must match surface Stokes drift)
            vspdws_phil2 = Stokes.phillips_profile(v0spdws_phil2[i0,j0,k], kws_phil2[i0,j0,k], zvec)
            veastws_phil2 = vspdws_phil2*sind(sdirws_phil2[i0,j0,k])
            vnorthws_phil2 = vspdws_phil2*cosd(sdirws_phil2[i0,j0,k])

            # Loop over vertical
            @printf(fout, "# z veast_comb [m/s]  vnorth_comb       veastsw      vnorthsw       veastws      vnorthws   veast_phil   vnorth_phil veastsw_phil2 vnorthsw_phil2 veastws_phil2 vnorthws_phil2\n")
            for (i,z) in enumerate(zvec)
                @printf(fout, "%5.2f %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n",
                abs(z), veastsw[i]+veastws[i], vnorthsw[i]+vnorthws[i], veastsw[i], vnorthsw[i], veastws[i], vnorthws[i],
                veast_phil[i], vnorth_phil[i], veastsw_phil2[i], vnorthsw_phil2[i], veastws_phil2[i], vnorthws_phil2[i])
            end # for z

            @printf(fout, "\n") # A blank line
        end # for t

    end # for infile

    # Plot maps
    if plt
        Vratiomean = mean(Vratio, dims=3)
        equaldepthmean = mean(equaldepth, dims=3)
        #equaldepthmean = mean(equaldepth[equaldepth.>0.0], dims=3)
        depthratiomean = mean(depthratio, dims=3)
        vcrossmean = mean(vcross, dims=3)
        # Compute native map projection coordinates of lat/lon grid.
        # First meshgrid the Julian way
        Lons = [lo for lo in lons, la in lats]
        Lats = [la for lo in lons, la in lats]

        fig=figure()
        # Colormap
        set_cmap("jet")
        mp = basemap.Basemap(projection="mill",llcrnrlon=0.0,llcrnrlat=-80.0,urcrnrlon=360.0,urcrnrlat=80.0,lon_0=180.0,resolution="c")
        #mp = basemap.Basemap(projection="ortho", lat_0=30, lon_0=-30, resolution="l")
        mp.drawcoastlines(linewidth=0.5)
        mp.drawcountries(linewidth=0.25)
        # Draw the edge of the map projection region (the projection limb)
        mp.drawmapboundary(fill_color="aqua")
        # Draw lat/lon grid lines
        mp.drawmeridians(collect(0:45:360), labels=[true,false,false,true])
        mp.drawparallels(collect(-90:30:90), labels=[true,false,false,false])
        x, y = mp(Lons, Lats)
        #x0, y0 = mp(340.0, 60.0)
        # Something funny going on with the lons and lats - off by one?
        x0, y0 = mp([340.0, 340.0], [60.0, 58.5])
        println("CCC x0 y0", x0[2], y0[2])
        xx = reshape(x, size(Lons))
        yy = reshape(y, size(Lats))
        mp.pcolor(xx, yy, Vratiomean[:,:,1])
        text(x0[2],y0[2],"x",color="red")
        cb = mp.colorbar()
        cb.set_label(L"$V_\mathrm{sw}/V_\mathrm{S}$ [~]")
        title("Ratio of swell to total Stokes transport")
        savefig("Fig/transpratio.png")
        savefig("Fig/transpratio.pdf")

        figure()
        mp2 = basemap.Basemap(projection="mill",llcrnrlon=0.0,llcrnrlat=-80.0,urcrnrlon=360.0,urcrnrlat=80.0,lon_0=180.0,resolution="c")
        mp2.pcolor(xx, yy, equaldepthmean[:,:,1])
        #plt.jet()
        mp2.drawcoastlines(linewidth=0.5)
        mp2.drawcountries(linewidth=0.25)
        #mp2.fillcontinents(color="coral", lake_color="aqua")
        # Draw lat/lon grid lines every 10 degrees.
        mp2.drawmeridians(collect(0:45:360), labels=[true,false,false,true])
        mp2.drawparallels(collect(-90:30:90), labels=[true,false,false,false])
        # Draw the edge of the map projection region (the projection limb)
        mp2.drawmapboundary(fill_color="aqua")
        #clim(0,30)
        clim(0,20)
        #clim(-30,30)
        cb = mp2.colorbar()
        cb.set_label("Balancing depth [m]")
        title("Balancing depth of swell and wind sea Stokes drift")
        savefig("Fig/baldepth.png")
        savefig("Fig/baldepth.pdf")

        fig=figure()
        #mp = basemap.Basemap(projection="mill", lon_0=0)
        mp = basemap.Basemap(projection="mill",llcrnrlon=0.0,llcrnrlat=-80.0,urcrnrlon=360.0,urcrnrlat=80.0,lon_0=180.0,resolution="c")
        #mp = basemap.Basemap(projection="ortho", lat_0=30, lon_0=-30, resolution="l")
        mp.drawcoastlines(linewidth=0.5)
        mp.drawcountries(linewidth=0.25)
        #mp.fillcontinents(color="coral", lake_color="aqua")
        # Draw the edge of the map projection region (the projection limb)
        mp.drawmapboundary(fill_color="aqua")
        # Draw lat/lon grid lines
        mp.drawmeridians(collect(0:45:360), labels=[true,false,false,true])
        mp.drawparallels(collect(-90:30:90), labels=[true,false,false,false])
        x, y = mp(Lons, Lats)
        xx = reshape(x, size(Lons))
        yy = reshape(y, size(Lats))
        mp.pcolor(xx, yy, depthratiomean[:,:,1])
        mp.colorbar()
        cb = mp.colorbar()
        cb.set_label(L"$k_\mathrm{ws}/k_\mathrm{sw}$ [~]")
        title("Ratio of swell to wind sea Stokes e-folding depth")
        savefig("Fig/depthratio.png")
        savefig("Fig/depthratio.pdf")

        figure()
        mp3 = basemap.Basemap(projection="mill",llcrnrlon=0.0,llcrnrlat=-80.0,urcrnrlon=360.0,urcrnrlat=80.0,lon_0=180.0,resolution="c")
        mp3.pcolor(xx, yy, vcrossmean[:,:,1])
        mp3.drawcoastlines(linewidth=0.5)
        mp3.drawcountries(linewidth=0.25)
        #mp2.fillcontinents(color="coral", lake_color="aqua")
        # Draw lat/lon grid lines
        mp3.drawmeridians(collect(0:45:360), labels=[true,false,false,true])
        mp3.drawparallels(collect(-90:30:90), labels=[true,false,false,false])
        # Draw the edge of the map projection region (the projection limb)
        mp3.drawmapboundary(fill_color="aqua")
        #clim(-0.5,0.5)
        clim(-0.1,0.1)
        mp3.colorbar()
        cb = mp3.colorbar()
        cb.set_label(L"$\sin(\Delta\theta) v_\mathrm{ws_0}v_\mathrm{sw_0}v_\mathrm{S_0}^{-2}$ [~]")
        title("Normalized cross product of swell and wind sea Stokes drift")
        savefig("Fig/vcrossmean.png")
        savefig("Fig/vcrossmean.pdf")

        gcf()
    end # if plt

    close(fout)

end # function read_stokes_write_combined_profile

function main()
    parsed_args = parse_commandline(ARGS)
    infiles = parsed_args["infiles"]
    outfile = parsed_args["outfile"]
    lon = parsed_args["lon"]
    lat = parsed_args["lat"]
    zmin = -parsed_args["dep"]
    dz = -parsed_args["dz"]
    strd = parsed_args["strd"]
    zvec = 0.0:dz:zmin
    read_stokes_write_combined_profile(infiles, outfile, lon, lat, zvec, strd)
end # main

main()

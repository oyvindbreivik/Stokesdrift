# Comparing the full 2D Stokes drift profile calculated from spectra to the combined profile.
#
# Now going to compare against full ERA-Interim profile calculated at 60N, 340E
# The full profile is calculated by Pro/MyWave/Stokes_profile/decode_point_spectra_to_stokes_profiles.F
# That code reads GRIB spectra. Python/read_stokes_profile.py plots it.
# Must recode this bit in Julia. Then compare the combined parametric profile to the full profile.
#
# run_stokes.jl generates the profile. This code is only used for plotting.
#
# Todo: add plotting of full profile against combined profile at selected time.

using Statistics
using NumericalIntegration
using Stokes
using PyPlot
using Printf

plotting = true
infiles = ["Stokesdrift/run_stokes_12h.asc", "../Pro/MyWave/Stokes_profile/output_stokes_profile_erai_natl_60N20W_2010.asc"]
NLEV = 300
NPAR_COMB = 7
NPAR_FULL = 5
NPOS = 4

irec = 0
irec0 = 189 # 2010-07-08

# Combined profile v full profile
# Velocity stats
meaneast_comb = []
rmseast_comb = []
stdeast_comb = []
meannorth_comb = []
rmsnorth_comb = []
stdnorth_comb = []
# Speed stats
meanspd_comb = []
rmsspd_comb = []
stdspd_comb = []

# Phllips profile v full profile
# Velocity stats
meaneast_phil = []
rmseast_phil = []
stdeast_phil = []
meannorth_phil = []
rmsnorth_phil = []
stdnorth_phil = []
# Speed stats
meanspd_phil = []
rmsspd_phil = []
stdspd_phil = []


# Open files
f_comb = open(infiles[1])
f_full = open(infiles[2])

zvec = zeros(NLEV)
profile_full = zeros(NLEV, NPAR_FULL)
profile_comb = zeros(NLEV, NPAR_COMB)
# Loop over records in files
while !eof(f_comb) && !eof(f_full)
    ### Profile from combined parametric spectrum
    # Read header lines
    global irec += 1
    ln = readline(f_comb)
    (lat_comb, lon_comb) = [parse(Float64, xs) for xs in split(ln)[2:3]]
    ln = readline(f_comb)
    ln = readline(f_comb)
    (Vspd, mwd) = [parse(Float64, xs) for xs in split(ln)[[2,9]]]

    # Loop over record
    #for k = 1:nz
    for k = 1:NLEV
        ln = readline(f_comb)
        profile_comb[k,:] .= [parse(Float64, xs) for xs in split(ln)]
    end # for
    if length(findall(isnan.(profile_full))) > 0
        println("CCC NaNs phil irec=$irec")
    end
    zvec = -profile_comb[:,1]
    veast_comb = profile_comb[:,2]
    vnorth_comb = profile_comb[:,3]
    vspd_comb = hypot.(veast_comb, vnorth_comb)
    v0spd = vspd_comb[1]

    # Compute Phillips profile from total parameters for comparison with combined and full profiles
    k_phil = Stokes.phillips_wavenumber(v0spd, Vspd)
    vspd_phil = Stokes.phillips_profile(v0spd, k_phil, zvec)
    sdir = atand.(veast_comb[1], vnorth_comb[1])
    veast_phil = vspd_phil*sind(sdir)
    vnorth_phil = vspd_phil*cosd(sdir)

    # Swell and wind sea profiles
    #veastsw = profile_comb[:,4]
    #vnorthsw = profile_comb[:,5]
    #veastws = profile_comb[:,6]
    #vnorthws = profile_comb[:,7]

    # Read last line (z == -30.0) and ignore it
    # CCC This will change
    ln = readline(f_comb)
    # Read separator line before next record
    ln = readline(f_comb)

    ### Profiles from full spectrum
    # Loop over npos locations
    for ipos = 1:NPOS
        # Read header lines
        ln = readline(f_full)
        (lat_full, lon_full) = [parse(Float64, xs) for xs in split(ln)[2:3]]
        ln = readline(f_full)
        ln = readline(f_full)
        # Select location matching that of profile_comb
        if lon_full≈lon_comb && lat_full≈lat_comb
            #println("CCC lon_full, lat_full irec $lon_full $lat_full $irec")
            # Loop over record
            for k = 1:NLEV
                ln = readline(f_full)
                profile_full[k,:] .= [parse(Float64, xs) for xs in split(ln)]
            end # for
            veast_full = profile_full[:,2]
            vnorth_full = profile_full[:,3]
            vspd_full = hypot.(veast_full, vnorth_full)

            # Velocity stats
            push!(meaneast_comb, mean(veast_comb - veast_full))
            push!(rmseast_comb, std(veast_comb - veast_full, mean=0.0))
            push!(stdeast_comb, std(veast_comb - veast_full))
            push!(meannorth_comb, mean(vnorth_comb - vnorth_full))
            push!(rmsnorth_comb, std(vnorth_comb - vnorth_full, mean=0.0))
            push!(stdnorth_comb, std(vnorth_comb - vnorth_full))

            push!(meaneast_phil, mean(veast_phil - veast_full))
            push!(rmseast_phil, std(veast_phil - veast_full, mean=0.0))
            push!(stdeast_phil, std(veast_phil - veast_full))
            push!(meannorth_phil, mean(vnorth_phil - vnorth_full))
            push!(rmsnorth_phil, std(vnorth_phil - vnorth_full, mean=0.0))
            push!(stdnorth_phil, std(vnorth_phil - vnorth_full))

            # Speed stats
            push!(meanspd_comb, mean(vspd_comb - vspd_full))
            push!(rmsspd_comb, std(vspd_comb - vspd_full, mean=0.0))
            push!(stdspd_comb, std(vspd_comb - vspd_full))
            push!(meanspd_phil, mean(vspd_phil - vspd_full))
            push!(rmsspd_phil, std(vspd_phil - vspd_full, mean=0.0))
            push!(stdspd_phil, std(vspd_phil - vspd_full))
        else
            # Skip other locations
            for k = 1:NLEV
                ln = readline(f_full)
            end # for        end # if

        end # if
        if plotting && irec==irec0 # CCC
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
        end # if

        ln = readline(f_full)
        # Read separator line
    end # for ipos

    # CCC Break after 5 while testing
    #if irec > 5
    #    break
    #end
end # while

# Remove NaNs. CCC this will change
deleteat!(meaneast_comb, findall(isnan.(meaneast_comb)))
deleteat!(rmseast_comb, findall(isnan.(rmseast_comb)))
deleteat!(stdeast_comb, findall(isnan.(stdeast_comb)))
deleteat!(meannorth_comb, findall(isnan.(meannorth_comb)))
deleteat!(rmsnorth_comb, findall(isnan.(rmsnorth_comb)))
deleteat!(stdnorth_comb, findall(isnan.(stdnorth_comb)))

deleteat!(meanspd_comb, findall(isnan.(meanspd_comb)))
deleteat!(rmsspd_comb, findall(isnan.(rmsspd_comb)))
deleteat!(stdspd_comb, findall(isnan.(stdspd_comb)))

deleteat!(meaneast_phil, findall(isnan.(meaneast_phil)))
deleteat!(rmseast_phil, findall(isnan.(rmseast_phil)))
deleteat!(stdeast_phil, findall(isnan.(stdeast_phil)))
deleteat!(meannorth_phil, findall(isnan.(meannorth_phil)))
deleteat!(rmsnorth_phil, findall(isnan.(rmsnorth_phil)))
deleteat!(stdnorth_phil, findall(isnan.(stdnorth_phil)))

deleteat!(meanspd_phil, findall(isnan.(meanspd_phil)))
deleteat!(rmsspd_phil, findall(isnan.(rmsspd_phil)))
deleteat!(stdspd_phil, findall(isnan.(stdspd_phil)))

if plotting
    n = length(stdeast_comb)
    bins = -0.05:0.001:0.05
    plusbins = 0:0.001:0.05
    pth = "Fig/"

    # East figures
    xstrvel = "Velocity [m s\$^{-1}\$]"
    xstrspd = "Speed [m s\$^{-1}\$]"
    numstr = "\$n\$ = $n"
    labels = ["Phillips profile", "Combined profile"]

    figure()
    totalmeaneast_phil = mean(meaneast_phil)
    totalmeaneast_comb = mean(meaneast_comb)
    hist(meaneast_phil, bins=bins)
    hist(meaneast_comb, bins=bins)
    xlabel(xstrvel)
    title("Mean diff, east comp, " * @sprintf("Phil: %7.5f ", mean(meaneast_phil)) * @sprintf("Comb: %7.5f", mean(meaneast_comb)) )
    legend(labels)
    gcf()
    fname = pth*"meaneast"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(rmseast_phil, bins=plusbins)
    hist(rmseast_comb, bins=plusbins)
    xlabel(xstrvel)
    title("RMS diff, east comp, " * @sprintf("Phil: %7.5f ", std(rmseast_phil, mean=0.0)) * @sprintf("Comb: %7.5f", std(rmseast_comb, mean=0.0)) )
    legend(labels)
    gcf()
    fname = pth*"rmseast"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(stdeast_phil)
    hist(stdeast_comb)
    xlabel(xstrvel)
    title("Stdev diff, east comp, " * @sprintf("Phil: %7.5f ", std(stdeast_phil)) * @sprintf("Comb: %7.5f", std(stdeast_comb)) )
    legend(labels)
    gcf()
    fname = pth*"stdeast"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(meannorth_phil, bins=bins)
    hist(meannorth_comb, bins=bins)
    xlabel(xstrvel)
    title("Mean diff, north comp, " * @sprintf("Phil: %7.5f ", mean(meannorth_phil)) * @sprintf("Comb: %7.5f", mean(meannorth_comb)) )
    legend(labels)
    gcf()
    fname = pth*"meannorth"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(rmsnorth_phil, bins=plusbins)
    hist(rmsnorth_comb, bins=plusbins)
    xlabel(xstrvel)
    title("RMS diff, north comp, " * @sprintf("Phil: %7.5f ", std(rmsnorth_phil, mean=0.0)) * @sprintf("Comb: %7.5f", std(rmsnorth_comb, mean=0.0)) )
    legend(labels)
    gcf()
    fname = pth*"rmsnorth"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(stdnorth_phil)
    hist(stdnorth_comb)
    xlabel(xstrvel)
    title("Stdev diff, north comp, " * @sprintf("Phil: %7.5f ", std(stdnorth_phil)) * @sprintf("Comb: %7.5f", std(stdnorth_comb)) )
    legend(labels)
    gcf()
    fname = pth*"stdnorth"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(meanspd_phil, bins=bins)
    hist(meanspd_comb, bins=bins)
    xlabel(xstrspd)
    title("Mean diff, speed, " * @sprintf("Phil: %7.5f ", mean(meanspd_phil)) * @sprintf("Comb: %7.5f", mean(meanspd_comb)) )
    legend(labels)
    gcf()
    fname = pth*"meanspd"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(rmsspd_phil, bins=plusbins)
    hist(rmsspd_comb, bins=plusbins)
    xlabel(xstrspd)
    title("RMS diff, speed, " * @sprintf("Phil: %7.5f ", std(rmsspd_phil, mean=0.0)) * @sprintf("Comb: %7.5f", std(rmsspd_comb, mean=0.0)) )
    legend(labels)
    gcf()
    fname = pth*"rmsspd"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    hist(stdspd_phil, bins=plusbins)
    hist(stdspd_comb, bins=plusbins)
    xlabel(xstrspd)
    title("Stdev diff, speed, " * @sprintf("Phil: %7.5f ", std(rmsspd_phil)) * @sprintf("Comb: %7.5f", std(rmsspd_comb)) )
    legend(labels)
    gcf()
    fname = pth*"stdspd"
    savefig(fname*".png")
    savefig(fname*".pdf")

end # if plotting

# Comparing the full 2D Stokes drift profile calculated from spectra to the combined profile.
#
# Now going to compare against full ERA-Interim profile calculated at 60N, 340E
# The full profile is calculated by Pro/MyWave/Stokes_profile/decode_point_spectra_to_stokes_profiles.F
# That code reads GRIB spectra. Python/read_stokes_profile.py plots it.
# Must recode this bit in Julia. Then compare the combined parametric profile to the full profile.
#
# run_stokes.jl generates the profile. This code is only used for plotting.
# Ignore 12:00 UTC as there is a discrepancy between the spectra and the integrated parameters.
#
# Todo: add plotting of full profile against combined profile at selected time.

using Statistics
using NumericalIntegration
using Stokes
using PyPlot
using PyCall
using Printf 

pyplt = pyimport("matplotlib.pyplot")
mplot3d = pyimport("mpl_toolkits.mplot3d")

profileplotting = false
statplotting = false
transplotting = true
transdiffplotting = true
infiles = ["Stokesdrift/run_stokes_12h.asc", "../Pro/MyWave/Stokes_profile/output_stokes_profile_erai_natl_60N20W_2010.asc"]
#infiles = ["Stokesdrift/run_stokes_short.asc", "Stokesdrift/output_stokes_profile_erai_natl_60N20W_short.asc"]
NREC = 365
NLEV = 300
NPAR_COMB = 13
NPAR_FULL = 5
NPOS = 4

irec = 0
irec0 = 189 # 2010-04-05 interesting, wide angle between swell and windsea
irec0 = 189*2 # 2010-07-08 less deviation between swell and windsea
irec0 = 7 # -01-04 Interesting curve, well represented by the combined profile
irec0 = 199*2+1 # -07-18

# Combined profile v full profile
# Velocity stats

meaneast_comb = Float64[]
rmseast_comb = similar(meaneast_comb)
stdeast_comb = similar(meaneast_comb)
meannorth_comb = similar(meaneast_comb)
rmsnorth_comb = similar(meaneast_comb)
stdnorth_comb = similar(meaneast_comb)
# Speed stats
meanspd_comb = similar(meaneast_comb)
rmsspd_comb = similar(meaneast_comb)
stdspd_comb = similar(meaneast_comb)

### Phillips profile vs full profile

# Stokes transport
veast_full_transp = similar(meaneast_comb)
vnorth_full_transp = similar(meaneast_comb)
vspd_full_transp = similar(meaneast_comb)

# Velocity stats
meaneast_phil = similar(meaneast_comb)
rmseast_phil = similar(meaneast_comb)
stdeast_phil = similar(meaneast_comb)
meaneast_comb_phil = similar(meaneast_comb)
rmseast_comb_phil = similar(meaneast_comb)
stdeast_comb_phil = similar(meaneast_comb)
meannorth_phil = similar(meaneast_comb)
rmsnorth_phil = similar(meaneast_comb)
stdnorth_phil = similar(meaneast_comb)
meannorth_comb_phil = similar(meaneast_comb)
rmsnorth_comb_phil = similar(meaneast_comb)
stdnorth_comb_phil = similar(meaneast_comb)

# Speed stats
meanspd_comb = similar(meaneast_comb)
rmsspd_comb = similar(meaneast_comb)
stdspd_comb = similar(meaneast_comb)
meanspd_phil = similar(meaneast_comb)
rmsspd_phil = similar(meaneast_comb)
stdspd_phil = similar(meaneast_comb)
meanspd_comb_phil = similar(meaneast_comb)
rmsspd_comb_phil = similar(meaneast_comb)
stdspd_comb_phil = similar(meaneast_comb)

# Transport stats
nmadspd_comb = similar(meaneast_comb)
nmadspd_phil = similar(meaneast_comb)
nmadspd_comb_phil = similar(meaneast_comb)
spd_comb_diff_transp = similar(meaneast_comb)
spd_phil_diff_transp = similar(meaneast_comb)
spd_comb_phil_diff_transp = similar(meaneast_comb)
east_comb_diff_transp = similar(meaneast_comb)
east_phil_diff_transp = similar(meaneast_comb)
east_comb_phil_diff_transp = similar(meaneast_comb)
north_comb_diff_transp = similar(meaneast_comb)
north_phil_diff_transp = similar(meaneast_comb)
north_comb_phil_diff_transp = similar(meaneast_comb)

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
    global irec += 1 # crazy scope rules
    ln = readline(f_comb)
    (lat_comb, lon_comb) = [parse(Float64, xs) for xs in split(ln)[2:3]]
    date_comb = split(ln)[4]
    ln2 = readline(f_comb)
    ln3 = readline(f_comb)
    (Vspd, Vspdsw, ksw, mwd, p1ps, mdts) = [parse(Float64, xs) for xs in split(ln3)[[2,3,5,9,10,11]]]
    ln4 = readline(f_comb)

    # Loop over record
    for k in 1:NLEV
        ln = readline(f_comb)
        profile_comb[k,:] .= [parse(Float64, xs) for xs in split(ln)]
    end # for

    zvec = -profile_comb[:,1]
    veast_comb = profile_comb[:,2]
    vnorth_comb = profile_comb[:,3]
    vspd_comb = hypot.(veast_comb, vnorth_comb)
    v0spd = vspd_comb[1]
    veast_phil = profile_comb[:,8]
    vnorth_phil = profile_comb[:,9]
    vspd_phil = hypot.(veast_phil, vnorth_phil)
    # Combined Phillips swell and wind sea profiles
    veast_comb_phil = profile_comb[:,10] + profile_comb[:,12]
    vnorth_comb_phil = profile_comb[:,11] + profile_comb[:,13]
    vspd_comb_phil = hypot.(veast_comb_phil, vnorth_comb_phil)

    # Read separator line before next record
    ln = readline(f_comb)

    ### Profiles from full spectrum

    # Loop over npos locations
    for ipos = 1:NPOS
        # Read header lines
        ln = readline(f_full)
        (lat_full, lon_full) = [parse(Float64, xs) for xs in split(ln)[2:3]]
        date_full = split(ln)[4]
        time_full = split(ln)[5]
        ln = readline(f_full)
        ln3f = readline(f_full)
        # Select location matching that of profile_comb, skip 12UTC due to inconsistency with spectra
        if lon_full≈lon_comb && lat_full≈lat_comb && irec%2==1
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

            push!(meaneast_comb_phil, mean(veast_comb_phil - veast_full))
            push!(rmseast_comb_phil, std(veast_comb_phil - veast_full, mean=0.0))
            push!(stdeast_comb_phil, std(veast_comb_phil - veast_full))
            push!(meannorth_comb_phil, mean(vnorth_comb_phil - vnorth_full))
            push!(rmsnorth_comb_phil, std(vnorth_comb_phil - vnorth_full, mean=0.0))
            push!(stdnorth_comb_phil, std(vnorth_comb_phil - vnorth_full))

            # Speed stats
            push!(meanspd_comb, mean(vspd_comb - vspd_full))
            push!(rmsspd_comb, std(vspd_comb - vspd_full, mean=0.0))
            push!(stdspd_comb, std(vspd_comb - vspd_full))
            push!(meanspd_phil, mean(vspd_phil - vspd_full))
            push!(rmsspd_phil, std(vspd_phil - vspd_full, mean=0.0))
            push!(stdspd_phil, std(vspd_phil - vspd_full))
            push!(meanspd_comb_phil, mean(vspd_comb_phil - vspd_full))
            push!(rmsspd_comb_phil, std(vspd_comb_phil - vspd_full, mean=0.0))
            push!(stdspd_comb_phil, std(vspd_comb_phil - vspd_full))

            # Integrate transport, note sign due to ordering of vectors
            push!(veast_full_transp, -integrate(zvec, veast_full))
            push!(vnorth_full_transp, -integrate(zvec, vnorth_full))
            transp = hypot.(veast_full_transp[end], vnorth_full_transp[end])
            push!(vspd_full_transp, transp)

            # Normalized transport speed stats
            #push!(nmadspd_comb, -integrate(zvec, abs.(vspd_comb - vspd_full))/transp)
            #push!(nmadspd_phil, -integrate(zvec, abs.(vspd_phil - vspd_full))/transp)
            #push!(nmadspd_comb_phil, -integrate(zvec, abs.(vspd_comb_phil - vspd_full))/transp)
            push!(nmadspd_comb, -integrate(zvec, hypot.(veast_comb-veast_full, vnorth_comb-vnorth_full))/transp)
            push!(nmadspd_phil, -integrate(zvec, hypot.(veast_phil-veast_full, vnorth_phil-vnorth_full))/transp)
            push!(nmadspd_comb_phil, -integrate(zvec, hypot.(veast_comb_phil-veast_full, vnorth_comb_phil-vnorth_full))/transp)

            # Transport stats, unnormalized
            push!(spd_comb_diff_transp, -integrate(zvec, hypot.(veast_comb-veast_full, vnorth_comb-vnorth_full)))
            push!(spd_phil_diff_transp, -integrate(zvec, hypot.(veast_phil-veast_full, vnorth_phil-vnorth_full)))
            push!(spd_comb_phil_diff_transp, -integrate(zvec, hypot.(veast_comb_phil-veast_full, vnorth_comb_phil-vnorth_full)))
            push!(east_comb_diff_transp, -integrate(zvec, veast_comb-veast_full))
            push!(east_phil_diff_transp, -integrate(zvec, veast_phil-veast_full))
            push!(east_comb_phil_diff_transp, -integrate(zvec, veast_comb_phil-veast_full))
            push!(north_comb_diff_transp, -integrate(zvec, vnorth_comb-vnorth_full))
            push!(north_phil_diff_transp, -integrate(zvec, vnorth_phil-vnorth_full))
            push!(north_comb_phil_diff_transp, -integrate(zvec, vnorth_comb_phil-vnorth_full))

            #push!(veast_comb_diff_transp, -integrate(zvec, abs.(veast_comb - veast_full)))
            #push!(vnorth_comb_diff_transp, -integrate(zvec, abs.(vnorth_comb - veast_full)))
            #push!(veast_phil_diff_transp, -integrate(zvec, abs.(veast_phil - veast_full)))
            #push!(vnorth_phil_diff_transp, -integrate(zvec, abs.(vnorth_phil - veast_full)))
            #push!(veast_comb_phil_diff_transp, -integrate(zvec, abs.(veast_comb_phil - veast_full)))
            #push!(vnorth_comb_phil_diff_transp, -integrate(zvec, abs.(vnorth_comb_phil - veast_full)))

           ### Plot profiles for selected time

            if profileplotting && irec==irec0
                # Swell and wind sea profiles
                veastsw = profile_comb[:,4]
                vnorthsw = profile_comb[:,5]
                veastws = profile_comb[:,6]
                vnorthws = profile_comb[:,7]
                veast = profile_comb[:,8]
                vnorth = profile_comb[:,9]
                veastsw_phil2 = profile_comb[:,10]
                vnorthsw_phil2 = profile_comb[:,11]
                veastws_phil2 = profile_comb[:,12]
                vnorthws_phil2 = profile_comb[:,13]
                #
                # 3D view
                fig=matplotlib.pyplot.figure()
                legendtexts = ("Phillips wind sea", "Combined", "Monochromatic swell", "Full 2D", "Phillips unidir total sea")
                titletext = "Date: $date_comb " * @sprintf("lat: %7.2f, lon: %7.2f", lat_comb, lon_comb)
                ax = fig.gca(projection="3d")
                plot(veastws, vnorthws, zvec)
                plot(veast_comb, vnorth_comb, zvec)
                plot(veastsw, vnorthsw, zvec)
                plot(veast_full, vnorth_full, zvec)
                plot(veast_phil, vnorth_phil, zvec)
                #plot(veastsw_phil2, vnorthsw_phil2, zvec)
                title(titletext)
                legend(legendtexts,loc="center left")
                xlabel(L"East velocity, $v_\mathrm{E}$ [m/s]")
                ylabel(L"North velocity, $v_\mathrm{N}$ [m/s]")
                zlabel(L"$z$ [m]")
                profile3dfig = "stokes_combined3d"
                savefig("Fig/$profile3dfig.pdf")
                savefig("Fig/$profile3dfig.png")
                gcf()

                # 2D bird's eye
                #fig2 = matplotlib.pyplot.figure()
                #plot(veastws, vnorthws)
                #plot(veast_comb, vnorth_comb)
                #plot(veastsw, vnorthsw)
                #plot(veast_full, vnorth_full)
                #plot(veast_phil, vnorth_phil)
                #axis("image")
                ##pyplt.axes().set_aspect("equal", "datalim")
                #title(titletext)
                #legend(legendtexts,loc="upper left")
                #xlabel(L"$u_{east}$ [m/s]")
                #ylabel(L"$u_{north}$ [m/s]")
                #grid(fig2)
                #profile2dfig = "stokes_combined2d"
                #savefig("Fig/$profile2dfig.pdf")
                #savefig("Fig/$profile2dfig.png")
                #gcf()

                # Speed profile
                #fig3 = matplotlib.pyplot.figure()
                #plot(hypot.(veastws, vnorthws), zvec)
                #plot(hypot.(veast_comb, vnorth_comb), zvec)
                #plot(hypot.(veastsw, vnorthsw), zvec)
                #plot(hypot.(veast_full, vnorth_full), zvec)
                #plot(hypot.(veast_phil, vnorth_phil), zvec)
                #title(titletext)
                #legend(legendtexts,loc="lower right")
                #xlabel(L"Speed, $||\mathbf{u}||$ [m/s]")
                #ylabel(L"$z$ [m]")
                #grid(fig3)
                #profilefig = "stokes_combined_speed_profile"
                #savefig("Fig/$profilefig.pdf")
                #savefig("Fig/$profilefig.png")
                #gcf()

                ### Phil2 figs below

                # 2D bird's eye phil2 (two Phillips profiles)
                fig4 = matplotlib.pyplot.figure()
                legendtexts = ("Phillips wind sea", "Combined", "Swell", "Full 2D", "Phillips unidir total sea")
                #legendtexts = ("Phillips (wind sea)", "Combined", "Combined Phillips", "Monochromatic (swell)", "Phillips (swell)", "Full 2D", "Phillips (total sea)")
                plot(veastws_phil2, vnorthws_phil2, color="tab:blue")
                plot(veastsw+veastws, vnorthsw+vnorthws, linestyle="-", color="tab:orange")
                plot(veastsw, vnorthsw, linestyle="-", color="tab:green")
                plot(veast_full, vnorth_full, color="tab:red")
                plot(veast_phil, vnorth_phil, color="tab:purple")
                plot(veastsw_phil2+veastws_phil2, vnorthsw_phil2+vnorthws_phil2, linestyle="-.", color="tab:orange")
                plot(veastsw_phil2, vnorthsw_phil2, linestyle="-.", color="tab:green")
                plot(veastws, vnorthws, linestyle="-.", color="tab:blue")
                axis("image")
                title(titletext)
                legend(legendtexts, loc="upper left")
                xlabel(L"East velocity, $v_\mathrm{E}$ [m/s]")
                ylabel(L"North velocity, $v_\mathrm{N}$ [m/s]")
                grid(fig4)
                profile2dphil2fig = "stokes_combined2d_phil2"
                savefig("Fig/$profile2dphil2fig.pdf")
                savefig("Fig/$profile2dphil2fig.png")
                gcf()

                # Speed profile
                fig5 = matplotlib.pyplot.figure()
                plot(hypot.(veastws_phil2, vnorthws_phil2), zvec, color="tab:blue")
                plot(hypot.(veastsw+veastws, vnorthsw+vnorthws), zvec, linestyle="-", color="tab:orange")
                plot(hypot.(veastsw, vnorthsw), zvec, linestyle="-", color="tab:green")
                plot(hypot.(veast_full, vnorth_full), zvec, color="tab:red")
                plot(hypot.(veast_phil, vnorth_phil), zvec, color="tab:purple")
                plot(hypot.(veastsw_phil2+veastws_phil2, vnorthsw_phil2+vnorthws_phil2), zvec, linestyle="-.", color="tab:orange")
                plot(hypot.(veastsw_phil2, vnorthsw_phil2), zvec, linestyle="-.", color="tab:green")
                plot(hypot.(veastws, vnorthws), zvec, linestyle="-.", color="tab:blue")
                title(titletext)
                legend(legendtexts,loc="lower right")
                xlabel(L"Speed, $||\mathbf{u}||$ [m/s]")
                ylabel(L"$z$ [m]")
                grid(fig5)
                profilephil2fig = "stokes_combined_speed_profile_phil2"
                savefig("Fig/$profilephil2fig.pdf")
                savefig("Fig/$profilephil2fig.png")
                gcf()

                # East profile
                fig6 = matplotlib.pyplot.figure()
                plot(veastws_phil2, zvec, color="tab:blue")
                plot(veastsw+veastws, zvec, linestyle="-", color="tab:orange")
                plot(veastsw, zvec, linestyle="-", color="tab:green")
                plot(veast_full, zvec, color="tab:red")
                plot(veast_phil, zvec, color="tab:purple")
                plot(veastsw_phil2+veastws_phil2, zvec, linestyle="-.", color="tab:orange")
                plot(veastsw_phil2, zvec, linestyle="-.", color="tab:green")
                plot(veastws, zvec, linestyle="-.", color="tab:blue")
                title(titletext)
                legend(legendtexts,loc="lower left")
                xlabel(L"East velocity, $v_\mathrm{E}$ [m/s]")
                ylabel(L"$z$ [m]")
                grid(fig6)
                profilephil2fig = "stokes_combined_east_profile_phil2"
                savefig("Fig/$profilephil2fig.pdf")
                savefig("Fig/$profilephil2fig.png")
                gcf()

                # North profile
                fig7 = matplotlib.pyplot.figure()
                plot(vnorthws_phil2, zvec, color="tab:blue")
                plot(vnorthsw+vnorthws, zvec, linestyle="-", color="tab:orange")
                plot(vnorthsw, zvec, linestyle="-", color="tab:green")
                plot(vnorth_full, zvec, color="tab:red")
                plot(vnorth_phil, zvec, color="tab:purple")
                plot(vnorthsw_phil2+vnorthws_phil2, zvec, linestyle="-.", color="tab:orange")
                plot(vnorthsw_phil2, zvec, linestyle="-.", color="tab:green")
                plot(vnorthws, zvec, linestyle="-.", color="tab:blue")
                title(titletext)
                legend(legendtexts,loc="lower left")
                xlabel(L"North velocity, $v_\mathrm{N}$ [m/s]")
                ylabel(L"$z$ [m]")
                grid(fig7)
                profilephil2fig = "stokes_combined_north_profile_phil2"
                savefig("Fig/$profilephil2fig.pdf")
                savefig("Fig/$profilephil2fig.png")
                gcf()

            end # if
        else
            # Skip other locations
            for k in 1:NLEV
                ln = readline(f_full)
            end # for        end # if

        end # if

        ln = readline(f_full)
        # Read separator line
    end # for ipos

    #= CCC Break after 5 while testing
    if irec > 5
        break
    end
    =#
end # while

if statplotting
    n = length(stdeast_comb)
    bins = -0.05:0.001:0.05
    plusbins = 0:0.001:0.03
    pth = "Fig/"

    # East figures
    xstrvel = "Velocity [m s\$^{-1}\$]"
    xstrspd = "Speed [m s\$^{-1}\$]"
    numstr = "\$n\$ = $n"
    labels = ["Phillips profile", "Combined profile", "Combined Phillips profile"]

    figure()
    hist(meaneast_phil, bins=bins)
    hist(meaneast_comb, bins=bins)
    hist(meaneast_comb_phil, bins=bins)
    xlabel(xstrvel)
    title("Mean diff, east comp, " * @sprintf("Phil: %7.5f ", mean(meaneast_phil)) * @sprintf("Comb: %7.5f ", mean(meaneast_comb)) * @sprintf("Comb Phil: %7.5f", mean(meaneast_comb_phil)) )
    legend(labels)
    gcf()
    fname = pth*"meaneast"
    savefig(fname*".png")
    savefig(fname*".pdf")

    figure()
    ylims = (0, 120)
    hist(rmseast_phil, bins=plusbins)
    hist(rmseast_comb, bins=plusbins)
    hist(rmseast_comb_phil, bins=plusbins)
    xlim(plusbins[1], plusbins[end])
    ylim(ylims...)
    xlabel(xstrvel)
    title("RMS diff, east comp, " * @sprintf("Phil: %7.5f ", mean(rmseast_phil)) * @sprintf("Comb: %7.5f ", mean(rmseast_comb)) * @sprintf("Comb Phil: %7.5f", mean(rmseast_comb_phil)) )
    legend(labels)
    gcf()
    fname = pth*"rmseast"
    savefig(fname*".png")
    savefig(fname*".pdf")

    #=
    figure()
    hist(meannorth_phil, bins=bins)
    hist(meannorth_comb, bins=bins)
    hist(meannorth_comb_phil, bins=bins)
    xlabel(xstrvel)
    title("Mean diff, north comp, " * @sprintf("Phil: %7.5f ", mean(meannorth_phil)) * @sprintf("Comb: %7.5f ", mean(meannorth_comb)) * @sprintf("Comb Phil: %7.5f", mean(meannorth_comb_phil)) )
    legend(labels)
    gcf()
    fname = pth*"meannorth"
    savefig(fname*".png")
    savefig(fname*".pdf")
    =#

    figure()
    hist(rmsnorth_phil, bins=plusbins)
    hist(rmsnorth_comb, bins=plusbins)
    hist(rmsnorth_comb_phil, bins=plusbins)
    xlim(plusbins[1], plusbins[end])
    ylim(ylims...)
    xlabel(xstrvel)
    title("RMS diff, north comp, " * @sprintf("Phil: %7.5f ", mean(rmsnorth_phil)) * @sprintf("Comb: %7.5f ", mean(rmsnorth_comb)) * @sprintf("Comb Phil: %7.5f", mean(rmsnorth_comb_phil)) )
    legend(labels)
    gcf()
    fname = pth*"rmsnorth"
    savefig(fname*".png")
    savefig(fname*".pdf")

    #=
    figure()
    hist(meanspd_phil, bins=bins)
    hist(meanspd_comb, bins=bins)
    hist(meanspd_comb_phil, bins=bins)
    xlabel(xstrspd)
    title("Mean diff, speed, " * @sprintf("Phil: %7.5f ", mean(meanspd_phil)) * @sprintf("Comb: %7.5f ", mean(meanspd_comb)) * @sprintf("Comb Phil: %7.5f", mean(meanspd_comb_phil)) )
    legend(labels)
    gcf()
    fname = pth*"meanspd"
    savefig(fname*".png")
    savefig(fname*".pdf")
    =#

    figure()
    hist(rmsspd_phil, bins=plusbins)
    hist(rmsspd_comb, bins=plusbins)
    hist(rmsspd_comb_phil, bins=plusbins)
    xlim(plusbins[1], plusbins[end])
    ylim(ylims...)
    xlabel(xstrspd)
    title("RMS diff, speed, " * @sprintf("Phil: %7.5f ", mean(rmsspd_phil)) * @sprintf("Comb: %7.5f ", mean(rmsspd_comb)) * @sprintf("Comb Phil: %7.5f", mean(rmsspd_comb_phil)) )
    legend(labels)
    gcf()
    fname = pth*"rmsspd"
    savefig(fname*".png")
    savefig(fname*".pdf")

end # if statplotting

if transplotting

    pth = "Fig/"
    xstrtransp = "Normalized transport"
    textpos = (0.05, 80.0)
    textpos2 = (0.75, 45.0)
    col = "k"
    ylims = (0, 100)
    transbins = 0:0.1:1.6
    figtransp = matplotlib.pyplot.figure()
    subplot(311)
    hist(nmadspd_phil, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    println("NMAD transport, " * @sprintf("Phil: %7.5f ", mean(nmadspd_phil)) * @sprintf("Comb: %7.5f ", mean(nmadspd_comb)) * @sprintf("Comb Phil: %7.5f", mean(nmadspd_comb_phil)) )
    title(L"Normalized difference from ERA-I profiles, $\delta V = V^{-1} \int_{-30 m}^0 \, |v_{mod}-v| \, dz$", fontsize=12)
    text(textpos..., "(a) Phillips unidirectional profile")
    text(textpos2..., "Normalized mean error: " * @sprintf("%7.5f ", mean(nmadspd_phil)))

    subplot(312)
    hist(nmadspd_comb, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    ylabel("Number of occurrences")
    text(textpos..., "(b) Phillips (wind sea) and monochromatic (swell) directional profile")
    text(textpos2..., "Normalized mean error: " * @sprintf("%7.5f ", mean(nmadspd_comb)))

    subplot(313)
    hist(nmadspd_comb_phil, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    text(0.3, textpos[2], "(c) Phillips wind sea and swell directional profile")
    text(textpos2..., "Normalized mean error: " * @sprintf("%7.5f ", mean(nmadspd_comb_phil)))
    xlabel(xstrtransp)
    gcf()
    fname = pth*"nmadspd"
    savefig(fname*".png")
    savefig(fname*".pdf")

end # if transplotting

if transdiffplotting
    pth = "Fig/"
    textpos = (0.00, 85.0)
    textpos3 = (0.12, 85.0)
    xstrtransp = L"Transport [m$^2/$s]"
    col = "k"
    ylims = (0, 100)
    transbins = 0:0.03:0.8
    figtransp = matplotlib.pyplot.figure()
    subplot(311)
    hist(spd_phil_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    #println("Transport stats, " * @sprintf("Phil: %7.5f ", mean(spd_phil_diff_transp)) * @sprintf("Comb: %7.5f ", mean(spd_comb_diff_transp)) * @sprintf("Comb Phil: %7.5f", mean(spd_comb_phil_diff_transp)) )
    statstxt="Transport stats:\n " * @sprintf("Phil: %7.5f ", mean(spd_phil_diff_transp)) * @sprintf("\n Comb: %7.5f ", mean(spd_comb_diff_transp)) * @sprintf("\n Comb Phil: %7.5f", mean(spd_comb_phil_diff_transp))
    println(statstxt)
    title(L"Difference from ERA-I profiles, $\Delta V = \int_{-30 m}^0 \, |v_{mod}-v| \, dz$", fontsize=12)
    text(textpos..., "(a) Phillips unidirectional profile")
    text(0.5, 30, "Mean error: " * @sprintf("%7.5f ", mean(spd_phil_diff_transp)))

    subplot(312)
    hist(spd_comb_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    ylabel("Number of occurrences")
    text(textpos..., "(b) Phillips (wind sea) and monochromatic (swell) directional profile")
    text(0.5, 30, "Mean error: " * @sprintf("%7.5f ", mean(spd_comb_diff_transp)))

    subplot(313)
    hist(spd_comb_phil_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    text(textpos3..., "(c) Phillips wind sea and swell directional profile")
    text(0.5, 30, "Mean error: " * @sprintf("%7.5f ", mean(spd_comb_phil_diff_transp)))
    xlabel(xstrtransp)
    gcf()
    fname = pth*"transpdiff"
    savefig(fname*".png")
    savefig(fname*".pdf")

    ### East transport
    xstrtransp = L"East transport [m$^2/$s]"
    col = "k"
    ylims = (0, 100)
    transbins = 0:0.03:0.8
    figtransp = matplotlib.pyplot.figure()
    subplot(311)
    hist(east_phil_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    title(L"Difference from ERA-I profiles, $\Delta V_E = \int_{-30 m}^0 \, (v_{mod,E}-v_E) \, dz$", fontsize=12)
    text(textpos..., "(a) Phillips unidirectional profile")
    text(0.5, 30, "Mean abs error: " * @sprintf("%7.5f ", mean(abs.(east_phil_diff_transp))))

    subplot(312)
    hist(east_comb_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    ylabel("Number of occurrences")
    text(textpos..., "(b) Phillips (wind sea) and monochromatic (swell) directional profile")
    text(0.5, 30, "Mean abs error: " * @sprintf("%7.5f ", mean(abs.(east_comb_diff_transp))))

    subplot(313)
    hist(east_comb_phil_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    text(textpos3..., "(c) Phillips wind sea and swell directional profile")
    text(0.5, 30, "Mean abs error: " * @sprintf("%7.5f ", mean(abs.(east_comb_phil_diff_transp))))
    xlabel(xstrtransp)
    gcf()
    fname = pth*"east_transpdiff"
    savefig(fname*".png")
    savefig(fname*".pdf")

    ### North transport
    xstrtransp = L"North transport [m$^2/$s]"
    col = "k"
    ylims = (0, 100)
    transbins = 0:0.03:0.8
    figtransp = matplotlib.pyplot.figure()
    subplot(311)
    hist(north_phil_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    title(L"Difference from ERA-I profiles, $\Delta V_E = \int_{-30 m}^0 \, (v_{mod,N}-v_N) \, dz$", fontsize=12)
    text(textpos..., "(a) Phillips unidirectional profile")
    text(0.5, 30, "Mean abs error: " * @sprintf("%7.5f ", mean(abs.(north_phil_diff_transp))))

    subplot(312)
    hist(north_comb_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    ylabel("Number of occurrences")
    text(textpos..., "(b) Phillips (wind sea) and monochromatic (swell) directional profile")
    text(0.5, 30, "Mean abs error: " * @sprintf("%7.5f ", mean(abs.(north_comb_diff_transp))))

    subplot(313)
    hist(north_comb_phil_diff_transp, bins=transbins, color=col)
    xlim(transbins[1], transbins[end])
    ylim(ylims...)
    text(textpos3..., "(c) Phillips wind sea and swell directional profile")
    text(0.5, 30, "Mean abs error: " * @sprintf("%7.5f ", mean(abs.(north_comb_phil_diff_transp))))
    xlabel(xstrtransp)
    gcf()
    fname = pth*"north_transpdiff"
    savefig(fname*".png")
    savefig(fname*".pdf")

end

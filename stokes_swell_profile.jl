using Stokes
using ArgParse
#using PyPlot # uncomment when needed

function parse_commandline()
    s = ArgParseSettings("""
    Compute Stokes shear and inverse Stokes depth from NetCDF input

    Example:
    time julia stokes_shear.jl -i /home/oyvindb/Data/Mdata/Stokes_shear/tst.nc -o kavg
    time julia stokes_shear.jl -i /home/oyvindb/Data/Mdata/Stokes_shear/stokes_shear.20100101.nc -o kavg
    time julia stokes_shear.jl -i /media/Elements/Stokes_shear2/stokes_shear.2010010?.nc -o kavg
    time julia stokes_swell_profile.jl -i /home/oyvindb/Huge2/Stokes_shear_ei/tstcombined.nc -o kavg
    """)

    @add_arg_table s begin
        "--infiles", "-i"
            nargs = '+'
            help = "input files (NetCDF)"
        "--outfilepostfix", "-o"
            help = "output file postfix (NetCDF)"
            nargs = 1
        "--plot", "-p"
            help = "plot figures"
            action = :store_true
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    plot = parsed_args["plot"]
    infiles = parsed_args["infiles"]
    outfilepostfix = parsed_args["outfilepostfix"][1]

    Stokes.read_stokes_write_combined_profile(infiles, outfile, plot)
end

main()


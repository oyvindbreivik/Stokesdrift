"""
Stokes.jl - A toolbox of Stokes drift functions.

Module Stokes contains a number of functions that return (vectorized) Stokes
drift variables such as the Stokes profile under a Phillips spectrum,
the transport to a given depth from the surface, and the inverse depth scale.

2018-09-01 - Conversion to Julia 1.0

Oyvind.Breivik@met.no
"""
module Stokes

using PyPlot, PyCall
using Printf
import SpecialFunctions.erfc
using Sphere

export fp_pm, hm0_pm
export surface_stokes_rascle09
export transport, phillips_wavenumber, phillips_profile
export phillips_shear, phillips_transportz

const GEARTH = 9.806
export GEARTH

"""
The peak frequency of the Pierson-Moskowitz spectrum.

# Argument:
 * `u10`: 10 m wind speed [``m/s``]

# Returns:
 * `fp`: Peak frequency [``Hz``]

2016-05-27
Oyvind.Breivik@met.no
"""
fp_pm(u10) = 0.877GEARTH./(2pi*u10)


"""
The significant wave height of the Pierson-Moskowitz spectrum.

# Argument
 * `u10`: 10 m wind speed [``m/s``]

# Returns
 * `hm0`: Significant wave height [``m``]

# Reference
 .. [1] WMO (1998): Guide to wave analysis and forecasting, report no 702.

2016-05-11
Oyvind.Breivik@met.no
"""
hm0_pm(u10) = 0.0246u10.^2


"""
The surface Stokes drift velocity formula of Ardhuin et al (2009).

# Argument:
 * `u10`: 10 m wind speed [``m/s``]
 * `hm0`: Significant wave height [``m``]
          [optional, set to negative if missing]. Uses hm0_pm if not provided
 * `fc`: Cutoff frequency [``Hz``]
          [default 0.5]

# Returns:
 * `v0`: Surface Stokes drift speed [``m/s``]

Reference:
Ardhuin, F et al (2009). Observation and estimation of Lagrangian, Stokes
and Eulerian currents induced by wind and waves at the sea surface, J Phys
Oceanogr, doi:10.1175/2009JPO4169.1

2016-05-11
Oyvind.Breivik@met.no
"""
function surface_stokes_ardhuin09(u10, hm0=-1, fc=0.5)
    if hm0 .< 0
        hm0 = hm0_pm(u10)
    end
    return 5.0e-4*(1.25 .- 0.25*(0.5./fc).^1.3).*u10.*min.(u10,14.5) .+ 0.027*(hm0 .- 0.4)
end


"""
Compute the transport from the significant wave height and the mean frequency
to infinite depth. All wave energy assumed to propagate in the same direction.

# Arguments:
 * `hm0`: Significant wave height [``m``]
 * `fm01`: Mean frequency [``Hz``]

# Returns:
 * `V`: Stokes transport speed [``m^2/s``]

2016-04-25
Oyvind.Breivik@met.no
"""
transport(hm0, fm01) = 2pi*fm01.*hm0.^2/16


"""
Compute the wavenumber under a monochromatic wave.

# Arguments:
 * `v0`: Surface Stokes drift speed [m/s]
 * `V`: Stokes transport speed [``m^2/s``]

# Returns:
 * `k`: Wavenumber [``rad/m``]

2016-04-25
Oyvind.Breivik@met.no
"""
mono_wavenumber(v0, V) = v0./2V


"""
Compute the inverse depth scale or wavenumber of a Phillips-type spectrum

# Arguments:
 * `v0`: Surface Stokes drift speed [m/s]
 * `V`: Stokes transport speed [``m^2/s``]
 * `beta`: Non-dimensional spectral shape parameter [~]

# Returns:
 * `k`: Inverse depth scale (or wavenumber) [rad/m]

2016-04-25
Oyvind.Breivik@met.no
"""
function phillips_wavenumber(v0, V; beta=1)
    k = v0.*(1 .- 2beta/3)./(2V)
end


"""
Compute the Stokes profile under a monochromatic wave.

# Arguments:
 * `v0`: Surface Stokes drift speed [m/s]
 * `k`: wavenumber [``rad/m``]
 * `z`: Depth (negative) below surface [``m``]

# Returns:
 * `v`: The Stokes profile [m/s]

2016-05-11
Oyvind.Breivik@met.no
"""
mono_profile(v0, k, z) = v0.*exp.(2k.*z)


"""
Compute the Stokes profile under a Phillips-type spectrum.

# Arguments:
 * `v0`: Surface Stokes drift speed [m/s]
 * `k`: Inverse depth scale (or wavenumber) [``rad/m``]
 * `z`: Depth (negative) below surface [``m``]
 * `beta`: Non-dimensional spectral shape parameter [``~``]

# Returns:
 * `v`: The Stokes profile [m/s]

2016-04-25
Oyvind.Breivik@met.no
"""
function phillips_profile(v0, k, z; beta=1)
    v0.*(exp.(2k.*z) - beta.*sqrt.(-2pi*k.*z).*erfc.(sqrt.(-2k.*z)))
end


"""
Compute the shear of the Stokes profile under the general Phillips-type
spectrum with a spectral shape parameter beta.

# Arguments:
 * `v0`: Surface Stokes drift speed [m/s]
 * `k`: Inverse depth scale (or wavenumber) [``rad/m``]
 * `z`: Depth (negative) below surface [``m``]
 * `beta`: Non-dimensional spectral shape parameter [~]

# Returns:
 * `shear`: The shear of the Stokes profile [1/s]

2016-04-25
Oyvind.Breivik@met.no
"""
function phillips_shear(v0, k, z; beta=1)
    shear = beta*v0.*sqrt(-pi*k./2z).*erfc.(sqrt(-2k.*z))
    if beta != 1
        shear += v0.*(2*(1-beta)*k.*exp.(2k.*z))
    end
    return shear
end


"""
Compute the transport under the Stokes profile of a general Phillips-type
spectrum with a spectral shape parameter `beta` to depth `z`.

# Arguments:
 * `v0`: Surface Stokes drift speed [m/s]
 * `k`: Inverse depth scale (or wavenumber) [``rad/m``]
 * `z`: Depth (negative) below surface [``m``]
 * `beta`: Non-dimensional spectral shape parameter [~]

# Returns:
 * `V`: Stokes transport speed [``m^2/s``]

2016-04-25
Oyvind.Breivik@met.no
"""
function phillips_transportz(v0, k, z; beta=1)
    (v0./2k).*(1 - exp.(2k.*z) -
      (2beta/3)*(1+sqrt(pi)*(-2k.*z).^1.5.*erfc.(sqrt(-2k.*z)) -
       (1-2k.*z).*exp.(2k.*z)))
end


"""
Computes the Stokes drift velocity vector from a 2D spectrum

# Arguments:
 * `spec2d`: Spectral density array nang x nfre [m^2/Hz/rad]
 * `th`: Direction array of length nang clockwise from north [deg]
 * `fr`: Frequency array of length nfre [Hz]
 * `z`: Vertical coordinate (positive upwards, ie negative below surface) [m]

# Returns:
 * `ust`: East component of Stokes drift velocity, array of length `z` [``m/s``]
 * `vst`: North component of Stokes drift velocity, array of length `z` [``m/s``]
 * `uhf`: East component of high frequency contribution to Stokes drift velocity, array of length `z` [``m/s``]
 * `vhf`: North component of high frequency contribution to Stokes drift velocity, array of length `z` [``m/s``]
 * `utransp`: East component of Stokes drift transport [``m^2/s``]
 * `vtransp`: North component of Stokes drift transport [``m^2/s``]
 * `utransphf`: East component of high frequency Stokes drift transport [``m^2/s``]
 * `vtransphf`: North component of high frequency Stokes drift transport [``m^2/s``]

# See also:
 * `momentspec2d`

"""
function stokesprofile2d(spec2d, th, fr, zlevs=[0,])

    FAC = 16pi^3/GEARTH

    # High-frequency cutoff
    fc = fr[end]
    fachf = FAC*fc^5
    delth = deg2rad.(minimum(diff(th)))
    theta = deg2rad.(th)

    dfr = Sphere.deltafr(fr)
    fact = FAC*fr.^3 .* dfr*delth

    xk = 4pi*fr.^2/GEARTH
    exp2k = exp.(2xk*zlevs')

    ust = zeros(zlevs)
    vst = zeros(zlevs)

    # Integrate Stokes drift vector for all vertical levels
    for (l,z) in enumerate(zlevs)
        ust[l] = sum(fact.*exp.(2xk*z)*sind.(th)'.*spec2d)
        vst[l] = sum(fact.*exp.(2xk*z)*cosd.(th)'.*spec2d)
    end

    # Add high-frequency tail
    mu = -8pi^2*zlevs/GEARTH
    tailint = exp.(-mu*fc^2)/fc - sqrt.(pi*mu)*(1.0-erf(fc*sqrt.(mu)))
    sinint = (spec2d[:,end]'*sind.(th))[1]
    cosint = (spec2d[:,end]'*cosd.(th))[1]
    uhf = fachf*sinint*tailint
    vhf = fachf*cosint*tailint
    utransphf = 2pi*sinint*(fc^2)/3
    vtransphf = 2pi*cosint*(fc^2)/3

    ust = ust+uhf
    vst = vst+vhf

    utransp = (2pi*sind.(th)'*spec2d*(fr.*dfr)*delth)[1]
    vtransp = (2pi*cosd.(th)'*spec2d*(fr.*dfr)*delth)[1]

    utransp = utransp+utransphf
    vtransp = vtransp+vtransphf

    return ust, vst, uhf, vhf, utransp, vtransp, vtransphf, utransphf
end


# Private functions
"""
Three approximate transport profiles and one exact.
"""
function transporty(y; b=1)
    c = 2b/3
    return 1-c - (1-c)*exp.(-y) + c*y.*exp.(-y) - c*sqrt.(pi)*y.^1.5.*erfc.(sqrt.(y))
end

function approx_transporty0(y; b=1)
    c = 2b/3
    return 1-c - (1-c)*exp.(-y)
end

function approx_transporty1(y; b=1)
    c = 2b/3
    return 1-c - (1-c)*exp.(-y) + c*y.*exp.(-y)
end

function approx_transporty2(y; b=1)
    c = 2b/3
    return 1-c - (1-c)*exp.(-y) - c*sqrt.(pi)*y.^1.5.*erfc.(sqrt.(y))
end

end # module

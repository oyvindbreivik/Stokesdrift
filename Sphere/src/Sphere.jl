"""
sphere.jl - A toolbox of spherical functions.

12 Aug 99, Oyvind.Breivik@nrsc.no.
2012-06-10, Converted to Python 2.7 by oyvind.breivik@ecmwf.int
2016-05-02 Converted to Julia v0.4, oyvind.breivik@met.no
2018-09-02 Converted to Julia v1.0, oyvind.breivik@met.no
"""
module Sphere
export ang360, ang180, ang2pi, angpi, spheredist, haversine, spherearcrad
export spherepos, spheredir

### Constants
REARTH = 6371001 # radius of the earth [m]
EPS = 1.0E-10

"""
    dec2sex(dec)

Convert from decimal degrees to sexagesimal (degrees, minutes and decimal seconds).

# Arguments:
 * `dec`: decimal degrees [``deg``]

# Returns:
 * `deg`: degrees [``deg``]
 * `min`: minutes [``[0,60)``]
 * `sec`: seconds [``[0,60)``]
 
# History:
 * 2005-03-14, Matlab, Oyvind.Breivik@met.no.
 * 2019-07-30, Julia v1.0, Oyvind.Breivik@met.no.
"""
function dec2sex(dec) 
    deg = trunc(dec)
    min = abs(dec-deg)*60
    sec = (min-trunc(min))*60
    min = trunc(min)

    return deg, min, sec
end # function

"""
    sex2dec(deg, min=0, sec=0.0)

Convert from the sexagesimal (degrees, minutes and decimal seconds) to decimal degrees.

# Arguments:
 * `deg`: degrees [``deg``]
 * `min`: minutes [``[0,60)``]
 * `sec`: seconds [``[0,60)``]

# Returns:
 * `dec`: decimal degrees [``deg``]
 
# History:
 * 2005-03-14, Matlab, Oyvind.Breivik@met.no.
 * 2019-07-30, Julia v1.0, Oyvind.Breivik@met.no.
"""
sex2dec(deg, min=0, sec=0.0) = deg + sign(deg)*min/60.0 + sign(deg)*sec/3600.0

"""
    ang360(ang)

Maps an array of angles [deg] to [0,360).
"""
ang360(ang) = mod.(ang, 360)


"""
    ang180(ang)

Maps an array of angles [deg] to [-180,180).
"""
ang180(ang) = mod.(ang.+180, 360) .- 180


"""
    ang2pi(theta)

Maps an array of angles [rad] to [0,2π).
"""
ang2pi(theta) = mod.(theta, 2π)


"""
    angpi(theta)

Maps an array of angles [deg] to [-180,180).
"""
angpi(theta) = mod.(theta.+π, 2π) .- π


"""
    spheredist(lon1, lat1, lon2, lat2; deg=false)

Returns the great circle (geodesic) distance measured in meters
between positions (lon1,lat1) and (lon2,lat2) [deg].

If deg=true, the distance in degrees is returned.

Arguments can be vectorized.

Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
Section 4.3.149.

Adopted from spherdist.

Requires haversine.

"""
function spheredist(lon1, lat1, lon2, lat2; deg=false)
    rlon1 = deg2rad.(lon1)
    rlat1 = deg2rad.(lat1)
    rlon2 = deg2rad.(lon2)
    rlat2 = deg2rad.(lat2)

    if deg
        return rad2deg.(haversine(rlon1,rlat1,rlon2,rlat2))
    else
        return REARTH*haversine(rlon1,rlat1,rlon2,rlat2)
    end
end # function

"""
    haversine(rlon1,rlat1,rlon2,rlat2)

Returns the great circle (geodesic) distance measured in radians [0,π)
between positions (rlon1,rlat1) and (rlon2,rlat2) [rad] using the haversine
formula.

    c = haversine(rlon1,rlat1,rlon2,rlat2)

Reference: http://en.wikipedia.org/wiki/Haversine_formula

"""
function haversine(rlon1,rlat1,rlon2,rlat2)
    dlat = rlat2 - rlat1  
    dlon = rlon2 - rlon1  
    a = sin.(dlat/2).^2 + cos.(rlat1).*cos.(rlat2).*sin.(dlon/2).^2  

    return 2asin.(sqrt.(a))
end # function


"""
    spherearcrad(rlon1,rlat1,rlon2,rlat2)

Returns the great circle (geodesic) distance measured in radians [0,π)
between positions (rlon1,rlat1) and (rlon2,rlat2) [rad].

    c = spherearcrad(rlon1,rlat1,rlon2,rlat2)

Arguments can be vectorized.

Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
Section 4.3.149.
"""
function spherearcrad(rlon1,rlat1,rlon2,rlat2)
    l1 = rlon1
    a  = π/2 .- rlat1

    l2 = rlon2
    c  = π/2 .- rlat2

    x1 = sin.(a).*cos.(l1)         # (x,y,z) of pos 1.
    y1 = sin.(a).*sin.(l1)
    z1 = cos.(a)

    x2 = sin.(c).*cos.(l2)         # (x,y,z) of pos 2.
    y2 = sin.(c).*sin.(l2)
    z2 = cos.(c)
    
    return acos.(x1.*x2+y1.*y2+z1.*z2) # Arc length [rad]
end # function


"""
    spherepos(lon1,lat1,dist,dr)

Finds the point (lon2,lat2) on the sphere separated from point (lon1,lat1) by 
dist [deg], a great circle path which crosses through points (lon1,lat1) and
(lon2,lat2). This great circle path has local direction dr [deg] relative to
north in position (lon1,lat1).

If any element of dist is negative it is taken as the distance in meters on the earth.

    lon2, lat2 = spherepos(lon1,lat1,dist,dr)

Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
Section 4.3.149.

Requires ang180 and ang360.

"""
function spherepos(lon1,lat1,dist,dr)
    if any(dist.<0)
        b = abs(dist)/REARTH
    else
        b = deg2rad.(dist)
    end # if

    a = deg2rad.(90.0 .- lat1)
    cc = deg2rad.(360.0 .- dr)

    c = acos.(cos.(a).*cos.(b)+sin.(a).*sin.(b).*cos.(cc))
    y = sin.(b).*sin.(cc)./(sin.(c)+EPS)
    x = (cos.(b)-cos.(c).*cos.(a))./(sin.(a).*sin.(c)+EPS)
    bb = atan.(y,x)

    lon2 = ang360(lon1-rad2deg.(bb))
    lat2 = ang180(90 .- rad2deg.(c))

    return lon2, lat2
end # function


"""
    spheredir(lon1, lat1, lon2, lat2)

Returns the local direction measured clockwise relative to north from
point (lon1,lat1) towards (lon2,lat2) following a great circle
(geodesic) path on a sphere. Arguments can be vectorized.

# Argument:
 * `lon1`: start longitude [``deg``]
 * `lat1`: start latitude [``deg``]
 * `lon2`: end longitude [``deg``]
 * `lat2`: end latitude [``deg``]

# Returns:
 * `dir`: direction [``deg``]
 
# References:
 *  Abramowitz and Stegun (1970): Handbook of mathematical functions,
    Section 4.3.149.
 
# Requires:
 * spherearcrad
 * ang360
 
# History:
 * 2005-03-14, Matlab, Oyvind.Breivik@met.no.
 * 2017-03-21, Julia v0.5, Oyvind.Breivik@met.no.
 * 2019-07-30, Julia v1.0, Oyvind.Breivik@met.no.


"""
function spheredir(lon1, lat1, lon2, lat2)
    a = deg2rad.(90.0 .- lat1)
    c = deg2rad.(90.0 .- lat2)
    bb = deg2rad.(lon1-lon2)

    b = spherearcrad(deg2rad.(lon1),deg2rad.(lat1),deg2rad.(lon2),deg2rad.(lat2))

    y = sin.(c).*sin.(bb)./(sin.(b)+EPS)
    x = (cos.(c)-cos.(a).*cos.(b))./(sin.(a).*sin.(b)+EPS)

    cc = atan.(y,x)
    return ang360(360.0 .- rad2deg.(cc))
end # function

"""
    rotspher(lon, lat, slon, slat, rot)

Projection for rotated spherical co-ordinates

# Argument:
 * `lon`: longitude [``deg``]
 * `lat`: latitude [``deg``]
 * `slon`: longitude of rotated south pole [``deg``]
 * `slat`: latitude of rotated south pole [``deg``]
 * `rot`: rotation angle [``deg, clockwise from north``]

# Returns:
 * `(rlon, rlat)`: rotated longitude and latitude [``deg``]
 

# References:
 *  Abramowitz and Stegun (1970): Handbook of mathematical functions,
    Section 4.3.149.
 
# Requires:
 * ang360
 
# History:
  * 2005-03-14, Matlab, Oyvind.Breivik@met.no.
  * 2017-03-21, Julia v0.5, Oyvind.Breivik@met.no.

"""
function rotspher(lon, lat, slon, slat, rot)

    a = π/2-deg2rad.(lat)  # co-latitude
    c = π/2+deg2rad.(slat) # co-latitude of rotated north pole
    B = deg2rad.(slon .+ 180 .- lon)

    cosb = cos.(c).*cos.(a)+sin.(c).*sin.(a).*cos.(B) # Saving cosb for precision
    b = acos.(cosb) # [0,π)
    sinA = sin.(a).*sin.(B)./sin.(b)
    cosA = (cos.(a) - cosb.*cos.(c))./(sin.(b).*sin.(c))
    A = atan.(sinA, cosA) # [-π,π)

    rlon = ang360(rad2deg.(A) - rot) # [0,360)
    rlat = rad2deg.(π/2 .- b)      # [-90,90)
    
    return rlon, rlat

end # function

"""
Inverse projection for rotated spherical co-ordinates

# Argument:
 * `rlon`: rotated longitude [``deg``]
 * `rlat`: rotated latitude [``deg``]
 * `slon`: longitude of rotated south pole [``deg``]
 * `slat`: latitude of rotated south pole [``deg``]
 * `rot`: rotation angle [``deg, clockwise from north``]

# Returns:
 * `(lon, lat)`: longitude and latitude [``deg``]
 
# Usage:
 * lon, lat = irotspher(rlon, rlat, slon, slat, rot)
 
# Requires:
 * ang360
 
# References:
  * Abramowitz and Stegun (1970): Handbook of mathematical functions,
    Section 4.3.149.
 
# History:
  2005-03-14, Matlab, Oyvind.Breivik@met.no.
  2017-03-21, Julia v0.5, Oyvind.Breivik@met.no.

"""
function irotspher(rlon, rlat, slon, slat, rot)

    b = π/2-deg2rad.(rlat)  # rotated co-latitude
    c = π/2+deg2rad.(slat)  # co-latitude of rotated north pole
    A = deg2rad.(rot+rlon)

    cosa = cos.(b).*cos.(c)+sin.(b).*sin.(c).*cos.(A) # Saving cosa for precision
    a = acos.(cosa)  # co-latitude [0,π)
    sinB = sin.(b).*sin.(A)./sin.(a)
    cosB = (cos.(b)-cos.(c).*cosa)./(sin.(c).*sin.(a))
    B = atan.(sinB, cosB) # [-π,π)

    lon = ang360(slon .+ 180 .- rad2deg.(B)) # [0,360)
    lat = rad2deg.(π/2 .- a)           # [-90,90)

    return lon, lat

end # function

end # module

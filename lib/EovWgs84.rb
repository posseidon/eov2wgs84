#!/usr/bin/ruby
require 'mathn.rb'
include Math

=begin
  * Name: EovWgs84 Converter
  * Description: Converting EOV Coordinates to GPS coordinates.
  * Author: Binh Nguyen - ELTE Department of Computer Science
  * Date: 2010.09.14
=end

#
# <h1>Converter Class</h1>
# <br>
# <h2>Usage:</h2>
# <br>
# <h3>EovWgs84.convert(LATITUDE, LONGITUDE)</h3>
#
module EovWgs84

  # Struct containing latitude and longitude for Coordinate format.
  Coordinates = Struct.new(:latitude, :longitude)
  
  #
  # Section: Reference Surface.
  #
  @ELIPSODIAL_X = 6378160
  @ELIPSODIAL_Y = 6356774.516
  @SQUARE_ELIPSODIAL_X = @ELIPSODIAL_X*@ELIPSODIAL_X
  @SQUARE_ELIPSODIAL_Y = @ELIPSODIAL_Y*@ELIPSODIAL_Y

  #
  # Section: Displacement
  #
  @DISPLACEMENT_X = 240318.128375
  @DISPLACEMENT_Y = 648196.841798

  #
  # Section: Error Factor
  #
  @ERROR_FACTOR = 1.000719704936

  #
  # Section: Projecting Sphere Radius
  #
  @PROJECTION_RADIUS = 6379296.41898993

  #
  # Section: Elipsodial geographic coordinates
  #
  @LAMBDA_0 = 19.048571778
  @PHI_0 = 47.1 * PI / 180
  
  #
  # Section: Additional
  #
  @ARCSEC = (180*3600) / PI

  @LAMBDA_0_RADIAN = @LAMBDA_0 * PI / 180
  @PHI_O_RADIAN = (47 + 1.0/6) * PI / 180
  @PHI_O_RADIAN_SPHERE = (47+7/60+20.0578/3600)*PI/180
  
  @TAU = (47+1.0/6) * (PI/180)
  @OMEGA = ((@SQUARE_ELIPSODIAL_X - @SQUARE_ELIPSODIAL_Y)/(@SQUARE_ELIPSODIAL_Y))*(cos(@TAU)**2)
  @ROOTS = sqrt(1+@OMEGA)
  @ZET = (1.5*@OMEGA*tan(@TAU))/@ARCSEC
  
  #=begin HD72 - WGS84
    # Bursa-Wolf Type Transformation.
  #=end
  @DX = 52.684 # (dX)
  @DY = -71.194 # (dY)
  @DZ = -13.975 # (dZ)
  @EX = 0.312 # (eX)
  @EY = 0.1063 # (eY)
  @EZ = 0.3729 # (eZ)
  @K = 0.0000010191 # (k)
  @H = 6356752.3142
  @G = 6378137.0
  @F = 6356774.516
  @E = 6378160.0
  @D = 0
  @i = (@E-@F)/@E
  @j = 2*@i-@i*@i
  
  #
  # Converting Coordinates from EOV to HD72.
  #
  def EovWgs84.to_hd72(coord_x, coord_y)
    # Shifting with 200000km and 650000 km
    x = coord_y.to_f - 200000
    y = coord_x.to_f - 650000
    
    fv = 2 * ( atan(E**(x/@PROJECTION_RADIUS)) - PI/4)
    lv = y / @PROJECTION_RADIUS
    
    
    f = asin(sin(fv)*cos(@PHI_0) + sin(@PHI_0) * cos(fv) * cos(lv))
    l = asin(sin(lv)*cos(fv) / cos(f))
    
    
    s = (f-@PHI_O_RADIAN_SPHERE) * @ARCSEC
    
    transform_a = ((0.5 * @OMEGA)/(@ROOTS*@ARCSEC*@ARCSEC))*(-1 + (tan(@TAU)**2) -@OMEGA + (5 * @OMEGA * (tan(@TAU)**2) ))    
    transform_b = (@TAU + (s*@ROOTS)/@ARCSEC - (s*s*@ZET)/@ARCSEC + (s*s*s*transform_a)/@ARCSEC )

    lat = transform_b * (180 / PI)
    lon = (@LAMBDA_0_RADIAN + l/@ERROR_FACTOR)*180/PI
    
    Coordinates.new(lat,lon)
  end
  
  #
  # Converting HD72 coordinates to WGS84.
  #
  def EovWgs84.to_wgs84(coord_x, coord_y)
    b = coord_x.to_f
    c = coord_y.to_f
    k = b * PI/180
    l = c * PI/180
    m = @E/sqrt(1-@j*sin(k)**2)
    n = (m+@D) * cos(k) * cos(l)
    o = (m+@D) * cos(k) * sin(l)
    p = (m * (1-@j)+@D)*sin(k)
    
    x = @DX + (1+@K) * ( n + ((@EZ/3600)*PI/180) * o - ((@EY/3600)*PI/180)* p)
    y = @DY + (1+@K) * (-n * ((@EZ/3600)*PI/180) + o + (p*(@EX/3600)*PI/180) )
    z = @DZ + (1+@K) * ( n * ((@EZ/3600)*PI/180) - o * ((@EX/3600)*PI/180) + p )

    # Latitude
    ah = atan2(y,x)  # koordinatacsere
    lat = ah * 180 / PI
    
    # Longitude
    aa = (@G-@H)/@G
    ab = 2*aa-aa*aa
    ac = (@G*@G-@H*@H)/@H/@H
    ae = sqrt(x*x+y*y)
    af = atan2(z*@G,ae*@H)
    ag = atan2(z + ac*@H*sin(af)**3, ae-ab*@G*cos(af)**3 )
    ad = @G/sqrt(1-ab*sin(ag)**2)
    lon = ag * 180/PI
    
    Coordinates.new(lon,lat)
  end
  
  #
  # Converting EOV to WGS84.
  # Using to_hd72 and to_wgs84
  # Returning Coordinates struct.
  #
  def EovWgs84.convert(coord_x, coord_y)
    hd72 = to_hd72(coord_x, coord_y)
    wgs84 = to_wgs84(hd72[:latitude], hd72[:longitude])
    return wgs84
  end
  
end

 
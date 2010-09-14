require 'test/unit'
require 'EovWgs84.rb'

class TestEov_to_Wgs84 < Test::Unit::TestCase
  
  def test_coord
    gps_coord = EovWgs84.convert(654892.103258,241774.006261)
    latitude = 47.5199036062382
    longitude = 19.1124036720113
    assert_equal(latitude.to_s, gps_coord[:latitude].to_s)
    assert_equal(longitude.to_s, gps_coord[:longitude].to_s)    
  end
  
end

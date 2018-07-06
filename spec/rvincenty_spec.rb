RSpec.describe RVincenty do
  # positions are from Wikipedia
  BAIKONUR_COSMODROME               = [45.965, 63.305         ]
  KENNEDY_SPACE_CENTER              = [28.585, -80.651        ]
  JIUQUAN_SATELLITE_LAUNCH_CENTER   = [40.957778, 100.291667  ]
  SAN_MARCO_PLATFORM                = [-2.938333, 40.2125     ]
  SOUTH_POLE                        = [-90, -0                ]
  NORTH_POLE                        = [90, -0                 ]
  AMUNDSEN_SCOTT_SOUTH_POLE_STATION = [-90, -139.266667       ]

  # distance values are from http://www.movable-type.co.uk/scripts/latlong-vincenty.html
  DISTANCES = [
    [BAIKONUR_COSMODROME, KENNEDY_SPACE_CENTER, 10_986_163.684, 328.551577777778],
    [SOUTH_POLE, JIUQUAN_SATELLITE_LAUNCH_CENTER, 14_537_850.120, 100.291667],
    [SOUTH_POLE, AMUNDSEN_SCOTT_SOUTH_POLE_STATION, 0],
    [SOUTH_POLE, NORTH_POLE, 20_003_931.459],
    [NORTH_POLE, SOUTH_POLE, 20_003_931.459],
    [BAIKONUR_COSMODROME, JIUQUAN_SATELLITE_LAUNCH_CENTER, 3_016_381.278, 87.1952055555555],
    [AMUNDSEN_SCOTT_SOUTH_POLE_STATION, SAN_MARCO_PLATFORM, 9_677_058.827, 179.479167],
    [AMUNDSEN_SCOTT_SOUTH_POLE_STATION, AMUNDSEN_SCOTT_SOUTH_POLE_STATION, 0],
    [JIUQUAN_SATELLITE_LAUNCH_CENTER, BAIKONUR_COSMODROME, 3016381.278, 293.168741],
  ]

  it "has a version number" do
    expect(RVincenty::VERSION).not_to be nil
  end

  describe 'inverse solution' do
    example do
      DISTANCES.each do |a, b, distance|
        calc_distance = RVincenty.distance(a, b)
        expect(distance).to be_within(0.0005).of(calc_distance)
      end
    end

    context 'when readme example' do
      example do
        point_a = 45.965, 63.305 # Baikonur Cosmodrome
        point_b = 28.585, -80.651 # John F. Kennedy Space Center
        distance = RVincenty.distance(point_a, point_b)
        expect(distance).to be_within(0.0005).of(10_986_163.684)
      end
    end
  end

  describe 'direct solution' do
    example do
      DISTANCES.reject do |_, _, _, bearing|
        bearing.nil?
      end.each do |a, (b_lat, b_lng), distance, bearing|
        lat, lng = RVincenty.direct(a, bearing, distance)
        # p "#{lat},#{lng} - #{b_lat}, #{b_lng}"
        expect(lat).to be_within(0.05).of(b_lat)
        expect(lng).to be_within(0.05).of(b_lng)
      end
    end
  end
end

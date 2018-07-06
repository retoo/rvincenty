require "rvincenty/version"

module RVincenty
  # Calculates the distance between two given points.
  # A point is a two-elmeent array containing the points coordinates (Latitude and Longitutde)
  # experessed as a floating point number.
  def self.distance(point_a, point_b)
    lat1, lon1 = point_a
    lat2, lon2 = point_b

    # WGS-84 ellipsiod
    a = 6378137.0
    b = 6356752.3142
    f = 1.0/298.257223563

    l = deg_to_rad(lon2 - lon1)
    u1 = Math.atan((1.0-f) * Math.tan(deg_to_rad(lat1)));
    u2 = Math.atan((1.0-f) * Math.tan(deg_to_rad(lat2)));

    sinU1 = Math.sin(u1)
    cosU1 = Math.cos(u1)
    sinU2 = Math.sin(u2)
    cosU2 = Math.cos(u2)

    lambda = l

    iterLimit = 100
    lambdaP = nil
    while true
      sinLambda = Math.sin(lambda)
      cosLambda = Math.cos(lambda)

      sinSigma = Math.sqrt((cosU2*sinLambda) * (cosU2*sinLambda) +
        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda))

      # co-incident points
      return 0 if (sinSigma==0)

      cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
      sigma = Math.atan2(sinSigma, cosSigma);
      sinalpha = cosU1 * cosU2 * sinLambda / sinSigma;
      cosSqalpha = 1.0 - sinalpha*sinalpha;
      cos2SigmaM = cosSigma - 2.0*sinU1*sinU2/cosSqalpha;
      if cos2SigmaM.nan?
        cos2SigmaM = 0
      end

      c = f/16*cosSqalpha*(4+f*(4-3*cosSqalpha));
      lambdaP = lambda;
      lambda = l + (1.0-c) * f * sinalpha * (sigma + c*sinSigma*(cos2SigmaM+c*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
      iterLimit -= 1

      break if (lambda-lambdaP).abs < 1e-12 || iterLimit <= 0
    end

    #formula failed to converge
    return NaN if (iterLimit==0)

    uSq = cosSqalpha * (a*a - b*b) / (b*b);
    va = 1 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175.0*uSq)))
    vb = uSq/1024.0 * (256.0+uSq*(-128.0+uSq*(74.0-47.0*uSq)))
    deltaSigma = vb*sinSigma*(cos2SigmaM+vb/4.0*(cosSigma*(-1.0+2.0*cos2SigmaM*cos2SigmaM)-
      vb/6*cos2SigmaM*(-3.0+4.0*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)))
    s = b*va*(sigma-deltaSigma)

    return s;
  end

  # Destination given distance & bearing from start point (direct solution)
  def self.direct(point_a, bearing, distance)
    a = 6378137.0
    b = 6356752.31414036
    f = 1 / 298.257222101

    lat1, lng1 = point_a
    az1 = bearing * Math::PI / 180.0
    sin_az1 = Math.sin(az1)
    cos_az1 = Math.cos(az1)

    u1 = Math.atan((1 - f) * Math.tan(lat1 * Math::PI / 180.0))
    sigma1 = Math.atan2((1 - f) * Math.tan(lat1 * Math::PI / 180.0), cos_az1)

    sin_u1 = Math.sin(u1)
    cos_u1 = Math.cos(u1)
    sin_alfa = cos_u1 * sin_az1
    cos_alfa2 = 1 - sin_alfa ** 2

    uu = cos_alfa2 * (a * a - b * b) / (b * b)
    a1 = 1 + (uu / 16384.0) * (4096 + uu * (-768 + uu * (320 - 175 * uu)))
    b1 = (uu / 1024.0) * (256 + uu * (-128 + uu * (74 - 47 * uu)))

    iter_limit = 100
    sigma = distance / (b * a1)
    loop do
      cos2sigmam = Math.cos(2 * sigma1 + sigma)
      sin_sigma = Math.sin(sigma)
      cos_sigma = Math.cos(sigma)
      d_sigma = b1 * sin_sigma * (cos2sigmam + (b1 / 4.0) * (cos_sigma * (-1 + 2 * cos2sigmam ** 2))) - (b1 / 6.0) * cos2sigmam * (-3 + 4 * sin_sigma ** 2) * (-3 + 4 * cos2sigmam ** 2)
      sigma0 = sigma
      sigma = distance / (b * a1) + d_sigma
      iter_limit -= 1
      break if (sigma - sigma0).abs <= 1.0e-12 || iter_limit <= 0
    end

    cos2sigmam = Math.cos(2 * sigma1 + sigma)
    sin_sigma = Math.sin(sigma)
    cos_sigma = Math.cos(sigma)

    phi2 = Math.atan2(sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_az1, (1 - f) * Math.sqrt(sin_alfa * sin_alfa + (sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_az1) ** 2))
    lamuda = Math.atan2(sin_sigma * sin_az1, cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_az1)
    c1 = (f / 16.0) * cos_alfa2 * (4 + f * (4 - 3 * cos_alfa2))
    omega = lamuda - (1 - c1) * f * sin_alfa * (sigma + c1 * sin_sigma * (cos2sigmam + c1 * cos_sigma * (-1 + 2 * cos2sigmam ** 2)))
    lamuda2 = lng1 * Math::PI / 180.0 + omega

    lat = phi2 * 180 / Math::PI
    lng = lamuda2 * 180 / Math::PI
    [lat, lng]
  end

  private

  def self.deg_to_rad(deg)
    (Math::PI * deg) / 180
  end
end
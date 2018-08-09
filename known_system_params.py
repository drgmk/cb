''' Parameters for systems with known circumbinary planets.'''

import numpy as np
import astropy.units as u
import funcs

# Kepler 16
k16 = funcs.CBSystem(m1 = 0.6897, f1 = 1., m2 = 0.20255, f2 = 0.01555,
                     r1 = 0.6489 * u.Rsun.to('au'), r2 = 0.22623 * u.Rsun.to('au'),
                     ab = 0.22431, eb = 0.15944, ib = np.deg2rad(90.3401),
                     wb = np.deg2rad(263.464), fb = np.deg2rad(186.53239),
                     mp = .03e-3, rp = 0.7538 * u.Rjupiter.to('au'),
                     ap = 0.7048, ep = 0.0069, ip = np.deg2rad(90.0322),
                     wp = np.deg2rad(318.0),fp = np.deg2rad(148.92),
                     Wp = np.deg2rad(0.003),
                     t0 = 2455212.12316)

# Kepler 34
k34 = funcs.CBSystem(m1 = 1.0479,f1 = 1.,m2 = 1.0208,f2 = 0.8475,
                     r1 = 1.1618 * u.Rsun.to('au'),r2 = 1.0927 * u.Rsun.to('au'),
                     ab = 0.22882, eb = 0.52087,ib = np.deg2rad(89.8584),
                     wb = 1.2468, fb = 3.4675,
                     mp = 0.220 * u.Mjup.to('Msun'), rp = 0.764 * u.Rjupiter.to('au'),
                     ap=1.0897743291, #ap = 1.0896,
                     ep = 0.182, ip = np.deg2rad(90.355),
                     wp = 0.1378,fp = 2.0623, Wp = np.deg2rad(-1.74),
                     t0 = 2454969.2000)

# Kepler 35
k35 = funcs.CBSystem(m1 = 0.8877,f1 = 1.,m2 = 0.8094,f2 = 0.3941,
                     r1 = 1.0284 * u.Rsun.to('au'),r2 = 0.7861 * u.Rsun.to('au'),
                     ab = 0.17617, eb = 0.1421,ib = np.deg2rad(90.4238),
                     wb = 1.507, fb = 0.06543,
                     mp = 0.127 * u.Mjup.to('Msun'), rp = 0.728 * u.Rjupiter.to('au'),
                     ap = 0.60347, ep = 0.042, ip = np.deg2rad(90.76),
                     wp = 1.1541, fp = 1.3069, Wp = np.deg2rad(-1.24),
                     t0 = 2454965.8500)

k47b = funcs.CBSystem(m1 = 1.043,f1 = 1.,m2 = 0.362,f2 = 0.00568,
                      r1 = 0.964 * u.Rsun.to('au'),r2 = 0.3506 * u.Rsun.to('au'),
                      ab = 0.0836, eb = 0.0234, ib = np.deg2rad(89.34),
                      wb = np.deg2rad(212.3), fb = np.deg2rad(-212.3+90),
                      mp = 0 * u.Mjup.to('Msun'), rp = 2.98 * u.Rearth.to('au'),
                      ap = 0.2956, ep = 0., ip = np.deg2rad(89.59),
                      wp = 0, fp = 2*np.pi*(31.367-29.306346)/49.514+np.pi/2, Wp = np.deg2rad(0.1),
                      t0 = -29.306346+2455000)

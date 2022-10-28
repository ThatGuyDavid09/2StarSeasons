import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import numpy as np

# Change these
planet_orbit_in_au = 1
star_semi_major_axis_solar_radii = 10
star_mass_in_solar_masses = 0.85
sun1_luminosity = 1
sun2_luminosity = 1
axial_tilt = 23.5
planet_mass_in_earth_masses = 1


# All based on this: https://worldbuilding.stackexchange.com/questions/62948/season-cycle-that-would-occur-on-a-habitable-planet-that-orbits-two-suns
# Alright, math time
semi_major_axis_m = star_semi_major_axis_solar_radii * 6.957e+8  # Convert solar radius to meters
planet_orbit_m = planet_orbit_in_au * 1.496e+11  # Convert AU to meters
G = 6.6743e-11
solar_mass = star_mass_in_solar_masses * 1.989e+30
planet_mass = planet_mass_in_earth_masses * 5.972e+24

orbital_period_s = 2 * np.pi * np.sqrt((semi_major_axis_m ** 3) / (G * 2 * solar_mass))
orbital_period = orbital_period_s / 86400
solar_orbit_rad_au = star_semi_major_axis_solar_radii / 2 / 215.032

planet_orbital_period_s = 2 * np.pi * np.sqrt((planet_orbit_m ** 3) / (G * (2 * solar_mass + planet_mass)))
planet_orbital_period = planet_orbital_period_s / 86400 # Days

def circle(radius, pos_rads):
    a = radius * np.cos(pos_rads)
    b = radius * np.sin(pos_rads)
    return a, b


fig, ax = plt.subplots(2, 2)
planet_diagram = ax[0, 0]
dist_diagram = ax[0, 1]
energy_diagram = ax[1, 0]
seasons_diagram = ax[1, 1]

theta = np.linspace(0, 2 * np.pi, 150)

planet_diagram.xaxis.set_ticks_position('bottom')
planet_diagram.yaxis.set_ticks_position('left')
planet_diagram.set_aspect(1)

planet_diagram.plot(*circle(solar_orbit_rad_au, theta), linewidth=1)

# Plot suns
planet_diagram.plot(*circle(solar_orbit_rad_au, 0), marker="o", markersize=15, markerfacecolor="green")
planet_diagram.plot(*circle(solar_orbit_rad_au, np.pi), marker="o", markersize=15, markerfacecolor="green")

# Plot planet and planet orbit
planet_diagram.plot(*circle(planet_orbit_in_au, theta), linewidth=1)
planet_diagram.plot(*circle(planet_orbit_in_au, np.pi / 2), marker="o", markersize=10, markerfacecolor="blue")
planet_diagram.set_title("Solar system diagram")


def sun_1(pos_rads):
    orbit_rad = semi_major_axis_m / 2
    a = orbit_rad * np.cos((planet_orbital_period / orbital_period) * pos_rads)
    b = orbit_rad * np.sin((planet_orbital_period / orbital_period) * pos_rads)
    return np.vstack((a, b))


def sun_2(pos_rads):
    orbit_rad = semi_major_axis_m / 2
    a = orbit_rad * np.cos((planet_orbital_period / orbital_period) * pos_rads)
    b = orbit_rad * np.sin((planet_orbital_period / orbital_period) * pos_rads)
    return np.vstack((-a, -b))


def planet(pos_rads):
    orbit_rad = planet_orbit_in_au * 1.496e+11
    a = orbit_rad * np.cos(pos_rads)
    b = orbit_rad * np.sin(pos_rads)
    return np.vstack((a, b))


def map_num_to_rads(num, max=planet_orbital_period, max_val=2 * np.pi):
    MAX_RADS = max_val
    scaled_val = (num / max) % 10
    return MAX_RADS * scaled_val


x_max = planet_orbital_period
x_lim = x_max
year_scale = np.linspace(0, x_max, 10000)  # 365, 1000)
sun_points_1 = sun_1(map_num_to_rads(year_scale))
plan_points = planet(map_num_to_rads(year_scale))
dist_to_sun1 = np.hypot(*(sun_points_1 - plan_points))
dist_diagram.plot(year_scale, dist_to_sun1)

sun_points = sun_2(map_num_to_rads(year_scale))
# plan_points = planet(map_num_to_rads(year_scale))
dist_to_sun2 = np.hypot(*(sun_points - plan_points))
dist_diagram.set_title("Distance to stars")
dist_diagram.plot(year_scale, dist_to_sun2)
dist_diagram.set_xlim([0, x_lim])


# Solar energy, need to account for star eclipsing
# Formula from https://worldbuilding.stackexchange.com/questions/87137/how-to-calculate-the-light-received-by-a-planet-during-a-binary-star-eclipse
def solar_energy(dist_to_sun1, dist_to_sun2):
    ang_diff = (dist_to_sun1 / dist_to_sun2) - 1
    total_b = (sun1_luminosity / (dist_to_sun1 ** 2)) + ang_diff * (sun2_luminosity / (dist_to_sun2 ** 2))
    return total_b


solar_vals = solar_energy(dist_to_sun1, dist_to_sun2)
energy_diagram.plot(year_scale, solar_vals / np.linalg.norm(solar_vals), color="black")
energy_diagram.set_title("Total energy received")
energy_diagram.set_xlim([0, x_lim])


# Assume that planet is at max axial tilt at 0 radians
def northern_energy(energy_recv, pos_num, north_degrees=45):
    pos_rads = map_num_to_rads(pos_num)
    mapping_1 = interp1d([0, np.pi], [np.radians(-axial_tilt), np.radians(axial_tilt)])
    mapping_2 = interp1d([0, np.pi], [np.radians(axial_tilt), np.radians(-axial_tilt)])
    rad_positions = np.where(pos_rads < np.pi, mapping_1(pos_rads % np.pi), mapping_2(pos_rads % np.pi))
    modifiers = np.cos(np.radians(north_degrees) + rad_positions)
    return modifiers * energy_recv


energy_vals = northern_energy(solar_vals, year_scale)
seasons_diagram.plot(year_scale, energy_vals, color="black")
seasons_diagram.set_title("Energy received at 45N")
# seasons_diagram.set_xlim([0, x_lim])

# show the plot
plt.show()

from astroquery.vizier import Vizier
import requests

# Create a custom session to disable SSL verification
session = requests.Session()
session.verify = False

# Attach the custom session to Vizier
Vizier._session = session

# Your existing code
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
from geopy.geocoders import Nominatim

# Parameters
num_stars = 700
canvas_width_cm = 120
canvas_height_cm = 180
output_pdf = "starry_sky_map.pdf"
city = "Donaueschingen, Germany"
date_time = "2002-07-28 12:00:00"

# Fetch location coordinates
geolocator = Nominatim(user_agent="star_map_creator")
location = geolocator.geocode(city)
latitude, longitude = location.latitude, location.longitude
print(f"Using coordinates: Latitude={latitude}, Longitude={longitude}")

# Convert location and time to astropy objects
observer_location = EarthLocation(lat=latitude, lon=longitude)
observation_time = Time(date_time)

# Fetch star data
vizier = Vizier(columns=["RAJ2000", "DEJ2000", "Vmag"])
result = vizier.query_constraints(catalog="I/239/hip_main", Vmag="<10")  # Bright stars
stars = result[0].to_pandas()

# Convert star coordinates to AltAz (observer's perspective)
star_coords = SkyCoord(ra=stars["RAJ2000"].values, dec=stars["DEJ2000"].values, unit="deg")
altaz_frame = AltAz(obstime=observation_time, location=observer_location)
altaz = star_coords.transform_to(altaz_frame)

# Filter stars above the horizon
stars = stars[altaz.alt.deg > 0]
alt = altaz.alt.deg
az = altaz.az.deg

# Select top N stars by brightness (Vmag)
stars = stars.iloc[np.argsort(stars["Vmag"])].head(num_stars)
alt = alt[np.argsort(stars["Vmag"])][:num_stars]
az = az[np.argsort(stars["Vmag"])][:num_stars]

# Normalize positions to fit canvas
az_norm = (az - az.min()) / (az.max() - az.min()) * canvas_width_cm
alt_norm = (alt - alt.min()) / (alt.max() - alt.min()) * canvas_height_cm

# Plot the star map
plt.figure(figsize=(canvas_width_cm / 2.54, canvas_height_cm / 2.54))
plt.scatter(az_norm, alt_norm, s=10, c="white")
plt.gca().set_facecolor("black")
plt.axis("off")

# Add title
plt.title(f"Star Map for {city} on {date_time}", color="white", fontsize=12, pad=20)

# Save as PDF
plt.savefig(output_pdf, format="pdf", dpi=300, bbox_inches="tight")
plt.close()

print(f"Starry sky map saved to {output_pdf}")

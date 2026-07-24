# %% code to extract some data from the HYCOM zarr files and export to csv
import xarray as xr
import numpy as np

# %% open and extract ssh and velocity data stores
dirin = '/Volumes/Salmacis/hycom/zstore/'
filein = dirin + 'hycom12-ssh-1-rechunked-corr.zarr'
ds_ssh = xr.open_zarr(filein).isel(time=[0,6,12,18])

filein_uv = dirin + 'hycom12-1-rechunked-corr.zarr'
ds_uv = xr.open_zarr(filein_uv).isel(time=[0, 6, 12, 18], Depth=0)

# %% South of Gulf Stream: 22-32N, 69-59W
dxy = 10
y1 = 22
y2 = y1+dxy-1
x1 = -69
x2 = x1+dxy

from clouddrift.sphere import distance
print(distance(0.5*(x1+x2),y1,0.5*(x1+x2),y2)/1e3)
print(distance(x1,0.5*(y1+y2),x2,0.5*(y1+y2))/1e3)

# %%
mask = (
    (ds_ssh["Latitude"] >= y1)
    & (ds_ssh["Latitude"] <= y2)
    & (ds_ssh["Longitude"] >= x1)
    & (ds_ssh["Longitude"] <= x2)
).compute()

# %% 281 by 251? now 251 by 251
ssh = ds_ssh["ssh"].where(mask, drop=True)
u = ds_uv["u"].where(mask, drop=True)
v = ds_uv["v"].where(mask, drop=True)
x = ssh["Longitude"].to_numpy()
y = ssh["Latitude"].to_numpy()

# %% flatten data

x_flat = np.ravel(x, order="F")
y_flat = np.ravel(y, order="F")

# %% create a list of flatten array for time steps 1 to 4

ssh_flats = []
u_flats = []
v_flats = []

for i in range(4):
    ssh_flat = np.ravel(ssh.isel(time=i).values, order="F")
    u_flat = np.ravel(u.isel(time=i).values, order="F")
    v_flat = np.ravel(v.isel(time=i).values, order="F")
    ssh_flats.append(ssh_flat)
    u_flats.append(u_flat)
    v_flats.append(v_flat)

# %% write to different csv files for each time step
# x, y, ssh 
# x, y, u
# x, y, v 
for i in range(4):
    np.savetxt(f"csv/ssh_time_{i}.csv", np.column_stack((x_flat, y_flat, ssh_flats[i])), delimiter=",", header="Longitude,Latitude,SSH", comments="")
    np.savetxt(f"csv/u_time_{i}.csv", np.column_stack((x_flat, y_flat, u_flats[i])), delimiter=",", header="Longitude,Latitude,U", comments="")
    np.savetxt(f"csv/v_time_{i}.csv", np.column_stack((x_flat, y_flat, v_flats[i])), delimiter=",", header="Longitude,Latitude,V", comments="")


# %% next, use pyproj to convert to distance on a plane, then export to csv
# this from ChatGPT ...
from pyproj import CRS, Transformer

lon0 = np.mean(x[0,:])
lat0 = np.mean(y[:,0])

crs_wgs84 = CRS.from_epsg(4326)  # WGS84

crs_local = CRS.from_proj4(f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=km +no_defs")

transformer = Transformer.from_crs(crs_wgs84, crs_local,always_xy=True)

x_flat_km, y_flat_km = transformer.transform(x_flat, y_flat)

# %% is this working?
import matplotlib.pyplot as plt
%matplotlib widget

plt.figure()
plt.plot(x_flat_km,y_flat_km,'.')
# plot the first and last points to see if they are in the right place
plt.plot(x_flat_km[0],y_flat_km[0],'ro')
plt.plot(x_flat_km[250],y_flat_km[250],'co')
plt.plot(x_flat_km[-1],y_flat_km[-1],'go')
plt.xlabel('X (km)')
plt.ylabel('Y (km)')
plt.title('Transformed Coordinates')
plt.show()

# %% now write the same csv files but with x_flat_km and y_flat_km instead of lon/lat
for i in range(4):
    np.savetxt(f"csv/ssh_time_{i}_km.csv", np.column_stack((x_flat_km, y_flat_km, ssh_flats[i])), delimiter=",", header="X_km,Y_km,SSH", comments="")
    np.savetxt(f"csv/u_time_{i}_km.csv", np.column_stack((x_flat_km, y_flat_km, u_flats[i])), delimiter=",", header="X_km,Y_km,U", comments="")
    np.savetxt(f"csv/v_time_{i}_km.csv", np.column_stack((x_flat_km, y_flat_km, v_flats[i])), delimiter=",", header="X_km,Y_km,V", comments="")


# %%

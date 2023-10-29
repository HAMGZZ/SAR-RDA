# SAR RDA
This is a simple simulation of SAR and RDA in matlab (can be easily changed to python with numpy).

This is almost the simplest implementation of the Range Doppler Algorithm out there.

`echo.mat` is a 4096 x4096 data matrix which contains real `echos` from RADARSAT-1 and should result in an image of Vancouver bay.

One can toggle `simulation` to either decode real data, or generated point targets. `NumberofSimTargets` determines the number of point targets in simulation.

`fast` toggles whether imshow or pcolor is used to display data. pcolor is arguably more beautiful, but greatly increases the processing time.

Please let me know if you are able to produce high quality images (one can make out roads or buildings) from the data.


Please give me credit if this is usefull!


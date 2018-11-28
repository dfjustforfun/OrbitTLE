##
# @file calculateOrbitTLE.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-11-28


# from datetime import datetime, timezone
import os
import argparse
import datetime
from pythonOrbitTools.Core.Tle import Tle, Field
from pythonOrbitTools.Orbit.Satellite import Satellite
from pythonOrbitTools.Core.Site import Site

# ROOT_DIR = "/home/cpf/Documents/OrbitTLE/"
ROOT_DIR = os.path.split(os.path.realpath(__file__))[0]

def readTLEData(filePath):
    with open(filePath, "r") as f:
        content = f.read()
        content = content.strip()
        lines = content.split('\n')
    return lines[0].strip(), lines[1].strip(), lines[2].strip()


parser = argparse.ArgumentParser(description= "calculateOrbitTLE.py")
parser.add_argument('--lon',        required=True, type=float,  help="Longitude, in degrees.")
parser.add_argument('--lat',        required=True, type=float,  help="Latitude, in degrees.")
parser.add_argument('--alt',        required=True, type=float,  help="The alitude of the site above the ellipsoid model, in kilometers.")
parser.add_argument('--startTime',  required=True, type=str,    help="YYYY-mm-dd/HH:MM:SS e.g. 2018-08-21/17:23:34")
parser.add_argument('--endTime',    required=True, type=str,    help="YYYY-mm-dd/HH:MM:SS e.g. 2018-08-21/18:04:05")
args = parser.parse_args()

lon         = args.lon
lat         = args.lat
alt         = args.alt
startTime   = datetime.datetime.strptime(args.startTime, "%Y-%m-%d/%H:%M:%S")
endTime     = datetime.datetime.strptime(args.endTime,   "%Y-%m-%d/%H:%M:%S")

# convert to UTC time
startTime_UTC   = startTime + datetime.timedelta(hours=-8.0)
endTime_UTC     = endTime   + datetime.timedelta(hours=-8.0)

satName, line1, line2 = readTLEData(os.path.join(ROOT_DIR, "tle.txt"))

tle = Tle(line1, line2, satName)
satellite = Satellite(tle)
siteEuqator = Site().InitializeByDegLatAndDegLonAndKmAltAndName(lat, lon, alt)

with open(os.path.join(ROOT_DIR, "output.eph"), "w") as f:
    dataTime_UTC = startTime_UTC

    while dataTime_UTC <= endTime_UTC:

        topoLook = siteEuqator.GetLookAngle(satellite.PositionEciByDateTime(dataTime_UTC))

        if topoLook.ElevationDeg >= 3.0:
            dataTime = dataTime_UTC + datetime.timedelta(hours=8.0)
            string = "{} {} {}\n".format(dataTime.strftime("%Y/%m/%d %H:%M:%S"), round(topoLook.ElevationDeg, 3), round(topoLook.AzimuthDeg, 3))
            f.write(string)

        dataTime_UTC = dataTime_UTC + datetime.timedelta(seconds=1.0)

print("Done!")

from matplotlib import pyplot as plt
from astropy import coordinates as coord
import astropy.units as u
import numpy as np
from urllib import parse
from urllib.request import urlopen
from urllib.error import HTTPError
from bs4 import BeautifulSoup as bs


class Flux:
    def __init__(self, resolvedFlux, unresolvedFlux):
        """Initialises a flux (contains unresolved and resolved flux)
        Inputs
        ------
        - resolvedFlux: the resolved flux of the source 
        - unresolvedFlux: the unresolved flux of the source 
        """
        self.resolved = resolvedFlux
        self.unresolved = unresolvedFlux
        assert isinstance(unresolvedFlux, float)
        assert isinstance(resolvedFlux, float)
    
class Source:

    def __init__(self, name, ivsname, coord, noObs, fluxes, calbool):
        """Initializes a source.

        Inputs
        ------
        - name: the J2000 name
        - ivsname: the ivsname
        - coord: astropy coord 
        - noObs: the number of times this source has been observed in the IVS cat (useful filter maybe)
        - fluxes: dict of fluxes with band as key names.
        """
        self.name = name
        self.ivsname = ivsname
        self.coord = coord
        self._noObs = noObs
        self.flux = fluxes
        self.isCal = calbool
        
    def plot_elevation(self, stations, obsTimes):
        f = plt.figure(1, figsize=(10,5))
        ax = f.add_subplot(111)
        ax.set_title("Elevation vs time for {}".format(self.name))
        ax.set_ylabel("Elevation (deg)")
        ax.set_xlabel("Time")
        for station in stations:
            els = station.source_elevation(self.coord, obsTimes).deg
            ax.plot(obsTimes.datetime, els, label=station.station)
    
        plt.legend()
        plt.show()

    def find_nmes(self):
        #scrape the ftp tests
        ftpPage = urlopen("http://www.evlbi.org/tog/ftp_fringes/ftp.html")
        bsFtpPage = bs(ftpPage.read(), features="html5lib")
        links=[]
        for link in bsFtpPage.findAll('a'):
            href = link.get('href')
            if "ftp_fringes" in href:
                try:
                    ftPage = urlopen(href)
                    ftPageTxt = str(ftPage.read())
                    if self.name in ftPageTxt:
                        links.append(link.get('href'))
                    elif self.ivsname in ftPageTxt:
                        links.append(link.get('href'))
                except HTTPError as err:
                    pass
        return links

    def get_astrogeo_link(self):
        #get a link for the astrogeo html section for this source (contains maps/uvrad etc)
        sourceCoordString = parse.quote("ra={:02.0f}:{:02.0f}:{:06.3f}&dec={:+03.0f}:{:02.0f}:{:06.3f}&num_sou=1&format=html".format(*self.coord.ra.hms, *self.coord.dec.dms), safe='=&')
        return("http://astrogeo.org/cgi-bin/calib_search_form.csh?{}".format(sourceCoordString))
        
    
# def load_rfc_cat(filename):
#     #loads the data initially from the rfc catalogue, note we have to take care of < in the fluxes. We ignore position error columns
#     sourceCat = np.loadtxt(filename,
#                            comments='#',
#                            dtype={'names': ('cal', 'ivsname', 'name', 'raH', 'raM', 'raS', 'decD', 'decM',
#                                             'decS', 'noObs', 'fluxSR', 'fluxSU', 'fluxCR', 'fluxCU',
#                                             'fluxXR', 'fluxXU', 'fluxUR', 'fluxUU', 'fluxKR', 'fluxKU'),
#                                   'formats': ('bool', '|S15', '|S15', 'int', 'int', 'float', 'int', 'int',
#                                               'float', 'int', 'float', 'float', 'float', 'float', 'float',
#                                               'float', 'float', 'float', 'float', 'float')},
#                            usecols=(0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19,20,21,22),
#                            converters={0: lambda cal: True if cal.decode()=='C' else False,
#                                        13: lambda f: 0.0 if '<' in str(f) else f,
#                                        14: lambda f: 0.0 if '<' in str(f) else f,
#                                        15: lambda f: 0.0 if '<' in str(f) else f,
#                                        16: lambda f: 0.0 if '<' in str(f) else f,
#                                        17: lambda f: 0.0 if '<' in str(f) else f,
#                                        18: lambda f: 0.0 if '<' in str(f) else f,
#                                        19: lambda f: 0.0 if '<' in str(f) else f,
#                                        20: lambda f: 0.0 if '<' in str(f) else f,
#                                        21: lambda f: 0.0 if '<' in str(f) else f,
#                                        22: lambda f: 0.0 if '<' in str(f) else f})

#     return sourceCat

def load_rfc_cat(filename, minFluxBand='c', minFlux=1.0):
    with open(filename, 'rt') as fin:
        coordStrings = []
        sources = []
        coord0 = coord.SkyCoord("00h00m00.0s +00d00m00.0s")
        for line in fin:
            if not (line.startswith('#') or line.startswith('N') or line.startswith('U')):
                cols = line.split()
                cols[13:23] = [float(f)  if '<' not in f else 0.0 for f in cols[13:23]]
                fluxes = {'s': Flux(cols[13], cols[14]),
                          'c': Flux(cols[15], cols[16]),
                          'x': Flux(cols[17], cols[18]),
                          'u': Flux(cols[19], cols[20]),
                          'k': Flux(cols[21], cols[22])}
                if fluxes[minFluxBand].unresolved > minFlux:
                    name = cols[2]
                    ivsname = cols[1]
                    #coords = coord.SkyCoord("{}h{}m{}s {}d{}m{}s".format(*cols[3:9]))
                    coordStrings.append("{}h{}m{}s {}d{}m{}s".format(*cols[3:9]))
                    sources.append(Source(name, ivsname, coord0, int(cols[12]), fluxes, True))
    coords = coord.SkyCoord(np.array(coordStrings))
    
    for c,s in zip(coords, sources):
        s.coord = c
    return sources



def get_up_sources(stationList, sourceList, obsTimes, minEl=20, minFlux=0.5, minFluxBand='c'):
    sources=[]
    for source in sourceList:
        for station in stationList:
            if not station.is_source_visible(source.coord, obsTimes, minEl*u.deg):
                break
        else:
            sources.append(source)

    
    sources.sort(key=lambda source:source.flux[minFluxBand].unresolved, reverse=True)
    return sources


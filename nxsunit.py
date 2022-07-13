# This program is public domain
# -*- coding: UTF-8 -*-
# Author: Paul Kienzle
"""
Define unit conversion support for NeXus style units.

The unit format is somewhat complicated.  There are variant spellings
and incorrect capitalization to worry about, as well as forms such as
"mili*metre" and "1e-7 seconds".

This is a minimal implementation.  It does not support the complete 
dimensional analysis provided by the package UDUnits on which NeXus is 
based, or even all the units used in the NeXus definition files.

Unlike other units modules, this module does not carry the units along
with the value, but merely provides a conversion function for
transforming values.

Usage example::

    import nxsunit
    # Simple call to convert a value from millimeters to meters
    v = nxsunit.convert(3000, 'mili*metre', 'm')

    # A reusable conversion function from hours to milliseconds
    u = nxsunit.converter('hours', 'ms')
    print(u(1))
    3600000

    # The underlying conversion object if you know the stored units
    # but don't yet know the target units.
    u = nxsunit.Converter('mili*metre')
    print(u.convert(3000, 'm')) # convert to meters
    3

NeXus example::

    # Load sample orientation in radians regardless of how it is stored.
    file = h5py.File('path/to/file.nxs')
    entry = file['/entry1']
    orientation = nxsunit.data_as(entry, 'sample/sample_orientation', 'radians')

The converter knows the dimension it is working with.  For example, if
the dimension is time, then u(300,'hour') will return the time in hours,
but if the dimension is angle, then u(300,'hour') will raise an error.
Using minutes will work in both cases.  When constructing a converter,
you may need to specify a dimension.  For example, Converter('second')
will create a time converter, but Converter('second', 'angle') will create
an angle converter.

The list of ambiguities, and the default dimension is given in the 
unit.AMBIGUITIES map.  The available dimensions and the conversion factors
are given in the unit.DIMENSIONS map.  Note that the temperature converters
have a scale and an offset rather than just a scale.

In addition to converting values :function:`data_as` provides some handy
features for loading nexus fields: (1) byte values are converted to unicode
if you give units as 'str'; (2) default values are filled in if the field
is missing; (3) data can be converted to the desired numpy type (single,
double, etc.); (4) scalar fields can be returned as scalars even if they are
stored as length one vectors; (5) scalar fields can be extended to the number
of points in the scan.

The nexus reader *data_as()* depends on h5py and numpy. If you are not using
this library to load HDF data then they are not required.
"""

# TODO: Parse the udunits database directly
# UDUnits:
#  http://www.unidata.ucar.edu/software/udunits/udunits-1/udunits.txt

from __future__ import division
import math

__all__ = ['Converter']
__version__ = "1.0"

# Limited form of units for returning objects of a specific type.
# Maybe want to do full units handling with e.g., pyre's
# unit class. For now lets keep it simple.  Note that
def _build_metric_units(unit, abbr):
    """
    Construct standard SI names for the given unit.
    Builds e.g.,
        s, ns
        second, nanosecond, nano*second
        seconds, nanoseconds
    Includes prefixes for femto through peta.

    Ack! Allows, e.g., Coulomb and coulomb even though Coulomb is not
    a unit because some NeXus files store it that way!

    Returns a dictionary of names and scales.
    """
    prefix = dict(
        peta=1e15, tera=1e12, giga=1e9, mega=1e6, kilo=1e3,
        deci=1e-1, centi=1e-2, milli=1e-3, mili=1e-3, micro=1e-6,
        nano=1e-9, pico=1e-12, femto=1e-15)
    short_prefix = dict(
        P=1e15, T=1e12, G=1e9, M=1e6, k=1e3,
        d=1e-1, c=1e-2, m=1e-3, u=1e-6,
        n=1e-9, p=1e-12, f=1e-15)
    # \u1D6DC is mathematical bold small mu
    # \u1D707 is mathematical italic small mu
    # \u1D741 is mathematical bold italic small mu
    # \u1D77B is mathematical sans-serif bold small mu
    # \u1D7B5 is mathematical sans-serif bold italic small mu
    short_prefix['\u00B5'] = 1e-6 # MICRO SIGN
    short_prefix['\u03BC'] = 1e-6 # Greek lowercase mu

    map = {abbr: 1}
    map.update([(P+abbr, scale) for (P, scale) in short_prefix.items()])
    for name in [unit, unit.capitalize()]:
        map.update({name: 1, name+'s': 1})
        map.update([(P+name, scale) for (P, scale) in prefix.items()])
        map.update([(P+'*'+name, scale) for (P, scale) in prefix.items()])
        map.update([(P+name+'s', scale) for (P, scale) in prefix.items()])
    return map

def _build_plural_units(**kw):
    """
    Construct names for the given units.  Builds singular and plural form.
    """
    map = {}
    map.update([(name, scale) for name, scale in kw.items()])
    map.update([(name+'s', scale) for name, scale in kw.items()])
    return map

def _build_degree_units(name, symbol, conversion):
    """
    Builds variations on the temperature unit name, including the degree
    symbol or the word degree.
    """
    map = {}
    map[symbol] = conversion
    for s in symbol, symbol.lower():
        map['deg'+s] = conversion
        map['deg_'+s] = conversion
        map['deg '+s] = conversion
        map['°'+s] = conversion
    for s in name, name.capitalize(), symbol, symbol.lower():
        map[s] = conversion
        map['degree_'+s] = conversion
        map['degree '+s] = conversion
        map['degrees '+s] = conversion
    return map

def _build_inv_units(names, conversion):
    """
    Builds variations on inverse units, including 1/x, invx and x^-1.
    """
    map = {}
    for s in names:
        map['1/'+s] = conversion
        map['inv'+s] = conversion
        map[s+'^-1'] = conversion
        map[s+'^{-1}'] = conversion
    return map

def _build_inv2_units(names, conversion):
    """
    Builds variations on inverse square units, including 1/x^2, invx^-2 and x^-2.
    """
    map = {}
    for s in names:
        map['1/'+s+'^2'] = conversion
        map['inv'+s+'^2'] = conversion
        map[s+'^-2'] = conversion
        map[s+'^{-2}'] = conversion
    return map

def _caret_optional(s):
    """
    Strip '^' from unit names.
    * WARNING * this will incorrectly transform 10^3 to 103.
    """
    stripped = [(k.replace('^',''),v) for k, v in s.items() if '^' in k]
    s.update(stripped)

def _build_all_units():
    """
    Fill in the global variables DIMENSIONS and AMBIGUITIES for all available
    dimensions.
    """
    # Gather all the ambiguities in one spot
    AMBIGUITIES['A'] = 'distance'     # distance (angstroms), current (amperes)
    AMBIGUITIES['second'] = 'time'    # time, angle
    AMBIGUITIES['seconds'] = 'time'
    AMBIGUITIES['sec'] = 'time'
    AMBIGUITIES['°'] = 'angle'        # temperature, angle
    AMBIGUITIES['minute'] = 'angle'   # time, angle
    AMBIGUITIES['minutes'] = 'angle'
    AMBIGUITIES['min'] = 'angle'
    AMBIGUITIES['C'] = 'charge'       # temperature (celsius), charge (coulombs)
    AMBIGUITIES['F'] = 'temperature'  # temperature
    AMBIGUITIES['R'] = 'radiation'    # temperature:rankines, radiation:roentgens
    AMBIGUITIES['rad'] = 'angle'      # angle, radiation
    AMBIGUITIES['mRad'] = 'angle'     # angle, radiation

    # Various distance measures
    distance = _build_metric_units('meter', 'm')
    distance.update(_build_metric_units('metre', 'm'))
    distance.update(_build_plural_units(micron=1e-6))
    # Unicode has a couple of different code points for Angstroms:
    #  \u00C5  Latin capital letter A with ring above (preferred)
    #  \u212B  Angstrom sign
    angstroms = (
        'A', # conflicts with amperes
        '\u00C5', '\u212B', 'Ang', 'ang',
        'angstrom', 'angstroms', 'ångström', 'ångströms',
        'Angstrom', 'Angstroms', 'Angstroem', 'Angstroems',
    )
    distance.update({s: 1e-10 for s in angstroms})
    DIMENSIONS['distance'] = distance

    # Various time measures.
    # Note: minutes are used for angle rather than time
    # Note: months/years are varying length so can't be used without date support
    time = _build_metric_units('second', 's')
    time.update(_build_plural_units(
        minute=60, hour=3600, day=24*3600, week=7*24*3600))
    time.update({'sec': 1, 'min': 60, 'hr': 3600})
    time.update({'1e-7 s': 1e-7, '1e-7 second': 1e-7, '1e-7 seconds': 1e-7})
    DIMENSIONS['time'] = time

    # Various angle measures.
    # Note: seconds are used for time rather than angle
    angle = _build_plural_units(
        degree=1, minute=1/60., second=1/3600.,
        arcdegree=1, arcminute=1/60., arcsecond=1/3600.,
        radian=180/math.pi,
        )
    angle.update(
        deg=1, min=1/60., sec=1/3600.,
        arcdeg=1, arcmin=1/60., arcsec=1/3600., 
        angular_degree=1, angular_minute=1/60., angular_second=1/3600., 
        rad=180/math.pi, mRad=0.180/math.pi, mrad=0.180/math.pi,
        )
    angle['\u00B0'] = 1  # unicode degree sign
    DIMENSIONS['angle'] = angle

    frequency = _build_metric_units('hertz', 'Hz')
    frequency.update(_build_inv_units(('s', 'sec', 'second'), 1))
    frequency.update(
        ('/'.join((c, s)), 1) 
        for c in ('count', 'counts')
        for s in ('s', 'sec', 'second'))
    frequency.update(cps=1)
    frequency.update(_build_plural_units(RPM=1/60., rpm=1/60., cpm=1/60.))
    frequency['rev/min'] = '1/60.'
    DIMENSIONS['frequency'] = frequency

    # Note: degrees are used for angle
    temperature = _build_metric_units('kelvin', 'K')
    for k, v in temperature.items():
        temperature[k] = (v, 0)  # add offset 0 to all kelvin temperatures
    temperature.update(
        _build_degree_units('celcius', 'C', (1, 273.15)))
    temperature.update(
        _build_degree_units('centigrade', 'C', temperature['degC']))
    temperature.update(
        _build_degree_units('fahrenheit', 'F', (5./9., 491.67-32)))
    temperature.update(
        _build_degree_units('rankine', 'R', (5./9., 0)))
    # special unicode symbols for fahrenheit and celcius
    temperature['\u2103'] = temperature['degC'] # unicode degrees C
    temperature['\u2109'] = temperature['degF'] # unicode degrees F
    DIMENSIONS['temperature'] = temperature

    pressure = _build_metric_units('pascal','Pa')
    pressure.update(
            bar=1e5, bars=1e5, mbar=1e2, 
            atm=101325, 
            mmHg=133.322387415,
            Torr=101325./760.,
            mTorr=101.325/760.,
            )
    DIMENSIONS['pressure'] = pressure

    charge = _build_metric_units('coulomb', 'C')
    charge.update({'microAmp*hour': 0.0036})
    DIMENSIONS['charge'] = charge

    resistance = _build_metric_units('ohm', 'Ω')
    DIMENSIONS['resistance'] = resistance

    sld = _build_inv2_units(angstroms, 1)
    sld.update(_build_inv2_units(('nm',), 100))
    sld['10^-6 Angstrom^-2'] = 1e-6
    sld['10^-6/Ang^2'] = 1e-6
    sld['10^-6/ang^2'] = 1e-6
    sld['10^-6/A^2'] = 1e-6
    sld['1E-6/Ang^2'] = 1e-6
    sld['1E-6/A^2'] = 1e-6
    _caret_optional(sld)
    DIMENSIONS['sld'] = sld

    Q = _build_inv_units(angstroms, 1)
    Q.update(_build_inv_units(('nm',), 10))
    Q.update(_build_inv_units(('cm',), 1e8))
    Q.update(_build_inv_units(('m',), 1e10))
    Q['10^-3 Angstrom^-1'] = 1e-3
    _caret_optional(Q)
    DIMENSIONS['Q'] = Q

    energy = _build_metric_units('electronvolt', 'eV')
    DIMENSIONS['energy'] = energy
    # Note: energy <=> wavelength <=> velocity requires a probe type

    count = _build_plural_units(count=1, step=1)
    DIMENSIONS['count'] = count

    # APS files may be using 'a.u.' for 'arbitrary units'.  Other
    # facilities are leaving the units blank, using ??? or not even
    # writing the units attributes.
    dimensionless = {None: 1, '': 1}
    dimensionless.update({'???': 1, 'a.u.': 1})
    DIMENSIONS['dimensionless'] = dimensionless

# Initialize DIMENSIONS and AMBIGUITIES
DIMENSIONS = {}
AMBIGUITIES = {}
_build_all_units()

def convert(value, source, target, dimension=None):
    """
    Return *value* converted from *source* units to *target* units.
    """
    return Converter(source, dimension).convert(value, target)

def converter(source, target, dimension=None):
    """
    Return function to convert value from *source* units to *target* units.
    """
    return Converter(source, dimension).conversion(target)

def data_as(group, fieldname, units, rep=None, NA=None, dtype=None, dimension=None):
    """
    Return value of hdf field in the desired units. This is a helper function
    for grabbing data from the NeXus file of the correct shape and units.

    *group* is the HDF file or entry loaded from h5py.

    *fieldname* is the relative path from group to field.

    *units* is the target units. Include dimension if units might be ambiguous.
    Use *units='str'* to convert bytes to unicode.

    *rep* is the number of points in the scan. Single values are extended to
    length *rep*. If *rep* is None then single values are returned as single
    values instead of 0-D arrays or length one vectors.

    *NA* is the default field value if the field is not available or None if
    missing fields should return None. If *rep* is not None then the default
    value will automatically be extended to length *rep*.

    *dtype* is the desired return type as any numpy dtype or dtype name.

    *dimension* is used if the source units might be ambiguous and the
    target units are not the dimension defined in *nxsunit.AMBIGUITIES*.
    """
    import numpy as np
    if fieldname not in group:
        if NA is not None and rep is not None:
            return np.repeat(NA, rep, axis=0)
        return NA
    # Get the value
    field = group[fieldname]
    value = field[()]
    # Convert the value.
    if units == 'str':
        # Bytes to unicode.
        utf8 = np.vectorize(lambda el: el.decode('utf-8'))
        value = utf8(value)
    else:
        # Unit conversion.
        source_units = field.attrs.get(b'units', b'')
        if isinstance(source_units, bytes):
            source_units = source_units.decode('utf-8')
        if source_units:
            value = convert(value, source_units, units, dimension=dimension)
        # Numpy type conversion.
        if dtype is not None:
            value = np.asarray(value, dtype=dtype)
    # TODO: How do we tell it to leave the dimensions alone?
    # Correct the dimensions.
    if rep is None:
        # Squeeze first dimension if rep is None.
        if value.ndim == 0:
            return value[()]
        elif value.shape[0] == 1:
            return value[0]
        else:
            return value
    elif value.ndim >= 1 and len(value) == rep:
        # Do nothing if length matches rep.
        return value
    elif value.ndim == 0 or value.shape[0] == 1:
        # Repeat value if length does not match rep.
        return np.repeat(value, rep, axis=0)
    else:
        # Raise error if length does not match rep.
        raise ValueError(
            f"Length of {field.name} does not match number of scan points"
            f" in {field.file.filename}")

class Converter(object):
    """
    Unit converter for NeXus style units.

    *units* gives the source units for the conversion.

    *dimension* (distance, time, angle, current, etc.) can be used
    for ambiguous units, such as *Converter('A', 'current')* to force amperes
    rather than angstroms. See *DIMENSIONS.keys()* for the list of dimensions
    and *AMBIGUITIES* for the default dimension for each ambiguous unit.

    The returned object is called with the source value and the target
    units. Use `converter.conversion(target)` to get a conversion function
    that can be applied to all source values.

    Raises TypeError if units or dimensions are unknown or inconsistent.
    """
    def __init__(self, units, dimension=None):
        self.units = units

        # Lookup dimension if not given
        if dimension:
            if dimension not in DIMENSIONS:
                raise TypeError(f"Unknown dimension {dimension}")
            self.dimension = dimension
        elif units in AMBIGUITIES:
            self.dimension = AMBIGUITIES[units]
        else:
            for k,v in DIMENSIONS.items():
                if units in v:
                    self.dimension = k
                    break
            else:
                raise TypeError(f"Unknown unit {units}")

        # Find the scale for the given units
        self.scalemap = DIMENSIONS[self.dimension]
        if units not in self.scalemap:
            raise TypeError(f"unable to find {units} in {self.dimension} units")
        self.scalebase = self.scalemap[self.units]

        # Assume correct units if stored units are arbitrary.
        # TODO: maybe have a 'strict' option
        if self.dimension == 'dimensionless':
            self.scalemap = None

    def convert(self, value, units=""):
        """
        Returns value in source units to a value in the target units.
        """
        # Note: calculating a*1 rather than simply returning a would produce
        # an unnecessary copy of the array, which in the case of the raw
        # counts array would be bad.  Sometimes copying and other times
        # not copying is also bad, but copy on modify semantics isn't
        # supported.
        if not units or self.scalemap is None:
            return value

        if units not in self.scalemap:
            raise KeyError(f"unable to find {units} in {self.dimension} units")
        if self.dimension == 'temperature':
            inscale, inoffset = self.scalebase
            outscale, outoffset = self.scalemap[units]
            return (value + inoffset)*(inscale/outscale) - outoffset
        else:
            inscale, outscale = self.scalebase, self.scalemap[units]
            return value * (inscale / outscale)
    __call__ = convert

    def conversion(self, units):
        """
        Returns a conversion function that takes values in the source units and
        produces values in the target units.
        """
        if not units or self.scalemap is None:
            return lambda value: value
        elif self.dimension == 'temperature':
            inscale, inoffset = self.scalebase
            outscale, outoffset = self.scalemap[units]
            return lambda value: (value + inoffset)*(inscale/outscale) - outoffset
        else:
            inscale, outscale = self.scalebase, self.scalemap[units]
            return lambda value: value * (inscale / outscale)


def _check(expect, get):
    if abs(expect - get) > 1e-10*(abs(expect)+abs(get)):
        raise ValueError("Expected %s but got %s" % (expect, get))
    #print expect,"==",get

def _checkstr(expect, get):
    if expect != get:
         raise ValueError("Expected %r but got %r"%(expect, get))

def test():
    _check(2, Converter('mm')(2000, 'm')) # 2000 mm -> 2 m
    _check(0.003, Converter('microseconds')(3, units='ms')) # 3 us -> 0.003 ms
    # If no units requested return the original units
    _check(45, Converter('nanokelvin')(45))  # 45 nK -> 45 nK
    _check(0.045, Converter('nanokelvin')(45, 'uK'))  # 45 nK -> 0.045 uK
    _check(0.5, Converter('seconds')(1800, units='hours')) # 1800 -> 0.5 hr
    # Temperature conversion with offset values
    _check(32, Converter('degC')(0, 'degF')) # 0 C -> 32 F
    _check(373.15, Converter('degF')(212, 'K')) #  212 F -> 373.15 K
    _check(-40, Converter('degF')(-40, 'degC')) # -40 F = -40 C
    # Unicode test
    _check(-273.147, Converter('mK')(3, '℃')) # 3 mK -> absolute zero in C + 0.003
    _check(2, Converter('1/A')(20, 'nm^-1'))
    # Check various forms of caller
    _check(2, Converter('1/A').convert(20, 'nm^-1'))
    _check(2, Converter('1/A').conversion('nm^-1')(20))
    _check(2, convert(20, '1/A', 'nm^-1'))
    # Check that arbitrary units are ignored
    _check(123, Converter('a.u.')(123, units=''))
    _check(123, Converter('a.u.')(123, units='s'))
    _check(123, Converter('a.u.')(123, units='mm'))
    # Check that strings are ignored if units are not set
    _checkstr('string', Converter(None)('string', None))
    _checkstr('string', Converter(None)('string', ''))
    _checkstr('string', Converter('')('string', None))
    _checkstr('string', Converter('')('string', ''))
    # TODO: more conversion tests

    try:
        Converter('help')
    except TypeError:
        pass
    else:
        raise Exception("unknown unit did not raise a type error")

    # TODO: more 
    try:
        import numpy as np
        import h5py
        file = h5py.File('test.nxs', 'w', driver='core', backing_store=False)
        file['/str'] = np.array(b'hello')
        file['/str1'] = np.array([b'hello'])
        file['/strn'] = np.array([b'hello', b'world'])
        _checkstr('hello', data_as(file, '/str', 'str'))
        _checkstr('hello', data_as(file, '/str1', 'str'))
        assert (np.array(['hello','world']) == data_as(file, '/strn', 'str')).all()
        file['/orientation'] = [35, 40]
        file['/orientation'].attrs['units'] = b'degrees'
        assert (np.radians([35, 40]) == data_as(file, '/orientation', 'radians')).all()
    except ImportError:
        pass

if __name__ == "__main__":
    test()

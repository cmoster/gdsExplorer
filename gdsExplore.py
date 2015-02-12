#

"""
    created 2015-1-27
"""

__author__ = 'maxinlong'

import struct
import gdspy
import warnings
import numpy

def _eight_byte_real(value):
    """
    Convert a number into the GDSII 8 byte real format.

    Parameters
    ----------
    value : number
        The number to be converted.

    Returns
    -------
    out : string
        The GDSII binary string that represents ``value``.
    """
    byte1 = 0
    byte2 = 0
    short3 = 0
    long4 = 0
    if value != 0:
        if value < 0:
            byte1 = 0x80
            value = -value
        exponent = int(numpy.floor(numpy.log2(value) * 0.25))
        mantissa = long(value * 16L**(14 - exponent))
        while mantissa >= 72057594037927936L:
            exponent += 1
            mantissa = long(value * 16L**(14 - exponent))
        byte1 += exponent + 64
        byte2 = (mantissa // 281474976710656L)
        short3 = (mantissa % 281474976710656L) // 4294967296L
        long4 = mantissa % 4294967296L
    return struct.pack(">HHL", byte1 * 256 + byte2, short3, long4)
def _eight_byte_real_to_float(value):
    """
    Convert a number from GDSII 8 byte real format to float.

    Parameters
    ----------
    value : string
        The GDSII binary string representation of the number.

    Returns
    -------
    out : float
        The number represented by ``value``.
    """
    short1, short2, long3 = struct.unpack('>HHL', value)
    exponent = (short1 & 0x7f00) // 256
    mantissa = (((short1 & 0x00ff) * 65536L + short2) * 4294967296L + long3) / 72057594037927936.0
    return (-1 if (short1 & 0x8000) else 1) * mantissa * 16L ** (exponent - 64)

class GdsExplore(gdspy.GdsImport):
    _rec_type_dict = { 0x00: 'version',
                       0x01: 'LIB',
                       0x02: 'LIBNAME',
                       0x03: 'Unit',
                       0x04: 'ENDSLIB',
                       0x05: 'Struct',
                       0x06: 'STRName',
                       0x07: 'ENDSTR',
                       0x08: 'Boundary',
                       0x09: 'PATH',
                       0x0A: 'SREF',
                       0x0B: 'AREF',
                       0x0C: 'TEXT',
                       0x0D: 'LAYER',
                       0x0E: 'DATATYPE',
                       0x0F: 'WIDTH',
                       0x10: 'XY',
                       0x11: 'ENDEL'}

    _dataType_recType_dict = {0x00: (0x04,0x07,0x08,0x09,0x0A,0x0B,0x11 ),
                              0x01: (0x17,0x1A,0x26),
                              0x02: (0x00,0x01,0x05,0x0d,0x0e,0x13,0x16),
                              0x03: (0x0f,0x10),
                              0x04: (),
                              0x05: (0x03,0x1b,0x1c),
                              0x06: (0x02,0x06,0x12,0x1f)}

    def __init__(self):
        self.cell_dict = {}
        self._incomplete = []





    def _read_record(self, stream):
        """
        Read a complete record from a GDSII stream file.

        Parameters
        ----------
        stream : file
            GDSII stream file to be imported.

        Returns
        -------
        out : 2-tuple
            Record type and data (as a numpy.array)
        """
        header = stream.read(4)
        if len(header) < 4:
            return None
        size, rec_type = struct.unpack('>HH', header)
        data_type = (rec_type & 0x00ff)
        rec_type = rec_type // 256
        rec_data = None

        if rec_type in self._dataType_recType_dict[data_type]:
            rec_data = self._read_raw_data(stream, data_type,size)
            print "rec {0:2x}{1:2x} [{2}] is read".format(rec_type,data_type, self._rec_type_dict[rec_type])
        else:
            warnings.warn("record type {0:2x} and data type {1:2x} not consistence, data(size {2}) dumped".format(rec_type,data_type,size))
            rec_dat = self._read_raw_data(stream, data_type,size)
        return [rec_type,data_type,rec_data]

    def _read_raw_data(self, stream,dataType, size):
        data = None
        if size > 4:
            if dataType == 0x01:
                data = numpy.array(struct.unpack('>{0}H'.format((size - 4) // 2), stream.read(size - 4)), dtype='uint')
            elif dataType == 0x02:
                data = numpy.array(struct.unpack('>{0}h'.format((size - 4) // 2), stream.read(size - 4)), dtype=int)
            elif dataType == 0x03:
                data = numpy.array(struct.unpack('>{0}l'.format((size - 4) // 4), stream.read(size - 4)), dtype=int)
            elif dataType == 0x05:
                data = numpy.array([_eight_byte_real_to_float(stream.read(8)) for _ in range((size - 4) // 8)])
            else:
                data = stream.read(size - 4)
                if data[-1] == '\0':
                    data = data[:-1]
        return  data



# test

gdsfile = '/home/maxinlong/Documents/work/design/gds/out_da.gds'

gdsStream = open(gdsfile)

gdsExp = GdsExplore()


data = gdsExp._read_record(gdsStream)

while data is not None:
    data = gdsExp._read_record(gdsStream)


gdsStream.close()
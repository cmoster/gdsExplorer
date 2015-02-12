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

class Polygon(gdspy.Polygon):
    def __init__(self, points=[], layer=0, datatype=0, verbose=True) :
        if len(points) > 199:
            if len(points) > 4094:
                raise ValueError("[GDSPY] Polygons with more than 4094 are not supported by the GDSII format.")
            elif verbose:
                warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)
        self.layer = layer
        self.points = numpy.array(points)
        self.datatype = datatype

    def setPoints(self,points):
        if len(points) > 199:
            if len(points) > 4094:
                raise ValueError("[GDSPY] Polygons with more than 4094 are not supported by the GDSII format.")
            elif verbose:
                warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)
        self.points = numpy.array(points[:-2].reshape((points.size // 2 - 1, 2)))

    def setLayer(self,layer):
        self.layer = layer

    def setDatatype(self,datatype):
        self.datatype = datatype

class CellReference(gdspy.CellReference):
    def setGeoPar(self,):
        print 'cell-ref'


class GdsExplore(gdspy.GdsImport):
    # _rec_type_dict = { 0x00: 'version',
    #                    0x01: 'LIB',
    #                    0x02: 'LIBNAME',
    #                    0x03: 'Unit',
    #                    0x04: 'ENDSLIB',
    #                    0x05: 'Struct',
    #                    0x06: 'STRName',
    #                    0x07: 'ENDSTR',
    #                    0x08: 'Boundary',
    #                    0x09: 'PATH',
    #                    0x0A: 'SREF',
    #                    0x0B: 'AREF',
    #                    0x0C: 'TEXT',
    #                    0x0D: 'LAYER',
    #                    0x0E: 'DATATYPE',
    #                    0x0F: 'WIDTH',
    #                    0x10: 'XY',
    #                    0x11: 'ENDEL',
    #                    0x12: 'SNAME' }

    _dataType_recType_dict = {0x00: (0x04,0x07,0x08,0x09,0x0A,0x0B,0x11 ),
                              0x01: (0x17,0x1A,0x26),
                              0x02: (0x00,0x01,0x05,0x0d,0x0e,0x13,0x16),
                              0x03: (0x0f,0x10),
                              0x04: (),
                              0x05: (0x03,0x1b,0x1c),
                              0x06: (0x02,0x06,0x12,0x1f)}

    def __init__(self, infile,unit=None, rename={}, layers={}, datatypes={}, texttypes={}, verbose=True):
        self.cell_dict = {}
        self._incomplete = []
        if infile.__class__ == ''.__class__:
            infile = open(infile, 'rb')
            close = True
        else:
            close = False
        emitted_warnings = []

        ascii_dec = lambda rec: rec[:-1].decode('ascii') if rec[-1] == 0 else rec.decode('ascii')

        record = self._read_record(infile)
        factor = 1e-6;
        element = None
        cell = None
        while record is not None:

            if cell is None:
                if record[0] == 0x03:
                    if unit is None:
                        factor = record[1][0]
                    else:
                        factor = record[1][1] / unit
                elif record[0] == 0x06:             #STRNAME
                    if not str is bytes:
                        # if record[1][-1] == 0:
                        #     record[1] = record[1][:-1].decode('ascii')
                        # else:
                        #     record[1] = record[1].decode('ascii')
                        record[1] = ascii_dec(record[1])
                    # name = rename.get(record[1], record[1])
                    name = record[1]
                    cell = gdspy.Cell(name, exclude_from_global=False)
                    self.cell_dict[name] = cell
                    print 'create cell [ ' + name + ' ]'
                # else:
                #     print 'other information'
            elif record[0] == 0x07:
                cell = None
                print 'cell closed'
            elif record[0] == 0x06:
                cell = None
                print 'cell closed un-normally'
            else:
                if element is None:
                    if record[0] == 0x08:           # create boundary element
                        element = Polygon()
                        cell.add(element)
                        print 'polygon created'
                    elif record[0] == 0x12:
                        # element = cellreference()
                        print 'non-polygon element '
                elif record[0] == 0x11:
                    element = None
                    print 'element closed normally'
                else:
                    if record[0] == 0x0d:
                        element.setLayer(layers.get(record[1][0], record[1][0]))
                    elif record[0] == 0x10:
                        element.setPoints( factor * record[1])
                    elif record[0] == 0x0e:
                        element.setDatatype(datatypes.get(record[1][0], record[1][0]))
                    else:
                        element = None
                        print 'illegal parameter, element closed un-normally'
            record = self._read_record(infile)



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
        rec_type = ( rec_type & 0xff00) // 256
        rec_data = None
        # print "rec {0:2x}{1:2x}  ".format(rec_type,data_type)

        if self.typeLegal(rec_type,data_type):
            print "                 rec {0:2x}{1:2x} [{2}] is read: size[{3}]".format(rec_type,data_type, self._record_name[rec_type],size-4)
            rec_data = self._read_raw_data(stream, data_type,size)
            return [rec_type,rec_data]
        else:
            # warnings.warn("record type {0:2x} and data type {1:2x} not consistence, data(size: {2}) dumped".format(rec_type,data_type,size-4))
            print "                 record_type {0:2x} & data_type {1:2x} not consistence, data(size: {2}) dumped".format(rec_type,data_type,size-4)
            self._dump_data(stream, size-4)
            return [0x3c,None]



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

    def _dump_data(self,stream,size):
        stream.read(size)
        return None

    def typeLegal(self,rec_type, data_type):
        if data_type <7:
            if rec_type in self._dataType_recType_dict[data_type]:
                return True
            else:
                return False
        else:
            return False



# test

# nameGdsfile = 'RFFE_9625p_dev_hfss'
nameGdsfile = 'va9628_top_analog_85_revb_hfss'
gdsfile = '/home/maxinlong/Documents/work/design/gds/' + nameGdsfile + '.gds'
gdsfileRev = '/home/maxinlong/Documents/work/design/gds/' + nameGdsfile + '_rev.gds'
# gdsfile = '/home/maxinlong/Documents/work/design/gds/out_da.gds'
# gdsStream = open(gdsfile)

gdsExp = GdsExplore(gdsfile)

print gdsExp

gdspy.gds_print(gdsfileRev)
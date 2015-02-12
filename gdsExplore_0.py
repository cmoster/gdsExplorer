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
        mantissa = long(value * 16L ** (14 - exponent))
        while mantissa >= 72057594037927936L:
            exponent += 1
            mantissa = long(value * 16L ** (14 - exponent))
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


#
# class Polygon(gdspy.Polygon):
# def __init__(self, points=[], layer=0, datatype=0, verbose=True) :
#         if len(points) > 199:
#             if len(points) > 4094:
#                 raise ValueError("[GDSPY] Polygons with more than 4094 are not supported by the GDSII format.")
#             elif verbose:
#                 warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)
#         self.layer = layer
#         self.points = numpy.array(points)
#         self.datatype = datatype
#
#     def setPoints(self,points):
#         if len(points) > 199:
#             if len(points) > 4094:
#                 raise ValueError("[GDSPY] Polygons with more than 4094 are not supported by the GDSII format.")
#             elif verbose:
#                 warnings.warn("[GDSPY] A polygon with more than 199 points was created (not officially supported by the GDSII format).", stacklevel=2)
#         self.points = numpy.array(points[:-2].reshape((points.size // 2 - 1, 2)))
#
#     def setLayer(self,layer):
#         self.layer = layer
#
#     def setDatatype(self,datatype):
#         self.datatype = datatype
#
# class CellReference(gdspy.CellReference):
#     def setGeoPar(self,):
#         print 'cell-ref'


class GdsExplore(gdspy.GdsImport):
    _ele_dataType_dict = {0x08: 'verbose',
                          0x0d: 'layer',
                          0x0e: 'datatype',
                          0x16: 'texttype',
                          0x10: 'xy',
                          # 0x04: 'ENDSLIB',
                          0x0f: 'width',
                          0x12: 'ref_cell',
                          0x13: 'colrow',
                          0x1a: 'x_reflection',
                          0x1b: 'magnification',
                          0x1c: 'rotation',
                          0x17: 'anchor',
                          0x21: 'ends'}

    _dataType_recType_dict = {0x00: (0x04, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x11 ),
                              0x01: (0x17, 0x1A, 0x26),
                              0x02: (0x00, 0x01, 0x05, 0x0d, 0x0e, 0x13, 0x16),
                              0x03: (0x0f, 0x10),
                              0x04: (),
                              0x05: (0x03, 0x1b, 0x1c),
                              0x06: (0x02, 0x06, 0x12, 0x1f)}


    def __init__(self, infile, unit=None, rename={}, layers={}, datatypes={}, texttypes={}, verbose=True):
        self.cell_dict = {}
        self._incomplete = []
        self._rec_element_dict = {0x08: self._create_polygon,
                                  0x09: self._create_path,
                                  0x0a: self._create_reference,
                                  0x0b: self._create_array,
                                  0x0c: self._create_label}

        if infile.__class__ == ''.__class__:
            infile = open(infile, 'rb')
            close = True
        else:
            close = False
        emitted_warnings = []
        self.factor = 1e-6
        create_element = None
        element = None
        cell_to_build = None
        kwargs = {}
        rec_method = None

        record = self._read_record(infile)

        while record is not None:
            rec_method = self._data_process(record[0])
            if record[0] == 0x03:  # UNIT
                if unit is None:
                    self.factor = record[1][0]
                else:
                    self.factor = record[1][1] / unit
            elif cell_to_build is None:
                if record[0] == 0x06:  #STRNAME  create structure
                    # if not str is bytes:
                    # if record[1][-1] == 0:
                    #     record[1] = record[1][:-1].decode('ascii')
                    # else:
                    #     record[1] = record[1].decode('ascii')
                    # record[1] = ascii_dec(record[1])
                    # name = rename.get(record[1], record[1])
                    name = rec_method(record[1])
                    cell_to_build = gdspy.Cell(name, exclude_from_global=True)
                    self.cell_dict[name] = cell_to_build
                    # else:
                    #     print 'other information'
                elif  record[0] == 0x04:         #ENDLIB
                    # print 'end of lib'
                    for ref in self._incomplete:
                        if ref.ref_cell in self.cell_dict:
                            ref.ref_cell = self.cell_dict[ref.ref_cell]
                            # print 'sref found in self.cell_dict'
                        elif ref.ref_cell in gdspy.Cell.cell_dict:
                            # print 'sref found in Cell.cell_dict'
                        # else:
                            ref.ref_cell = gdspy.Cell.cell_dict.get(ref.ref_cell, ref.ref_cell)
            elif record[0] == 0x07:         #ENDSTR
                print 'create cell [{} ]'.format(self.cell_dict.get(name))
                cell_to_build = None
                # print 'cell closed'

            else:
                if create_element is None:
                    create_element = self._rec_element_dict.get(record[0], record[0])
                elif record[0] == 0x11:     #END_ELEMENT
                    element = create_element(**kwargs)
                    cell_to_build.add(element)
                    # print 'element {} created'.format(self._record_name[record[0]])
                    element = None
                    create_element = None
                    kwargs = {}
                else:
                    ele_data_type = GdsExplore._ele_dataType_dict.get(record[0])
                    if not (ele_data_type is None or rec_method is None):
                        kwargs[ele_data_type] = rec_method(record[1])

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

        if self.typeLegal(rec_type, data_type):
            # print "                 rec {0:2x}{1:2x} [{2}] is read: size[{3}]".format(rec_type, data_type,self._record_name[rec_type],size - 4)
            rec_data = self._read_raw_data(stream, data_type, size)
            return [rec_type, rec_data]
        else:
            # warnings.warn("record type {0:2x} and data type {1:2x} not consistence, data(size: {2}) dumped".format(rec_type,data_type,size-4))
            print "                 record_type {0:2x} & data_type {1:2x} not consistence, data(size: {2}) dumped".format(
                rec_type, data_type, size - 4)
            self._dump_data(stream, size - 4)
            return [0x3c, None]

    def _read_raw_data(self, stream, dataType, size):
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
        return data

    def _dump_data(self, stream, size):
        stream.read(size)
        return None

    def typeLegal(self, rec_type, data_type):
        if data_type < 7:
            if rec_type in GdsExplore._dataType_recType_dict[data_type]:
                return True
            else:
                return False
        else:
            return False

    def _data_process(self, rec_type):

        if rec_type in (0x0d, 0x0e, 0x16, 0x1b, 0x1c):
            return (lambda rec: rec[0])  # record[1][0]
        elif rec_type == 0x13:
            return (lambda rec: rec)
        elif rec_type in (0x10, 0x1f):
            return (lambda rec: self.factor * (rec ))
        elif rec_type in (0x12, 0x19, 0x06):
            return (lambda rec: rec[:-1].decode('ascii') if rec[-1] == 0 else rec.decode('ascii'))
        elif rec_type == 0x1a:                      # x_reflection
            return (lambda rec: ((int(rec[0]) & 0x8000) > 0))
        elif rec_type == 0x17:                      # anchor
            return (lambda rec: ['nw', 'n', 'ne', None, 'w', 'o', 'e', None, 'sw', 's', 'se'][int(rec[0]) & 0x000f])
        else:
            # print 'no method for such type data'
            return None

    def _create_array(self, **kwargs):
        col_row = kwargs.pop('colrow')
        kwargs['columns'] = col_row[0]
        kwargs['rows'] = col_row[1]
        xy = kwargs.pop('xy')
        kwargs['origin'] = xy[0:2]
        if 'x_reflection' in kwargs:
            if 'rotation' in kwargs:
                sa = -numpy.sin(kwargs['rotation'] * numpy.pi / 180.0)
                ca = numpy.cos(kwargs['rotation'] * numpy.pi / 180.0)
                x2 = (xy[2] - xy[0]) * ca - (xy[3] - xy[1]) * sa + xy[0]
                y3 = (xy[4] - xy[0]) * sa + (xy[5] - xy[1]) * ca + xy[1]
            else:
                x2 = xy[2]
                y3 = xy[5]
            if kwargs['x_reflection']:
                y3 = 2 * xy[1] - y3
            kwargs['spacing'] = ((x2 - xy[0]) / kwargs['columns'], (y3 - xy[1]) / kwargs['rows'])
        else:
            kwargs['spacing'] = ((xy[2] - xy[0]) / kwargs['columns'], (xy[5] - xy[1]) / kwargs['rows'])
        ref = gdspy.CellArray(**kwargs)
        if not isinstance(ref.ref_cell, gdspy.Cell):
            self._incomplete.append(ref)
        return ref


    def flatten_by_layer(self):
        single_polygon = None
        single_cell = None
        cell_dict ={}
        layer = 0
        elem_by_layer ={}
        for single_cell in self.top_level():
            # print isinstance(single_cell,gdspy.Cell)
            polygon_dic = single_cell.get_polygons(by_spec=True)
            for key in polygon_dic.iterkeys():
                layer = key[0]
                for num_polygon in range(len(polygon_dic[key])):
                    if elem_by_layer.has_key(layer):
                        elem_by_layer[layer].append(gdspy.Polygon(polygon_dic[key][num_polygon],layer))
                    else:
                        elem_by_layer[layer] = [gdspy.Polygon(polygon_dic[key][num_polygon],layer)]

        for layer in elem_by_layer.iterkeys():
            print 'layer {}, {} elements'.format(layer, len(elem_by_layer[layer]))

        for layer in elem_by_layer.iterkeys():
            # if layer > 3: pass
            # else:
                cell_name  = 'layer' + str(layer)
                union_cell = gdspy.Cell(cell_name, exclude_from_global=True)
                cell_dict[cell_name] = union_cell

                elem_list = elem_by_layer[layer]
                print cell_name +  '  ' + str(len(elem_list))
                elem_test = elem_list.pop(0)
                united = False
                while len(elem_list) > 0:
                    for num_elem in range(len(elem_list)):
                        # print ' num_elem ' + str(num_elem)
                        union_polygonset = gdspy.boolean([elem_test,elem_list[num_elem]],lambda x,y: x or y, layer=layer,eps=5e-9)
                        if union_polygonset is None or len(union_polygonset.polygons) >1:
                            united = False
                            # union_cell.add(elem_test)
                            # print union_cell
                        else:
                             elem_list[num_elem] = gdspy.Polygon(union_polygonset.polygons[0], layer=layer)
                             united = True
                             break
                    if not united:
                        union_cell.add(elem_test)
                        print union_cell
                    elem_test = elem_list.pop(0)
                    # print str(len(elem_list))
                union_cell.add(elem_test)
                print union_cell

        # union_cell = gdspy.Cell('layer1',exclude_from_global=True)
        # cell_dict['layer1'] = union_cell
        # for element in elem_by_layer[2]:
        #     union_cell.add(element)

        return  cell_dict









# def checkCell(cells).:
#     for cell in cells.itervalues():
#
#         if isinstance(cell, gdspy.Cell):
#             pass
#         else:
#             print cell + '  ' + repr(isinstance(cell,gdspy.Cell))
#
#
# def checkCellRef(cells):
#     for cell in cells.itervalues():
#         for c in cell.get_dependencies():
#             # if not isinstance(c, gdspy.Cell):
#                 print c
#

# test

# nameGdsfile = 'RFFE_9625p_dev_hfss'
nameGdsfile = 'va9628_top_analog_85_revb_hfss'
# nameGdsfile = 'Layout1'
gdsfile = '/home/maxinlong/Documents/work/design/gds/' + nameGdsfile + '.gds'
gdsfileRev = '/home/maxinlong/Documents/work/design/gds/' + nameGdsfile + '_uni.gds'
# gdsfile = '/home/maxinlong/Documents/work/design/gds/out_da.gds'
# gdsStream = open(gdsfile)

gdsExp = GdsExplore(gdsfile)
print '***********************'

# for cell in gdsExp._incomplete:
#     print cell.ref_cell

# print gdspy.Cell.cell_dict

# checkCellRef(gdsExp.cell_dict)

cell_union = gdsExp.flatten_by_layer()



gdspy.gds_print(gdsfileRev,{'layer15': cell_union['layer15']})


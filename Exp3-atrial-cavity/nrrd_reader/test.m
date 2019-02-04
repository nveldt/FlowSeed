% Testing nrrdread
% -------------------------------------------------------------------------
[data, metadata] = nrrdread('data/test1d_ascii.nrrd');

assert(all(data == (1:27)'), 'Invalid data matrix for test1d');
assert(metadata.dimension == 1, 'Dimension is not 1 for test1d');
assert(strcmp(metadata.type, 'uint8'), 'Type is not uint8 for test1d');
assert(metadata.sizes == 27, 'Vector length is not 27 for test1d');
assert(strcmp(metadata.encoding, 'ascii'), 'Not ASCII encoding for test1d');


[data, metadata] = nrrdread('data/test2d_ascii.nrrd');

assert(all(all(data == reshape(1:27, [3 9]))), 'Invalid data matrix for test2d');
assert(metadata.dimension == 2, 'Dimension is not 2 for test2d');
assert(strcmp(metadata.type, 'uint16'), 'Type is not uint16 for test2d');
assert(all(metadata.sizes == [3 9]), 'Sizes not right for test2d');
assert(strcmp(metadata.encoding, 'ascii'), 'Not ASCII encoding for test2d');
assert(metadata.spacedimension == 2, 'Space dimension is not 2 for test2d');
assert(all(metadata.spacings == [1.0458 2]), 'Spacing not correct for test2d');
assert(all(isnan(metadata.spacedirections)), 'Space directions not correct for test2d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test2d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test2d');
assert(all(metadata.spaceorigin == [100 200]), 'Space origin not correct for test2d');


[data, metadata] = nrrdread('data/test3d_ascii.nrrd');

assert(all(all(all(data == reshape(1:27, [3 3 3])))), 'Invalid data matrix for test3d');
assert(metadata.dimension == 3, 'Dimension is not 3 for test3d');
assert(strcmp(metadata.type, 'uint32'), 'Type is not uint32 for test3d');
assert(all(metadata.sizes == [3 3 3]), 'Sizes not right for test3d');
assert(strcmp(metadata.encoding, 'ascii'), 'Not ASCII encoding for test3d');
assert(strcmp(metadata.space, 'left-posterior-superior'), 'Space not right for test3d');
assert(all(all(metadata.spacedirections == [1 0 0; 0 1 0; 0 0 1])), 'Space directions not correct for test3d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test3d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test3d');
assert(all(metadata.spaceorigin == [100.1 200.3 -500]), 'Space origin not correct for test3d');


[data, metadata] = nrrdread('data/test1d_raw.nrrd');

assert(all(data == (1:27)'), 'Invalid data matrix for test1d');
assert(metadata.dimension == 1, 'Dimension is not 1 for test1d');
assert(strcmp(metadata.type, 'uint8'), 'Type is not uint8 for test1d');
assert(metadata.sizes == 27, 'Vector length is not 27 for test1d');
assert(strcmp(metadata.encoding, 'raw'), 'Not raw encoding for test1d');


[data, metadata] = nrrdread('data/test2d_raw.nrrd');

assert(all(all(data == reshape(1:27, [3 9]))), 'Invalid data matrix for test2d');
assert(metadata.dimension == 2, 'Dimension is not 2 for test2d');
assert(strcmp(metadata.type, 'uint16'), 'Type is not uint16 for test2d');
assert(all(metadata.sizes == [3 9]), 'Sizes not right for test2d');
assert(strcmp(metadata.encoding, 'raw'), 'Not raw encoding for test2d');
assert(metadata.spacedimension == 2, 'Space dimension is not 2 for test2d');
assert(all(metadata.spacings == [1.0458 2]), 'Spacing not correct for test2d');
assert(all(isnan(metadata.spacedirections)), 'Space directions not correct for test2d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test2d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test2d');
assert(all(metadata.spaceorigin == [100 200]), 'Space origin not correct for test2d');


[data, metadata] = nrrdread('data/test3d_raw.nrrd');

assert(all(all(all(data == reshape(1:27, [3 3 3])))), 'Invalid data matrix for test3d');
assert(metadata.dimension == 3, 'Dimension is not 3 for test3d');
assert(strcmp(metadata.type, 'uint32'), 'Type is not uint32 for test3d');
assert(all(metadata.sizes == [3 3 3]), 'Sizes not right for test3d');
assert(strcmp(metadata.encoding, 'raw'), 'Not raw encoding for test3d');
assert(strcmp(metadata.space, 'left-posterior-superior'), 'Space not right for test3d');
assert(all(all(metadata.spacedirections == [1 0 0; 0 1 0; 0 0 1])), 'Space directions not correct for test3d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test3d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test3d');
assert(all(metadata.spaceorigin == [100.1 200.3 -500]), 'Space origin not correct for test3d');


% Test big endian for ASCII - should have no effect
[data, metadata] = nrrdread('data/test3d_ascii.nrrd');
metadata.endian = 'big';

assert(all(all(all(data == reshape(1:27, [3 3 3])))), 'Invalid data matrix for test3d');
assert(metadata.dimension == 3, 'Dimension is not 3 for test3d');
assert(strcmp(metadata.type, 'uint32'), 'Type is not uint32 for test3d');
assert(all(metadata.sizes == [3 3 3]), 'Sizes not right for test3d');
assert(strcmp(metadata.encoding, 'ascii'), 'Not ASCII encoding for test3d');
assert(strcmp(metadata.space, 'left-posterior-superior'), 'Space not right for test3d');
assert(all(all(metadata.spacedirections == [1 0 0; 0 1 0; 0 0 1])), 'Space directions not correct for test3d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test3d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test3d');
assert(all(metadata.spaceorigin == [100.1 200.3 -500]), 'Space origin not correct for test3d');


% Test big endian for raw encoding
[data, metadata] = nrrdread('data/test3d_bigendian_raw.nrrd');

assert(all(all(all(data == reshape(1:27, [3 3 3])))), 'Invalid data matrix for test3d');
assert(metadata.dimension == 3, 'Dimension is not 3 for test3d');
assert(strcmp(metadata.type, 'uint32'), 'Type is not uint32 for test3d');
assert(all(metadata.sizes == [3 3 3]), 'Sizes not right for test3d');
assert(strcmp(metadata.encoding, 'raw'), 'Not raw encoding for test3d');
assert(strcmp(metadata.endian, 'big'), 'Not big endian for test3d');
assert(strcmp(metadata.space, 'left-posterior-superior'), 'Space not right for test3d');
assert(all(all(metadata.spacedirections == [1 0 0; 0 1 0; 0 0 1])), 'Space directions not correct for test3d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test3d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test3d');
assert(all(metadata.spaceorigin == [100.1 200.3 -500]), 'Space origin not correct for test3d');


% Test big endian for raw encoding, but this time no endian metadata field
% will be present. Need to force the endianness
[data, metadata] = nrrdread('data/test3d_bigendian_raw_noendianfield.nrrd', 'Endian', 'big');

assert(all(all(all(data == reshape(1:27, [3 3 3])))), 'Invalid data matrix for test3d');
assert(metadata.dimension == 3, 'Dimension is not 3 for test3d');
assert(strcmp(metadata.type, 'uint32'), 'Type is not uint32 for test3d');
assert(all(metadata.sizes == [3 3 3]), 'Sizes not right for test3d');
assert(strcmp(metadata.encoding, 'raw'), 'Not raw encoding for test3d');
assert(strcmp(metadata.endian, 'big'), 'Not big endian for test3d');
assert(strcmp(metadata.space, 'left-posterior-superior'), 'Space not right for test3d');
assert(all(all(metadata.spacedirections == [1 0 0; 0 1 0; 0 0 1])), 'Space directions not correct for test3d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test3d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test3d');
assert(all(metadata.spaceorigin == [100.1 200.3 -500]), 'Space origin not correct for test3d');


[data, metadata] = nrrdread('data/test1d_gzip.nrrd');

assert(all(data == (1:27)'), 'Invalid data matrix for test1d');
assert(metadata.dimension == 1, 'Dimension is not 1 for test1d');
assert(strcmp(metadata.type, 'uint8'), 'Type is not uint8 for test1d');
assert(metadata.sizes == 27, 'Vector length is not 27 for test1d');
assert(strcmp(metadata.encoding, 'gzip'), 'Not gzip encoding for test1d');


[data, metadata] = nrrdread('data/test2d_gzip.nrrd');

assert(all(all(data == reshape(1:27, [3 9]))), 'Invalid data matrix for test2d');
assert(metadata.dimension == 2, 'Dimension is not 2 for test2d');
assert(strcmp(metadata.type, 'uint16'), 'Type is not uint16 for test2d');
assert(all(metadata.sizes == [3 9]), 'Sizes not right for test2d');
assert(strcmp(metadata.encoding, 'gzip'), 'Not gzip encoding for test2d');
assert(metadata.spacedimension == 2, 'Space dimension is not 2 for test2d');
assert(all(metadata.spacings == [1.0458 2]), 'Spacing not correct for test2d');
assert(all(isnan(metadata.spacedirections)), 'Space directions not correct for test2d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test2d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test2d');
assert(all(metadata.spaceorigin == [100 200]), 'Space origin not correct for test2d');


[data, metadata] = nrrdread('data/test3d_gzip.nrrd');

assert(all(all(all(data == reshape(1:27, [3 3 3])))), 'Invalid data matrix for test3d');
assert(metadata.dimension == 3, 'Dimension is not 3 for test3d');
assert(strcmp(metadata.type, 'uint32'), 'Type is not uint32 for test3d');
assert(all(metadata.sizes == [3 3 3]), 'Sizes not right for test3d');
assert(strcmp(metadata.encoding, 'gzip'), 'Not gzip encoding for test3d');
assert(strcmp(metadata.space, 'left-posterior-superior'), 'Space not right for test3d');
assert(all(all(metadata.spacedirections == [1 0 0; 0 1 0; 0 0 1])), 'Space directions not correct for test3d');
assert(all(strcmp(metadata.kinds, 'domain')), 'Kinds not correct for test3d');
assert(all(strcmp(metadata.spaceunits, 'mm')), 'Space units not correct for test3d');
assert(all(metadata.spaceorigin == [100.1 200.3 -500]), 'Space origin not correct for test3d');


% Test custom fields without a custom field map
[data, metadata] = nrrdread('data/test_customFields.nrrd');

assert(all(data == (1:27)'), 'Invalid data matrix for test_customFields');
assert(metadata.dimension == 1, 'Dimension is not 1 for test_customFields');
assert(strcmp(metadata.type, 'uint8'), 'Type is not uint8 for test_customFields');
assert(metadata.sizes == 27, 'Vector length is not 27 for test_customFields');
assert(strcmp(metadata.encoding, 'ascii'), 'Not ASCII encoding for test_customFields');

% Check the fieldMap to make sure the custom fields with spacings are
% correctly handled.
x = find(strcmp(metadata.fieldMap(:, 1), 'intlist'));
assert(~isempty(x), 'intlist field is not present in fieldMap');
assert(strcmp(metadata.fieldMap(x, 2), 'int list'), 'Incorrect field mapping for int list');

x = find(strcmp(metadata.fieldMap(:, 1), 'doublelist'));
assert(~isempty(x), 'doublelist field is not present in fieldMap');
assert(strcmp(metadata.fieldMap(x, 2), 'double list'), 'Incorrect field mapping for double list');

x = find(strcmp(metadata.fieldMap(:, 1), 'stringlist'));
assert(~isempty(x), 'stringlist field is not present in fieldMap');
assert(strcmp(metadata.fieldMap(x, 2), 'string list'), 'Incorrect field mapping for string list');

x = find(strcmp(metadata.fieldMap(:, 1), 'intvector'));
assert(~isempty(x), 'intvector field is not present in fieldMap');
assert(strcmp(metadata.fieldMap(x, 2), 'int vector'), 'Incorrect field mapping for int vector');

x = find(strcmp(metadata.fieldMap(:, 1), 'doublevector'));
assert(~isempty(x), 'doublevector field is not present in fieldMap');
assert(strcmp(metadata.fieldMap(x, 2), 'double vector'), 'Incorrect field mapping for double vector');

x = find(strcmp(metadata.fieldMap(:, 1), 'intmatrix'));
assert(~isempty(x), 'intmatrix field is not present in fieldMap');
assert(strcmp(metadata.fieldMap(x, 2), 'int matrix'), 'Incorrect field mapping for int matrix');

x = find(strcmp(metadata.fieldMap(:, 1), 'doublematrix'));
assert(~isempty(x), 'doublematrix field is not present in fieldMap');
assert(strcmp(metadata.fieldMap(x, 2), 'double matrix'), 'Incorrect field mapping for double matrix');

% Check that all of the custom fields are parsed as strings because no
% custom field map is given
assert(ischar(metadata.int), 'Custom field int is not a string and should be');
assert(ischar(metadata.double), 'Custom field double is not a string and should be');
assert(ischar(metadata.string), 'Custom field string is not a string and should be');
assert(ischar(metadata.intlist), 'Custom field int list is not a string and should be');
assert(ischar(metadata.doublelist), 'Custom field double list is not a string and should be');
assert(ischar(metadata.stringlist), 'Custom field string list is not a string and should be');
assert(ischar(metadata.intvector), 'Custom field int vector is not a string and should be');
assert(ischar(metadata.doublevector), 'Custom field double vector is not a string and should be');
assert(ischar(metadata.intmatrix), 'Custom field int matrix is not a string and should be');
assert(ischar(metadata.doublematrix), 'Custom field double matrix is not a string and should be');

% Test custom fields with a custom field map
fieldMap = {'int' 'int'; 'double' 'double'; 'string' 'string'; ...
            'intlist' 'int list'; 'doublelist' 'double list'; ...
            'stringlist' 'string list'; 'intvector' 'int vector'; ...
            'doublevector' 'double vector'; 'intmatrix' 'int matrix'; ...
            'doublematrix' 'double matrix'};
[data, metadata] = nrrdread('data/test_customFields.nrrd', 'CustomFieldMap', fieldMap);

assert(all(data == (1:27)'), 'Invalid data matrix for test_customFields');
assert(metadata.dimension == 1, 'Dimension is not 1 for test_customFields');
assert(strcmp(metadata.type, 'uint8'), 'Type is not uint8 for test_customFields');
assert(metadata.sizes == 27, 'Vector length is not 27 for test_customFields');
assert(strcmp(metadata.encoding, 'ascii'), 'Not ASCII encoding for test_customFields');

% Check that all of the custom fields are parsed as strings because no
% custom field map is given
assert(metadata.int == 24, 'Custom field int is incorrect');
assert(metadata.double == 25.5566, 'Custom field double is incorrect');
assert(strcmp(metadata.string, 'this is a long string of information that is important.'), 'Custom field string is incorrect');
assert(all(metadata.intlist == [1 2 3 4 5 100]), 'Custom field intlist is incorrect');
assert(all(metadata.doublelist == [0.2 0.502 0.8]), 'Custom field doublelist is incorrect');
assert(all(strcmp(metadata.stringlist, {'words' 'are' 'split' 'by' 'space' 'in' 'list'})), 'Custom field string list is not a string and should be');
assert(all(metadata.intvector == [100 200 -300]), 'Custom field intvector is incorrect');
assert(all(metadata.doublevector == [100.5 200.3 -300.99]), 'Custom field doublevector is incorrect');
assert(all(all(metadata.intmatrix == [1 0 0; 0 1 0; 0 0 1])), 'Custom field intmatrix is incorrect');
assert(all(all(metadata.doublematrix == [1.2 0.3 0; 0 1.5 0; 0 -0.55 1.6])), 'Custom field doublematrix is incorrect');

% Testing nrrdwrite
% -------------------------------------------------------------------------
filename = 'data/test.nrrd';

[data, metadata] = nrrdread('data/test1d_ascii.nrrd');
nrrdwrite(filename, data, metadata);

% Open file
[fid, msg] = fopen(filename, 'r');
assert(fid > 3, ['Could not open file: ' msg]);

assert(strcmp(fgetl(fid), 'NRRD0005'), 'Invalid magic string');
fgetl(fid); fgetl(fid); fgetl(fid);
assert(strcmp(fgetl(fid), 'type: uint8'), 'Invalid type');
assert(strcmp(fgetl(fid), 'dimension: 1'), 'Invalid dimension');
assert(strcmp(fgetl(fid), 'sizes: 27'), 'Invalid sizes');
assert(strcmp(fgetl(fid), 'kinds: domain'), 'Invalid domain');
assert(strcmp(fgetl(fid), 'endian: little'), 'Invalid endian');
assert(strcmp(fgetl(fid), 'encoding: ascii'), 'Invalid encoding');
assert(strcmp(fgetl(fid), 'spacings: 1.0458'), 'Invalid spacings');

% TODO Still need to write documentation for nrrdwrite
% TODO Still need to fix some quotes thing for nrrdwrite


[data, metadata] = nrrdread('data/test2d_ascii.nrrd');
nrrdwrite(filename, data, metadata);

% Open file
fclose(fid);
[fid, msg] = fopen(filename, 'r');
assert(fid > 3, ['Could not open file: ' msg]);

fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
assert(strcmp(fgetl(fid), 'type: uint16'), 'Invalid type');
assert(strcmp(fgetl(fid), 'dimension: 2'), 'Invalid dimension');
assert(strcmp(fgetl(fid), 'space dimension: 2'), 'Invalid space dimension');
assert(strcmp(fgetl(fid), 'sizes: 3 9'), 'Invalid sizes ');
assert(strcmp(fgetl(fid), 'space directions: none none'), 'Invalid space directions');
assert(strcmp(fgetl(fid), 'kinds: domain domain'), 'Invalid domain');
assert(strcmp(fgetl(fid), 'endian: little'), 'Invalid endian');
assert(strcmp(fgetl(fid), 'encoding: ascii'), 'Invalid encoding');
assert(strcmp(fgetl(fid), 'spacings: 1.0458 2'), 'Invalid spacings');
assert(strcmp(fgetl(fid), 'space units: mm mm'), 'Invalid space units');
assert(strcmp(fgetl(fid), 'space origin: (100,200)'), 'Invalid space origin');
fgetl(fid);
assert(strcmp(fgetl(fid), '1 2 3'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '4 5 6'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '7 8 9'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '10 11 12'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '13 14 15'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '16 17 18'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '19 20 21'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '22 23 24'), 'Invalid ASCII data');
assert(strcmp(fgetl(fid), '25 26 27'), 'Invalid ASCII data');


[data, metadata] = nrrdread('data/test3d_ascii.nrrd');
nrrdwrite(filename, data, metadata);

% Open file
fclose(fid);
[fid, msg] = fopen(filename, 'r');
assert(fid > 3, ['Could not open file: ' msg]);

fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
assert(strcmp(fgetl(fid), 'type: uint32'), 'Invalid type');
assert(strcmp(fgetl(fid), 'dimension: 3'), 'Invalid dimension');
assert(strcmp(fgetl(fid), 'space: left-posterior-superior'), 'Invalid space');
assert(strcmp(fgetl(fid), 'sizes: 3 3 3'), 'Invalid sizes ');
assert(strcmp(fgetl(fid), 'space directions: (1,0,0) (0,1,0) (0,0,1)'), 'Invalid space directions');
assert(strcmp(fgetl(fid), 'kinds: domain domain domain'), 'Invalid domain');
assert(strcmp(fgetl(fid), 'endian: little'), 'Invalid endian');
assert(strcmp(fgetl(fid), 'encoding: ascii'), 'Invalid encoding');
assert(strcmp(fgetl(fid), 'space units: mm mm mm'), 'Invalid space units');
assert(strcmp(fgetl(fid), 'space origin: (100.1,200.3,-500)'), 'Invalid space origin');
fgetl(fid);

for i = 1:27
    assert(strcmp(fgetl(fid), num2str(i)), 'Invalid ASCII data');
end


% No custom field map is given so these are read/written as strings
[data, metadata] = nrrdread('data/test_customFields.nrrd');
nrrdwrite(filename, data, metadata);

% Open file
fclose(fid);
[fid, msg] = fopen(filename, 'r');
assert(fid > 3, ['Could not open file: ' msg]);

fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
assert(strcmp(fgetl(fid), 'type: uint8'), 'Invalid type');
assert(strcmp(fgetl(fid), 'dimension: 1'), 'Invalid dimension');
assert(strcmp(fgetl(fid), 'sizes: 27'), 'Invalid sizes ');
assert(strcmp(fgetl(fid), 'kinds: domain'), 'Invalid domain');
assert(strcmp(fgetl(fid), 'endian: little'), 'Invalid endian');
assert(strcmp(fgetl(fid), 'encoding: ascii'), 'Invalid encoding');
assert(strcmp(fgetl(fid), 'spacings: 1.0458'), 'Invalid spacings');
assert(strcmp(fgetl(fid), 'int:= 24'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double:= 25.5566'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'string:= this is a long string of information that is important.'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'int list:= 1 2 3 4 5 100'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double list:= 0.2 0.502 0.8'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'string list:= words are split by space in list'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'int vector:= (100, 200, -300)'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double vector:= (100.5,200.3,-300.99)'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'int matrix:= (1,0,0) (0,1,0) (0,0,1)'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double matrix:= (1.2,0.3,0) (0,1.5,0) (0,-0.55,1.6)'), 'Invalid custom fields');


% Custom fields with field map
fieldMap = {'int' 'int'; 'double' 'double'; 'string' 'string'; ...
            'intlist' 'int list'; 'doublelist' 'double list'; ...
            'stringlist' 'string list'; 'intvector' 'int vector'; ...
            'doublevector' 'double vector'; 'intmatrix' 'int matrix'; ...
            'doublematrix' 'double matrix'};
[data, metadata] = nrrdread('data/test_customFields.nrrd', 'CustomFieldMap', fieldMap);
nrrdwrite(filename, data, metadata, 'CustomFieldMap', fieldMap);

% Open file
fclose(fid);
[fid, msg] = fopen(filename, 'r');
assert(fid > 3, ['Could not open file: ' msg]);

fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
assert(strcmp(fgetl(fid), 'type: uint8'), 'Invalid type');
assert(strcmp(fgetl(fid), 'dimension: 1'), 'Invalid dimension');
assert(strcmp(fgetl(fid), 'sizes: 27'), 'Invalid sizes ');
assert(strcmp(fgetl(fid), 'kinds: domain'), 'Invalid domain');
assert(strcmp(fgetl(fid), 'endian: little'), 'Invalid endian');
assert(strcmp(fgetl(fid), 'encoding: ascii'), 'Invalid encoding');
assert(strcmp(fgetl(fid), 'spacings: 1.0458'), 'Invalid spacings');
assert(strcmp(fgetl(fid), 'int:= 24'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double:= 25.5566'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'string:= this is a long string of information that is important.'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'int list:= 1 2 3 4 5 100'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double list:= 0.2 0.502 0.8'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'string list:= words are split by space in list'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'int vector:= (100,200,-300)'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double vector:= (100.5,200.3,-300.99)'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'int matrix:= (1,0,0) (0,1,0) (0,0,1)'), 'Invalid custom fields');
assert(strcmp(fgetl(fid), 'double matrix:= (1.2,0.3,0) (0,1.5,0) (0,-0.55,1.6)'), 'Invalid custom fields');

fclose(fid);
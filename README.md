
Simple nastran to OF converter.  
Now only works correctly with SI units (meters).  
Tested with OpenFOAM+v2106, and an NX nastran file.

Patches and cell zones are generated based on the property card IDs.

TODO: There are some quirky solutions in the file parsing and probably some bugs...

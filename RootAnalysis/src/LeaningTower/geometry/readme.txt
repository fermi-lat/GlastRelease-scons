This directory contains several files used by the LeaningTower scripts:

geometry files:

A geometry file has up to 7 columns: name posZ shiftY shiftX rotZ rotY rotX.
The first versions of geometry files had only two columns, name and posZ.
Later, the other columns were added.  For backward compatibility, the order was
kept, i.e. z first.  If there are less than 7 entries in a row, the missing
values are assumed to be 0.

Z, the vertical position, is centered in the material (the substrate has a
thickness of 400um).

********************************************************************************

stack2geometry.txt
    describes the geometry of "stack2".
    http://glastserver.pi.infn.it/TrackerElectricalTests/Stacks/stack2.html

MCTowergeometry.txt
    file created by Nicola Omodei.  Don't know when and where it was used.  For
    sure, it's not a real tower, as X and Y planes are alternating.

Generic tower: ****************************************************************

TowerGeometryGleamv5r8.txt
    cluster positions extracted from Gleam v5r8.

TowerGeometrySaggini.txt
    numbers extracted by Nicola Saggini (engineer at Pisa) from the drawings of
    the tower.

Tower 0: **********************************************************************

Tower0Geometry.txt
    Tower0 after alignment.

Tower A: **********************************************************************

TowerAgeometry306000357.txt
    TowerA upside-down, before mounting of the side walls.

TowerAgeometry306000384.txt
    TowerA upside-down, after mounting of the side walls.

Tower B: **********************************************************************

TowerBgeometry306000510.txt
    TowerB upside-down, after mounting of the side walls.

fitting files: *****************************************************************

Tower0FittingPlanes.txt
    contains a setup where the residuals of a particular plane is being
    determined by hits in the trays (with same orientation) adjacent above and
    below.  Because several planes of tower 0 were bad or not existent, a
    canonical setup couldn't be used.

This directory contains several files used by the LeaningTower scripts:

geometry files:

A geometry file has up to 7 columns: name posZ shiftY shiftX rotZ rotY rotX

The first versions had only name and posZ.  For backward compatibility, I kept
the order this way.  If there are less than 7 entries in a column, the values
are assumed to be 0.

stack2geometry.txt:  describes the geometry of "stack2".
        http://glastserver.pi.infn.it/TrackerElectricalTests/Stacks/stack2.html

MCTowergeometry.txt:  file created by Nicola Omodei.  Don't know when and where
                      it was used.  For sure, it's not a real tower, as X and Y
                      planes are alternating.

TowerGeometryGleamv5r8.txt:  cluster positions extracted from Gleam v5r8.

TowerGeometrySaggini.txt:  numbers were extracted by Nicola Saggini (engineer
                           at Pisa) from the drawings of the tower.

Tower0Geometry.txt:  contains the plane positions (center of material) of Tower0
 after alignment.

TowerAgeometry357.txt:  contains the plane positions (center of material) of
                        TowerA before mounting of the side walls, in upside-down
                        position (alignment not yet done).   Run #306000357

TowerAgeometry384.txt:  contains the plane positions (center of material) of
                        TowerA after mounting of the side walls, in upside-down
                        position.   Run #306000384

TowerAgeometry.txt:  contains the plane positions (center of material) of TowerA
                     (alignment not yet done).

fitting files:

Tower0FittingPlanes.txt:  contains a setup where the residuals of a particular
                          plane is being determined by hits in the trays (with
                          same orientation) adjacent above and below.  Because
                          several planes of tower 0 were bad or not existent, a
                          canonical setup couldn't be used.

This directory contains several files used by the LeaningTower scripts:

geometry files:

A geometry file has up to 7 columns: name posZ shiftY shiftX rotZ rotY rotX.
The first versions of geometry files had only two columns, name and posZ.
Later, the other columns were added.  For backward compatibility, the order was
kept, i.e. z first.  If there are less than 7 entries in a row, the missing
values are assumed to be 0.

Z, the vertical position, is centered in the material (the substrate has a
thickness of 400um).

For historical reasons, what is called tower here, is in reality a TkrFM.
Instead, bay describes the position of a real tower (TkrFM + CAL) in the grid.
The corresponces are:

Bay position   TkrFM
       0         A
       1         2
       2        14
       3        15
       4         B
       5         1
       6        12
       7        13
       8         5
       9         3
      10         7
      11         9
      12         6
      13         4
      14        10
      15        11

Full LAT **********************************************************************

Bay*geometry135005524.txt are retrieved from the cat of runs 135005518, 20, 22,
and 24.

*******************************************************************************

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
    upside-down, before mounting of the side walls.

TowerAgeometry306000382.txt
TowerAgeometry306000384.txt
TowerAgeometry306000385.txt
    upside-down, after mounting of the side walls.

TowerAgeometry303002054.txt
TowerAgeometry308000281.txt
TowerAgeometry308000326.txt
TowerAgeometry308000333.txt
TowerAgeometry308000448.txt
TowerAgeometry308000449.txt
TowerAgeometry308000576.txt
TowerAgeometry308000583.txt
    runs taken at Alenia.  The conditions can be retrieved from
    http://www.pi.infn.it/~kuss/public/aleniaHistory.html.

TowerAgeometry306000445.txt
    upright, after being back from Alenia

TowerAgeometry306000475.txt
    upright, after removing and remounting one side-wall to fix cable 6.

Tower B: **********************************************************************

TowerBgeometry306000483.txt
    upside-down, before mounting of the side walls, with 3 cables still not
    connected.

TowerBgeometry306000510.txt
TowerBgeometry306000517.txt
    upside-down, after mounting of the side walls.




fitting files: ****************************************************************

Tower0FittingPlanes.txt
    contains a setup where the residuals of a particular plane is being
    determined by hits in the trays (with same orientation) adjacent above and
    below.  Because several planes of tower 0 were bad or not existent, a
    canonical setup couldn't be used.

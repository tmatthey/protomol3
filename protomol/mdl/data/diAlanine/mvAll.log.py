self.loadMacro('showHideDejaVuGUI', '/afs/nd.edu/coursesp.06/cse/cse40531.01/Public/pmv/share/lib/python2.3/site-packages/Pmv/Macros/DejaVuMac', cascade='DejaVuMac', menuEntry='showHideDejaVuGUI', menuBar='menuRoot', log=0, menuButton='Macros')
self.setUserPreference(('transformationLogging', 'final'))
self.setUserPreference(('showProgressBar', 'show'))
self.browseCommands('mdlCommands', commands=None, log=0, package='Pmv')
self.readMolecule('/afs/nd.edu/coursesp.06/cse/cse40531.01/Public/mdl/examples/diAlanine/blockdialanine_eq.pdb', ask=0, parser=None, log=0)
##
## Saving State for Viewer
##

## Light Model
## End Light Model

## Light sources
## End Light sources 7

## Cameras
## Camera Number 0
state = {'color': (0.0, 0.0, 0.0, 1.0), 'height': 394, 'lookAt': [0.0, 0.0, 0.0], 'rootx': 506, 'pivot': [0.0, 0.0, 0.0], 'translation': [0.0, 0.0, 0.0], 'sideBySideTranslation': 0.0, 'fov': 40.0, 'scale': [1.0, 1.0, 1.0], 'stereoMode': 'MONO', 'width': 865, 'sideBySideRotAngle': 3.0, 'boundingbox': 0, 'projectionType': 0, 'contours': False, 'direction': [0.0, 0.0, -30.0], 'far': 50.0, 'lookFrom': [0.0, 0.0, 30.0], 'antialiased': False, 'rotation': [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], 'near': 0.10000000000000001, 'rooty': 52}
apply(self.GUI.VIEWER.cameras[0].Set, (), state)

state = {'end': 40, 'density': 0.10000000000000001, 'color': (0.0, 0.0, 0.0, 1.0), 'enabled': 1, 'start': 25, 'mode': 'GL_LINEAR'}
apply(self.GUI.VIEWER.cameras[0].fog.Set, (), state)

## End Cameras

## Clipping planes
## End Clipping planes

## Root object
state = {'blendFunctions': ('GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'), 'stippleLines': False, 'scissorAspectRatio': 1.0, 'immediateRendering': False, 'shading': 'smooth', 'pivot': [-1.1411449782219505, 12.272929319385195, -2.2483950219511506], 'translation': [1.1411449002464715, -12.272929879176772, 2.248394879441121], 'scissorH': 200, 'frontPolyMode': 'line', 'inheritFrontPolyMode': False, 'inheritLineStipple': 0, 'inheritShadeModel': False, 'instanceMatrices': [[1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]], 'scissorX': 0, 'scissorY': 0, 'listed': True, 'inheritPointWidth': 0, 'pickable': 1, 'pointWidth': 3.0, 'scissorW': 200, 'cull': 'back', 'stipplePolygons': False, 'pickableVertices': False, 'inheritMaterial': False, 'depthMask': 1, 'scale': [1.6384301195131299, 1.6384301195131299, 1.6384301195131299], 'lighting': False, 'inheritCulling': False, 'inheritPolygonStipple': 0, 'rotation': [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], 'transparent': False, 'outline': False, 'name': 'root', 'backPolyMode': 'line', 'visible': True, 'inheritBackPolyMode': False, 'scissor': 0, 'protected': 1, 'antiAliased': False, 'inheritLineWidth': 0, 'lineWidth': 1.0, 'inheritXform': 0}
apply(self.GUI.VIEWER.rootObject.Set, (), state)

## End Root Object

## Material for root
if self.GUI.VIEWER.rootObject:
    pass  ## needed in case there no modif
## End Materials for root

## Clipping Planes for root
if self.GUI.VIEWER.rootObject:
    self.GUI.VIEWER.rootObject.clipP = []
    self.GUI.VIEWER.rootObject.clipPI = []
    pass  ## needed in case there no modif
## End Clipping Planes for root

##
## Saving State for objects in Viewer
##

## Object root|pickSpheres
## Object root|misc
## Object root|blockdialanine_eq
## Object root|misc|pickSpheresGeom
## Object root|blockdialanine_eq|selection
## Object root|blockdialanine_eq|sticks
## Object root|blockdialanine_eq|balls
## Object root|blockdialanine_eq|cpk
## Object root|blockdialanine_eq|lines
## Object root|misc|pickSpheresGeom|pickSpheres
## Object root|blockdialanine_eq|lines|bonded
## Object root|blockdialanine_eq|lines|nobnds
## Object root|blockdialanine_eq|lines|bondorder
## End Object root|blockdialanine_eq|lines|bondorder


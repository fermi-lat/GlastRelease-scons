
# Author(s):
# Luca Baldini (luca.baldini@pi.infn.it)
#
# Utility class to display clusters in the CAL.


import ROOT
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetFrameBorderMode(0)

TEXT_STYLE   = 43
TEXT_SIZE    = 18
TITLE_SIZE   = 18
LABEL_SIZE   = 16
ROOT.gStyle.SetTextFont(TEXT_STYLE)
ROOT.gStyle.SetTextSize(TEXT_SIZE)
ROOT.gStyle.SetLabelFont(TEXT_STYLE, 'XY')
ROOT.gStyle.SetLabelSize(LABEL_SIZE, 'XY')
ROOT.gStyle.SetTitleFont(TEXT_STYLE, 'XY')
ROOT.gStyle.SetTitleSize(TITLE_SIZE, 'XY')



from math import sqrt


# Miscellanea dimensions.
CSI_LENGTH = 326.0
CSI_WIDTH = 26.7
CSI_HEIGHT = 19.9
CELL_VERT_PITCH = 21.35
CELL_HOR_PITCH = 27.84
CAL_MODULE_WIDTH = 363.0
CAL_MODULE_HEIGHT = 222.0
TOWER_PITCH = 374.5
CAL_TOWER_GAP = TOWER_PITCH - CAL_MODULE_WIDTH
CAL_MODULE_TOP = -20.0
CAL_MODULE_BOTTOM = CAL_MODULE_TOP - CAL_MODULE_HEIGHT
CAL_MODULE_VCENTER = 0.5*(CAL_MODULE_TOP + CAL_MODULE_BOTTOM)
CAL_MODULE_CSI_TOP = -58.5 # This is wrong, use the CORR value below.
CAL_MODULE_CSI_TOP_CORR = CAL_MODULE_CSI_TOP + 0.5*CSI_HEIGHT
CAL_MODULE_CSI_BOTTOM_CORR = CAL_MODULE_CSI_TOP_CORR - 7*CELL_VERT_PITCH -\
                             CSI_HEIGHT


CAL_LAYER_Z_DICT = {0: -58.07,
                    1: -79.42,
                    2: -100.77,
                    3: -122.12,
                    4: -143.47,
                    5: -164.82,
                    6: -186.17,
                    7: -207.52
                    }

def getCalLayerZ(layer):
    return CAL_LAYER_Z_DICT[layer]


X_MAX = 2.1*TOWER_PITCH
X_MIN = -X_MAX
Z_MAX = 0
Z_MIN = -250
MODULE_COLOR = ROOT.kGray
LOG_COLOR = ROOT.kGray


def getSideCanvas(name, title = None, width = 1200, height = 480):
    title = title or name
    rightMargin = 0.01
    leftMargin = 0.05
    c = ROOT.TCanvas(name, title, width, height)
    c.Divide(1, 2)
    c.GetPad(1).SetRightMargin(rightMargin)
    c.GetPad(1).SetLeftMargin(leftMargin)
    c.GetPad(2).SetRightMargin(rightMargin)
    c.GetPad(2).SetLeftMargin(leftMargin)
    return c

def getTopCanvas(name, title = None, width = 600, height = 600):
    title = title or name
    c = ROOT.TCanvas(name, title, width, height)
    c.SetRightMargin(0.01)
    c.SetBottomMargin(0.08)
    c.SetTopMargin(0.075)
    return c



class Box(ROOT.TBox):

    def __init__(self, xcenter, ycenter, width, height, color = None):
        x1 = xcenter - 0.5*width
        y1 = ycenter - 0.5*height
        x2 = xcenter + 0.5*width
        y2 = ycenter + 0.5*height
        self.XCenter = xcenter
        self.YCenter = ycenter
        self.Width   = width
        self.Height  = height
        ROOT.TBox.__init__(self, x1, y1, x2, y2)
        self.SetFillStyle(0)
        if color is not None:
            self.SetLineColor(color)

    def draw(self):
        self.Draw()



class CalSideModule:

    def __init__(self, xcenter, view):
        self.ModuleBox = Box(xcenter, CAL_MODULE_VCENTER, CAL_MODULE_WIDTH,
                             CAL_MODULE_HEIGHT, MODULE_COLOR)
        self.Logs = []
        if view == 'yz':
            for i in range(8):
                yc = CAL_MODULE_CSI_TOP - i*CELL_VERT_PITCH
                if i % 2 != 0:
                    xc = xcenter
                    log = Box(xc, yc, CSI_LENGTH, CSI_HEIGHT, LOG_COLOR)
                    self.Logs.append(log)
                else:
                    for j in range(12):
                        xc = xcenter + (j - 5.5)*CELL_HOR_PITCH
                        log = Box(xc, yc, CSI_WIDTH, CSI_HEIGHT, LOG_COLOR)
                        self.Logs.append(log)
        elif view == 'xz':
            for i in range(8):
                yc = CAL_MODULE_CSI_TOP - i*CELL_VERT_PITCH
                if i % 2 == 0:
                    xc = xcenter
                    log = Box(xc, yc, CSI_LENGTH, CSI_HEIGHT, LOG_COLOR)
                    self.Logs.append(log)
                else:
                    for j in range(12):
                        xc = xcenter + (j - 5.5)*CELL_HOR_PITCH
                        log = Box(xc, yc, CSI_WIDTH, CSI_HEIGHT, LOG_COLOR)
                        self.Logs.append(log)

    def draw(self):
        self.ModuleBox.draw()
        for log in self.Logs:
            log.draw()



class CalTopModule:

    def __init__(self, xcenter, ycenter):
        self.ModuleBox = Box(xcenter, ycenter, CAL_MODULE_WIDTH,
                             CAL_MODULE_WIDTH, MODULE_COLOR)

    def draw(self):
        self.ModuleBox.draw()
        
        

class CalLayout:

    def __init__(self):
        self.BaseHistogramXZ = ROOT.TH1F('hxz', 'XZ view', 100, X_MIN, X_MAX)
        self.BaseHistogramXZ.SetMinimum(Z_MIN)
        self.BaseHistogramXZ.SetMaximum(Z_MAX)
        self.BaseHistogramXZ.GetXaxis().SetTitle('x (mm)')
        self.BaseHistogramXZ.GetYaxis().SetTitle('z (mm)')
        self.BaseHistogramXZ.GetYaxis().SetTitleOffset(0.60)
        self.BaseHistogramYZ = self.BaseHistogramXZ.Clone('hyz')
        self.BaseHistogramYZ.SetTitle('YZ view')
        self.BaseHistogramYZ.GetXaxis().SetTitle('y (mm)')
        self.BaseHistogramXY = self.BaseHistogramXZ.Clone('hxy')
        self.BaseHistogramYZ.SetTitle('XY view')
        self.BaseHistogramXY.SetMinimum(X_MIN)
        self.BaseHistogramXY.SetMaximum(X_MAX)
        self.BaseHistogramXY.GetXaxis().SetTitle('x (mm)')
        self.BaseHistogramXY.GetYaxis().SetTitle('y (mm)')
        self.BaseHistogramXY.GetYaxis().SetTitleOffset(1.5)
        self.Modules = {'xz': [],
                        'yz': [],
                        'xy': []
                        }
        for i in range(4):
            xcenter = (i - 1.5)*TOWER_PITCH
            for view in ['xz', 'yz']:
                self.Modules[view].append(CalSideModule(xcenter, view))
            for j in range(4):
                ycenter = (j - 1.5)*TOWER_PITCH
                self.Modules['xy'].append(CalTopModule(xcenter, ycenter))

    def draw(self, view = 'xz'):
        if view == 'xz':
            self.BaseHistogramXZ.Draw()
            for module in self.Modules['xz']:
                module.draw()
        elif view == 'yz':
            self.BaseHistogramYZ.Draw()
            for module in self.Modules['yz']:
                module.draw()
        elif view == 'xy':
            self.BaseHistogramXY.Draw()
            for module in self.Modules['xy']:
                module.draw()



if __name__ == '__main__':
    cal = CalLayout()
    c1 = getSideCanvas('Side')
    c1.cd(1); cal.draw('xz')
    c1.cd(2); cal.draw('yz')
    c1.cd(); c1.Update()
    c2 = getTopCanvas('Top')
    cal.draw('xy')
    c2.Update()

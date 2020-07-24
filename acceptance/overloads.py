import ROOT
import os


class niceTLegend(ROOT.TLegend):
    def __init__(self, xlow=0.6, ylow=0.6, xhigh=0.85, yhigh=0.85, textsize=0.05):
        super(niceTLegend, self).__init__()
        self.SetX1NDC(xlow)
        self.SetX2NDC(xhigh)
        self.SetY1NDC(ylow)
        self.SetY2NDC(yhigh)
        self.SetTextSize(textsize)
        self.SetMargin(0.2)
        self.SetFillColor(ROOT.TColor.GetColorTransparent(0, 0.2))
        self.SetBorderSize(0)

    def UpdateCoords(self, x1, y1, x2, y2):
        self.Draw()
        ROOT.gPad.Update()
        self.SetX1NDC(x1)
        self.SetY1NDC(y1)
        self.SetX2NDC(x2)
        self.SetY2NDC(y2)


class ExtendedCanvas(ROOT.TCanvas):
    def __init__(self, name="c", title="c", ww=800, wh=600):
        super(ExtendedCanvas, self).__init__(name, title, ww, wh)

    def save(self, name):
        """ Save with name in PDF and PNG format """
        self.update()
        dirname = os.path.dirname(name)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        self.SaveAs(f"{name}.pdf")
        self.SaveAs(f"{name}.png")

    @property
    def empty(self):
        return len(self.GetListOfPrimitives()) == 0


class SinglePadCanvas(ExtendedCanvas):
    _tile = None
    _xtitle = None
    _ytitle = None
    _xaxis = None
    _yaxis = None
    _xrange = None
    _yrange = None

    def __init__(self, name="c", title="c", ww=800, wh=600):
        super(SinglePadCanvas, self).__init__(name, title, ww, wh)
        self.title = ""
        ROOT.gPad.SetTicks()

    def nostats(self):
        """ Remove all stat report boxes """
        for obj in self.GetListOfPrimitives():
            if not isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.TGraph):
                continue
            obj.SetStats(0)

    def update(self):
        """ Update all pads """
        ROOT.gPad.Modified()
        ROOT.gPad.Update()

    def paper_style(self):
        self.SetTopMargin(0.05)
        self.SetRightMargin(0.05)
        self.SetBottomMargin(0.11)
        self.SetLeftMargin(0.11)

        self._find_axes()
        if self._xaxis:
            self._xaxis.SetTitleSize(0.045)
            self._xaxis.SetTitleOffset(1)
            self._xaxis.SetLabelSize(0.045)
        if self._yaxis:
            self._yaxis.SetTitleSize(0.045)
            self._yaxis.SetTitleOffset(1.15)
            self._yaxis.SetLabelSize(0.045)

        self.update()

    @property
    def title(self):
        """ Get current title """
        return self._title

    @title.setter
    def title(self, title):
        """ Set a new title """
        self._title = title
        for obj in self.GetListOfPrimitives():
            if hasattr(obj, "SetTitle"):
                obj.SetTitle(title)
        self.update()

    @property
    def xtitle(self):
        """ Get x axis label """
        return self._xtitle

    @xtitle.setter
    def xtitle(self, label):
        """ Set a new x axis label """
        self._xtitle = label
        if not self._xaxis:
            self._find_axes()
        self._xaxis.SetTitle(label)
        self.update()

    @property
    def ytitle(self):
        """ Get y axis label """
        return self._ytitle

    @ytitle.setter
    def ytitle(self, label):
        """ Set a new x axis label """
        self._ytitle = label
        if not self._yaxis:
            self._find_axes()
        self._yaxis.SetTitle(label)
        self.update()

    def _find_axes(self):
        """ Find all x and y axes painted to the TPads """
        self.update()
        self._xaxis = None
        self._yaxis = None
        for obj in self.GetListOfPrimitives():
            if not (isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.TGraph)):
                continue
            ax = obj.GetXaxis()
            if not isinstance(ax, ROOT.TAxis):
                continue
            self._xaxis = ax
            ax = obj.GetYaxis()
            self._yaxis = ax
            break
        assert self._xaxis is not None, "Could not find x-axes! Are objects drawn?"
        assert self._yaxis is not None, "Could not find y-axes! Are objects drawn?"

    @property
    def xrange(self):
        """ Get the current xrange """
        return self._xrange

    @xrange.setter
    def xrange(self, range_user):
        """ SetRangeUser wrapper """
        self._xrange = range_user
        if not self._xaxis:
            self._find_axes()
        self._xaxis.SetRangeUser(*range_user)
        self.update()

    @property
    def yrange(self):
        """ Get the current xrange """
        return self._xrange

    @yrange.setter
    def yrange(self, range_user):
        """ SetRangeUser wrapper """
        self._yrange = range_user
        if not self._yaxis:
            self._find_axes()
        self._yaxis.SetRangeUser(*range_user)
        self.update()


class noMarginCanvas(SinglePadCanvas):
    def __init__(self, name="c", title="c", ww=800, wh=600):
        super(noMarginCanvas, self).__init__(name, title, ww, wh)
        ROOT.gPad.SetTopMargin(0.045)
        ROOT.gPad.SetRightMargin(0.054)
        ROOT.gPad.SetLeftMargin(0.075)
        ROOT.gPad.SetBottomMargin(0.075)
        ROOT.gPad.SetTicks()

    def presentation_style(self):
        self.SetBottomMargin(0.15)
        self.SetLeftMargin(0.125)
        if self._xaxis:
            self._xaxis.SetTitleSize(0.07)
            self._xaxis.SetTitleOffset(0.8)
        if self._yaxis:
            self._yaxis.SetTitleSize(0.07)
            self._yaxis.SetTitleOffset(0.8)
        self.update()

    def paper_style(self):
        self.SetBottomMargin(0.11)
        self.SetLeftMargin(0.11)
        self._find_axes()
        if self._xaxis:
            self._xaxis.SetTitleSize(0.045)
            self._xaxis.SetTitleOffset(1)
        if self._yaxis:
            self._yaxis.SetTitleSize(0.045)
            self._yaxis.SetTitleOffset(1.15)

        self.update()


class topMarginCanvas(noMarginCanvas):
    def __init__(self, name="c", title="c", ww=800, wh=600):
        super(topMarginCanvas, self).__init__(name, title, ww, wh)
        ROOT.gPad.SetTopMargin(0.09)
        ROOT.gPad.SetTicks()


class optimizedLRMarginCanvas(SinglePadCanvas):
    def __init__(self, name="c", title="c", ww=800, wh=600):
        super(optimizedLRMarginCanvas, self).__init__(name, title, ww, wh)
        ROOT.gPad.SetTopMargin(0.045)
        ROOT.gPad.SetRightMargin(0.1)
        ROOT.gPad.SetLeftMargin(0.1)
        ROOT.gPad.SetBottomMargin(0.1)
        ROOT.gPad.SetTicks()

    def presentation_style(self):
        pass

    def paper_style(self):
        self.SetBottomMargin(0.11)
        self.SetLeftMargin(0.11)
        self.SetRightMargin(0.1)
        self._find_axes()
        if self._xaxis:
            self._xaxis.SetTitleSize(0.045)
            self._yaxis.SetTitleOffset(0.9)
        if self._yaxis:
            self._yaxis.SetTitleSize(0.045)
            self._yaxis.SetTitleOffset(1.15)
        self.update()


class optimizedLBMarginCanvas(ROOT.TCanvas):
    def __init__(self, name="c", title="c", ww=800, wh=600):
        super(optimizedLBMarginCanvas, self).__init__(name, title, ww, wh)
        ROOT.gPad.SetTopMargin(0.045)
        ROOT.gPad.SetRightMargin(0.054)
        ROOT.gPad.SetLeftMargin(0.1)
        ROOT.gPad.SetBottomMargin(0.1)
        ROOT.gPad.SetTicks()


class splitCanvas(ExtendedCanvas):
    def __init__(self, name="c", title="c", ratio="2to1", ww=800, wh=800):
        super(splitCanvas, self).__init__(name, title, ww, wh)
        self._title = title
        self._xrange = []
        self._yranges = []
        self._xtitle = ""
        self._ytitles = ["", ""]
        self.__xaxes = []
        self.__yaxes = []
        self.__logx = False
        ratio_frst = float(ratio.split("to")[0])
        ratio_scnd = float(ratio.split("to")[1])
        padBorder = ratio_frst / (ratio_frst + ratio_scnd)
        padBorder = 1.0 - padBorder
        self.pad_ratio = ratio_frst / ratio_scnd
        self.pad1 = ROOT.TPad("pad1", "upper pad", 0.0, padBorder, 1.0, 1.0)
        self.pad2 = ROOT.TPad("pad2", "lower pad", 0.0, 0.0, 1.0, padBorder)
        name1 = "pad1_at_0x{0:x}".format(id(self.pad1))
        self.pad1.SetName(name1)
        self.pad1.SetTitle(name1)
        name2 = "pad1_at_0x{0:x}".format(id(self.pad2))
        self.pad2.SetName(name2)
        self.pad2.SetTitle(name2)
        for pad in [self.pad1, self.pad2]:
            pad.UseCurrentStyle()
            pad.SetBorderMode(0)
            pad.SetFrameBorderMode(0)
            pad.SetBorderSize(0)
            pad.SetTicks()
            pad.Draw()
            pad.SetRightMargin(0.05)
        self.pad1.SetBottomMargin(0.02)
        self.pad2.SetTopMargin(0.02)
        self.pad2.SetBottomMargin(0.25)

    def update(self):
        """ Update all pads """
        self.pad1.Modified()
        self.pad1.Update()
        self.pad2.Modified()
        self.pad2.Update()
        ROOT.gPad.Modified()
        ROOT.gPad.Update()

    def adjust_size(self):
        """ Adjust the tick and label sizes as required """
        self.update()
        if not self.__xaxes:
            self.__find_axes()
        self.__xaxes[0].SetLabelSize(0)

        self.__xaxes[1].SetTickLength(1.0 / self.pad_ratio * self.__xaxes[0].GetTickLength())
        self.__xaxes[1].SetLabelOffset(0.02)
        self.__xaxes[1].SetLabelSize()  # reset it
        self.__xaxes[1].SetLabelSize(1.0 / self.pad_ratio * self.__xaxes[1].GetLabelSize())
        self.__xaxes[1].SetTitleSize(1.0 / self.pad_ratio * 0.08)
        self.__xaxes[1].SetTitleOffset(0.8)

        self.__yaxes[0].SetTitleSize(1.0 / self.pad_ratio * 0.08)
        self.__yaxes[0].SetTitleOffset(0.5)

        self.__yaxes[1].SetLabelSize(1.0 / self.pad_ratio * self.__yaxes[0].GetLabelSize())
        self.__yaxes[1].SetTitleOffset(0.5)
        self.__yaxes[1].SetTitleSize(1.0 / self.pad_ratio * 0.08)
        self.__yaxes[1].CenterTitle()

        # redraw top axis
        self.pad1.RedrawAxis()
        self.update()

    def __find_axes(self):
        """ Find all x and y axes painted to the TPads """
        self.update()
        self.__xaxes = []
        self.__yaxes = []
        for obj in self.pad1.GetListOfPrimitives():
            if not (isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.TGraph)):
                continue
            ax = obj.GetXaxis()
            if not isinstance(ax, ROOT.TAxis):
                continue
            self.__xaxes.append(ax)
            ax = obj.GetYaxis()
            self.__yaxes.append(ax)
            break
        assert len(self.__xaxes) == 1, f"Wrong number of x-axes: {len(self.__xaxes)}"
        assert len(self.__yaxes) == 1, f"Wrong number of y-axes: {len(self.__yaxes)}"
        for obj in self.pad2.GetListOfPrimitives():
            if not (isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.TGraph)):
                continue
            ax = obj.GetXaxis()
            if not isinstance(ax, ROOT.TAxis):
                continue
            self.__xaxes.append(ax)
            ax = obj.GetYaxis()
            self.__yaxes.append(ax)
            break
        assert len(self.__xaxes) == 2, f"Wrong number of x-axes: {len(self.__xaxes)}"
        assert len(self.__yaxes) == 2, f"Wrong number of y-axes: {len(self.__yaxes)}"

    @property
    def title(self):
        """ Get current title """
        return self._title

    @title.setter
    def title(self, title):
        """ Set a new title """
        self._title = title
        for obj in self.pad1.GetListOfPrimitives():
            if hasattr(obj, "SetTitle"):
                obj.SetTitle(title)
        for obj in self.pad2.GetListOfPrimitives():
            if hasattr(obj, "SetTitle"):
                obj.SetTitle("")
        self.update()

    @property
    def xtitle(self):
        """ Get x axis label """
        return self._xtitle

    @xtitle.setter
    def xtitle(self, label):
        """ Set a new x axis label """
        self._xtitle = label
        if not self.__xaxes:
            self.__find_axes()
        self.__xaxes[1].SetTitle(label)

    @property
    def ytitles(self):
        """ Get the two y axis titles for the upper and lower TPads """
        return self._ytitles

    @ytitles.setter
    def ytitles(self, titles):
        if not (isinstance(titles, list) or isinstance(titles, tuple)):
            raise ValueError("Must specify a tuple/list of titles")
        self._ytitles = titles
        if not self.__yaxes:
            self.__find_axes()
        self.__yaxes[0].SetTitle(titles[0])
        self.__yaxes[1].SetTitle(titles[1])

    @property
    def xrange(self):
        """ Get the current xrange """
        return self._xrange

    @xrange.setter
    def xrange(self, range_user):
        """ SetRangeUser wrapper """
        self._xrange = range_user
        if not self.__xaxes:
            self.__find_axes()
        self.__xaxes[0].SetRangeUser(*range_user)
        self.__xaxes[1].SetRangeUser(*range_user)
        self.update()

    @property
    def yranges(self):
        """ Get the current xrange """
        if not self._yranges:
            if not self.__yaxes:
                self.__find_axes()
            self._yranges = (
                (self.__yaxes[0].GetBinLowEdge(0), self.__yaxes[0].GetBinLowEdge(self.__yaxes[0].GetNbins() + 1),),
                (self.__yaxes[1].GetBinLowEdge(0), self.__yaxes[1].GetBinLowEdge(self.__yaxes[1].GetNbins() + 1),),
            )
        return self._yranges

    @yranges.setter
    def yranges(self, list_of_ranges):
        """ SetRangeUser wrapper """
        if len(list_of_ranges) != 2:
            raise ValueError("Must proved a list/tuple for top and bottom TPad")
        self._yranges = list_of_ranges
        if not self.__yaxes:
            self.__find_axes()
        if len(list_of_ranges[0]) > 1:
            self.__yaxes[0].SetRangeUser(*list_of_ranges[0])
        if len(list_of_ranges[1]) > 1:
            self.__yaxes[1].SetRangeUser(*list_of_ranges[1])
        self.update()

    def nostats(self):
        """ Remove all stat report boxes """
        for obj in self.pad1.GetListOfPrimitives():
            if not isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.TGraph):
                continue
            obj.SetStats(0)
        for obj in self.pad2.GetListOfPrimitives():
            if not isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.TGraph):
                continue
            obj.SetStats(0)

    @property
    def logx(self):
        return self.__logx

    @logx.setter
    def logx(self, val):
        self.__logx = val
        self.SetLogx(val)
        self.pad1.SetLogx(val)
        self.pad2.SetLogx(val)
        if not self.__xaxes:
            self.__find_axes()
        self.__xaxes[0].SetMoreLogLabels()
        self.update()

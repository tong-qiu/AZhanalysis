"""
Plot signal acceptance x efficiency from CxAODReader histograms files
"""

import asyncio
import os
import json
import numpy as np
import re

from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

from rootIO import safeTFile
from overloads import ROOT, SinglePadCanvas, niceTLegend
from helpers import ATLASInternal, EnergyAndLumi, DrawText

ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gROOT.SetBatch(1)

_executor = ThreadPoolExecutor(64)
_lock = Lock()

# some static normalisation factors from auxiliary files
xsections_file = (
    "/ptmp/mpp/hyperion/dev/CxAODFrameWork/CxAODVHreso-2019-09-24/source/CxAODOperations/data/XSections_13TeV.txt"
)
with open(xsections_file) as xs_read:
    XS_RAW = xs_read.readlines()
SIGMA_EFF = {
    line.split()[4]: float(line.split()[1]) * float(line.split()[2]) * float(line.split()[3])
    for line in XS_RAW
    if not line.startswith("#") and len(line.split()) > 4
}


class Info:
    signal_regions = (
        "1tag2pjet",
        "2tag2pjet",
        "1tag1pfat_noaddbjetsr",
        "2tag1pfat_noaddbjetsr",
    )
    __target_histograms = {
        "1tag2pjet": "1tag2pjet/{signal}_1tag2pjet_0ptv_SR_mVH",
        "2tag2pjet": "2tag2pjet/{signal}_2tag2pjet_0ptv_SR_mVH",
        "1tag1pfat_noaddbjetsr": "1tag1pfat0pjet/{signal}_1tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH",
        "2tag1pfat_noaddbjetsr": "2tag1pfat0pjet/{signal}_2tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH",
    }

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.backup_fname = f"{self.name}_{self.signalname}.json"
        self.integrals = {region: 0.0 for region in self.signal_regions}
        self.entries = {region: 0 for region in self.signal_regions}
        self.target_histograms = {
            region: hist.format(signal=self.signalname) for region, hist in self.__target_histograms.items()
        }
        self.valid = False
        self.mass = int(re.search("([1-9][0-9]+)", self.signalname).group(1))
        self.sigma_eff = SIGMA_EFF[self.signalname]

    def load_fetch_dir(self):
        """ load all signal files in the fetch directory of the specified path """
        basepath = os.path.join(self.path, "fetch")
        files = [os.path.join(basepath, f) for f in os.listdir(basepath)]
        files = [f for f in files if os.path.isfile(f)]

        # filter by signal name
        def filt(name):
            for search in ("HVT", "bbA", "ggA"):
                if search in self.signalname and search in os.path.basename(name):
                    return True
                if search == "ggA" and self.signalname.startswith("AZh") and search in os.path.basename(name):
                    return True
            return False

        files = [f for f in files if filt(f)]
        for file in files:
            self.read(file)

    def read(self, file):
        """ read a file and add to self """
        print(f"[{self.signalname}] scanning file {file}")
        with safeTFile(file) as fread:
            for region in self.signal_regions:
                hist = fread.Get(self.target_histograms[region])
                if not isinstance(hist, ROOT.TH1):
                    continue
                self.integrals[region] += hist.Integral()
                self.entries[region] += hist.GetEntries()
        self.valid = np.sum(list(self.integrals.values()) > 0)

    def __str__(self):
        val = f"{self.name} - {self.signalname}\n"
        for region in self.signal_regions:
            val += f"{region}: {self.entries[region]:.0f} entries -> {self.integrals[region]:.2f} integral\n"
        return val.rstrip("\n")

    def save(self):
        with open(self.backup_fname, "w") as fout:
            json.dump(self.__dict__, fout)

    def load_from_file(self):
        if not os.path.isfile(self.backup_fname):
            return False
        with open(self.backup_fname) as fread:
            self.__dict__.update(json.load(fread))
            self.valid = np.sum(list(self.entries.values())) > 0  # having at least one entry makes info valid
        return True

    def calculate_efficiency(self, region=None):
        if region is None:
            # sum over regions
            integral = np.sum(np.fromiter(self.integrals.values(), float))
            efficiency = integral / self.lumi / self.sigma_eff
        else:
            efficiency = float(self.integrals[region]) / self.lumi / self.sigma_eff
        return efficiency


class MergedInfo:
    """ Merges MC16 a+d+e into a single info """

    def __init__(self, infos):
        values = list(infos.values())
        self.signalname = values[0].signalname
        self.signal_regions = values[0].signal_regions
        self.mass = values[0].mass
        self.sigma_eff = SIGMA_EFF[self.signalname]

        self.lumi = np.sum(list(info.lumi for info in values))
        self.integrals = {
            region: np.sum(np.fromiter([info.integrals[region] for info in values], float))
            for region in self.signal_regions
        }
        self.entries = {
            region: np.sum(np.fromiter([info.entries[region] for info in values], float))
            for region in self.signal_regions
        }

    def calculate_efficiency(self, region=None):
        if region is None:
            # sum over regions
            integral = np.sum(np.fromiter(self.integrals.values(), float))
            efficiency = integral / self.lumi / self.sigma_eff
        else:
            efficiency = float(self.integrals[region]) / self.lumi / self.sigma_eff
        return efficiency


class InfoPlotter:
    legend_title = {
        "1tag2pjet": "1 #font[52]{b}-tag, resolved",
        "2tag2pjet": "2 #font[52]{b}-tags, resolved",
        "1tag1pfat_noaddbjetsr": "1 #font[52]{b}-tag, merged",
        "2tag1pfat_noaddbjetsr": "2 #font[52]{b}-tags, merged",
        "sum": "All signal regions",
    }
    graph_style = {
        "1tag2pjet": {"color": ROOT.kTeal - 8, "markerstyle": 8, "markersize": 0.75, "linestyle": 1, "linewidth": 2},
        "2tag2pjet": {"color": ROOT.kTeal - 8, "markerstyle": 8, "markersize": 0.75, "linestyle": 2, "linewidth": 2},
        "1tag1pfat_noaddbjetsr": {
            "color": ROOT.kAzure - 8,
            "markerstyle": 8,
            "markersize": 0.75,
            "linestyle": 1,
            "linewidth": 2,
        },
        "2tag1pfat_noaddbjetsr": {
            "color": ROOT.kAzure - 8,
            "markerstyle": 8,
            "markersize": 0.75,
            "linestyle": 2,
            "linewidth": 2,
        },
        "sum": {"color": ROOT.kBlack, "markerstyle": 8, "markersize": 1.2, "linestyle": 1, "linewidth": 3},
    }
    text = {
        "HVTZHvvqq": {"title": "0L channel, HVT Z'#rightarrow Zh", "xtitle": "#font[52]{Z'} resonance mass [GeV]"},
        "AZhvvbb": {
            "title": "0L channel, gg #rightarrow A #rightarrow Zh",
            "xtitle": "#font[52]{A} resonance mass [GeV]",
        },
        "bbAZhvvbb": {
            "title": "0L channel, gg #rightarrow b#bar{b}A, A #rightarrow Zh",
            "xtitle": "#font[52]{A} resonance mass [GeV]",
        },
    }

    def __init__(self, infos):
        infos = filter(lambda x: x and x.signalname != "AZhvvbb700", infos)  # buggy sample
        self.infos = sorted(infos, key=lambda x: x.mass)
        assert len(self.infos) > 0, "Trying to initialise InfoPlotter with 0 infos"
        self.signal_regions = list(self.infos[0].signal_regions)
        self.signalname = self.infos[0].signalname.replace(str(self.infos[0].mass), "")
        self.canv = SinglePadCanvas()
        self.__root_objects = []  # stupid ROOT deletes stuff otherwise

    def apply_style(self, graph, region):
        graph.SetLineStyle(self.graph_style[region]["linestyle"])
        graph.SetLineWidth(self.graph_style[region]["linewidth"])
        graph.SetLineColor(self.graph_style[region]["color"])

        graph.SetMarkerStyle(self.graph_style[region]["markerstyle"])
        graph.SetMarkerSize(self.graph_style[region]["markersize"])
        graph.SetMarkerColor(self.graph_style[region]["color"])

    def plot_by_channel(self):
        masses = [info.mass for info in self.infos]
        efficiencies = {
            region: [info.calculate_efficiency(region) for info in self.infos] for region in self.signal_regions
        }
        efficiencies["sum"] = [info.calculate_efficiency() for info in self.infos]

        self.canv.cd()
        xax, yax = None, None
        legend = niceTLegend()
        self.__root_objects.append(legend)
        for region in ["sum"] + self.signal_regions:
            graph = ROOT.TGraph()
            self.__root_objects.append(graph)
            graph.SetTitle(region)
            for idx, (mass, eff) in enumerate(zip(masses, efficiencies[region])):
                graph.SetPoint(idx, mass, eff)
            if xax is None:
                graph.Draw("APL")
                xax = graph.GetXaxis()
                yax = graph.GetYaxis()
                graph.SetMinimum(0)
                graph.SetMaximum(1.8 * np.max(efficiencies["sum"]))
            else:
                graph.Draw("PLsame")
            self.apply_style(graph, region)
            legend.AddEntry(graph, self.legend_title[region], "PL")

        self.canv.paper_style()
        self.canv.xtitle = self.text[self.signalname]["xtitle"]
        self.canv.ytitle = "Acceptance #times efficiency"
        self.canv.title = ""
        self.canv.xrange = (100, 5100)
        if "HVT" not in self.signalname:
            self.canv.xrange = (250, 1650)

        yax.SetNdivisions(505)

        legend.UpdateCoords(0.6, 0.65, 0.9, 0.93)
        legend.Draw()

        ATLASInternal(x=0.15, y=0.88)
        EnergyAndLumi(x=0.15, y=0.82, lumi=139, size=0.035)
        DrawText(self.text[self.signalname]["title"], x=0.15, y=0.78, size=0.035)

        self.canv.save(f"signal_efficiency_{self.signalname}")


def main(args):
    pattern = re.compile(f"[0-9]+.*({args.signalname}[0-9]+)")
    matches = filter(None, (pattern.search(line) for line in XS_RAW))
    signals = list(set(match.group(1) for match in matches))
    if args.run_async:
        loop = asyncio.get_event_loop()
        tasks = [loop.run_in_executor(_executor, calc, args, signal) for signal in signals]
        infos = loop.run_until_complete(asyncio.gather(*tasks))
    else:
        infos = [calc(args, signal) for signal in signals]

    plotter = InfoPlotter(infos)
    plotter.plot_by_channel()


def calc(args, signal):
    infos = {
        "a": Info(name="mc16a", path=args.signals[0], signalname=signal, lumi=36.2077e3),
        "d": Info(name="mc16d", path=args.signals[1], signalname=signal, lumi=44.3074e3),
        "e": Info(name="mc16e", path=args.signals[2], signalname=signal, lumi=58.4500e3),
    }
    for info in infos.values():
        if not info.load_from_file():
            info.load_fetch_dir()
            info.save()
        if not info.valid:
            return

    return MergedInfo(infos)
    # eff = np.sum(np.fromiter((info.calculate_efficiency() * info.lumi for info in infos.values()), float))
    # total_lumi = np.sum(np.fromiter((info.lumi for info in infos.values()), float))
    # eff /= total_lumi

    # print(f"Total efficiency for {signal}: {eff:1.4e}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--signals", required=True, nargs="+", help="CxAODReader run directory for signal histograms",
    )
    parser.add_argument("--signalname", default="HVTZHvvqq", help="Signal to be investigated")
    parser.add_argument(
        "--run-async", action="store_true", help="Run async, faster if reading new info from ROOT files",
    )
    main(parser.parse_args())

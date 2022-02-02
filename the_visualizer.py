import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
from . import f_translation as trt
from . import adjacent_biopy as adj
from dna_features_viewer import GraphicFeature, GraphicRecord
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLineEdit, QInputDialog, QFormLayout
from PyQt5.QtWidgets import QMainWindow, QLabel, QGridLayout, QDesktopWidget
from PyQt5.QtGui import *
from PyQt5.QtCore import *


def find_nth(string, substring, n):
    """ Finds the nth occurrence of substring in string and returns its index. """
    start = string.find(substring)
    while start >= 0 and n > 1:
        start = string.find(substring, start + len(substring))
        n -= 1
    return start

def get_info(filetype):
    df = pd.read_csv("%s_results_with_rfs.ods" % filetype, sep='\t')
    accession = df["ORF Number"].tolist()
    orf_list = []
    orf_name = []
    for orfises in range(len(accession)):
        # if accession[orfises] not in orf_name:
        #     orf_name.append(accession[orfises])
        if str(accession[orfises]) not in orf_list:
            orf_list.append(str(accession[orfises]))
    orf_seqs = df["ORF Sequence"].tolist()
    peptides = df["pepSeq"].tolist()
    specs = df["Total Spec Counts"].tolist()
    coords = df["Genome Coordinates"].tolist()
    orf_frames = df["Reading Frame"].tolist()
    spec_dic = {}
    for i in range(len(orf_seqs)):
        if orf_seqs[i] not in spec_dic:
            spec_dic[orf_seqs[i]] = ""
            for j in range(len(peptides)):
                if peptides[j] in orf_seqs[i]:
                    spec_dic[orf_seqs[i]] += "%s," % peptides[j]
    # print(spec_dic)
    orfs = []
    for i in range(len(orf_seqs)):
        if orf_seqs[i] not in orfs:
            orfs.append(orf_seqs[i])
    return orf_seqs, coords, spec_dic, orf_list, orf_frames, accession


class ORF(object):
    def __init__(self, frame, sequence, number):
        self.frame = frame
        self.seq = sequence
        self.number = number
        pass


class SpecCount(object):
    def __init__(self):
        pass


def visualization(orf_numbers, filetype):
    orf_seqs, coords, spec_dic, orf_list, orf_frames, accession = get_info(filetype)
    orf_number = accession.index(orf_numbers)
    # sequence = orf_seqs[orf_number]
    # orf_frame = orf_frames[orf_number]
    # orf_coord = coords[orf_number]
    genome = "MSMEG_genome.fasta"
    # orf_x = trt.ReadingFrame(genome, str(sequence), orf_coord)
    # rf = orf_x.find_frame()
    orf = ORF(orf_frames[orf_number], orf_seqs[orf_number], orf_number)

    # print(sequence)
    peps = spec_dic[orf.seq].split(",")[:-1]
    # print(peps)
    # print(sequence)
    # coord = coords[orf_number].split(" - ")
    orf.coords = coords[orf_number].split(" - ")
    # print(coord)
    if int(orf.coords[0]) > int(orf.coords[1]):
        start_c = int(orf.coords[1])
        end_c = int(orf.coords[0])
        orf.orientation = -1
    elif int(orf.coords[0]) < int(orf.coords[1]):
        start_c = int(orf.coords[0])
        end_c = int(orf.coords[1])
        orf.orientation = +1

    # finding adjacent genes
    gb_info = adj.Adjacent("RefSeq_Reading_Frames.ods", start_c, end_c)
    gb_features = gb_info.find_neighbours()
    if len(gb_features) > 0:
        gb_starts = []
        gb_real_starts= []
        gb_real_ends = []
        gb_ends = []
        gb_names = []
        gb_seqs = []
        gb_frames = []
        gb_index = 0
        while gb_index < len(gb_features):
            gb_starts.append(gb_features[gb_index])
            gb_ends.append(gb_features[gb_index + 1])
            if gb_features[gb_index] < gb_features[gb_index+1]:
                gb_real_starts.append(gb_features[gb_index])
                gb_real_ends.append(gb_features[gb_index+1])
            elif gb_features[gb_index] > gb_features[gb_index+1]:
                gb_real_starts.append(gb_features[gb_index+1])
                gb_real_ends.append(gb_features[gb_index])
            gb_names.append(gb_features[gb_index+2][2:-2])
            gb_seqs.append(gb_features[gb_index+3][2:-2])
            gb_frames.append(gb_features[gb_index+4])
            gb_index += 5
    # print(gb_starts)
    # print(gb_features)


    features = [
        GraphicFeature(start=start_c, end=end_c, strand=orf.orientation, color="#ffd700", label="ORF %s | %s" %
                                                                                            (orf_numbers, orf.frame)),
    ]

    features[0].data["fixed_level"] = int(orf.frame[2:])
    # features[0].strand = 1
    # adding spec counts to figure
    spec_number = 0
    pep_dic = {

    }
    for m in range(len(peps)):
        pep_start = orf.seq.find(peps[m])*3
        # print(pep_start)
        pep_end = (len(peps[m])*3+pep_start)
        if peps[m] not in pep_dic:
            if int(orf.coords[0]) > int(orf.coords[1]):
                features.append(
                    GraphicFeature(start=start_c + pep_start, end=start_c + pep_end, color="#ccccff", thickness=3))
            if int(orf.coords[0]) < int(orf.coords[1]):
                features.append(
                    GraphicFeature(start=start_c + pep_end, end=start_c + pep_start, color="#ccccff", thickness=3))
            pep_dic[peps[m]] = 1
            spec_number += 1
            features[0 + spec_number].data["fixed_level"] = int(orf.frame[2:]) + spec_number
        elif peps[m] in pep_dic:
            if pep_dic[peps[m]] <= 10:
                if int(orf.coords[0]) > int(orf.coords[1]):
                    features.append(
                        GraphicFeature(start=start_c + pep_start, end=start_c + pep_end, color="#ccccff", thickness=3))
                if int(orf.coords[0]) < int(orf.coords[1]):
                    features.append(
                        GraphicFeature(start=start_c + pep_end, end=start_c + pep_start, color="#ccccff", thickness=3))
                pep_dic[peps[m]] += 1
                spec_number += 1
                features[0 + spec_number].data["fixed_level"] = int(orf.frame[2:]) + spec_number
    pep_dic = {}


    # adding adjacent genes to figure
    if len(gb_features) > 0:
        for k in range(len(gb_starts)):
            orient = ""
            if int(gb_starts[k]) < int(gb_ends[k]):
                gb_coords = "%s - %s" % (gb_starts[k], gb_ends[k])
                orient = (+1)
            elif int(gb_starts[k]) > int(gb_ends[k]):
                gb_coords = "%s - %s" % (gb_ends[k], gb_starts[k])
                orient = (-1)
            # gb_coords = "%s - %s" % (gb_starts[k], gb_ends[k])
            # print(gb_seqs[k])
            # gb_x = trt.ReadingFrame(genome, str(gb_seqs[k]), gb_coords)
            # gb_rf = gb_x.find_frame()
            features.append(GraphicFeature(start=int(gb_starts[k]), end=int(gb_ends[k]), color="lightblue",
                                           label="%s | %s" % (gb_names[k], gb_frames[k]), strand=orient))

    # features[3].open_left
    # features[3].strand = 1
    # print(vars(features[3]))
    # print(gb_names)
    for i in range(len(features)):
        if features[i].label is not None:
            finder = find_nth(features[i].label, " |", 1)
            # print(features[i].label)
            # print(features[i].label[:finder])
            # print(features[i].label[:finder])
            # if "fixed_level" in features[i].data:
                # print("RF:" + str(features[i].data['fixed_level']))
                #
            if len(gb_features) > 0:
                if features[i].label[:finder] in gb_names:
                    # print("label: %s" % features[i].label[:finder])
                    ind = gb_names.index(features[i].label[:finder])
                    # print(gb_names[ind])
                    features[i].data['fixed_level'] = int(gb_frames[ind][2:])
                # print(features[i].data['fixed_level'])
        # if features[i].strand == 1:
        #     features[i].open_left
        # elif features[i].strand == -1:
        #     features[i].open_right
            print(features[i].label)
            if "RF -" in features[i].label:
                print(True)
                features[i].strand = 1
            # elif "RF +" in features[i].label:
            #     features[i].strand = -1
        print(features[i].strand)
    if "ORF" in features[0].label and "RF -" in features[0].label:
        features[0].strand = -1
        # features[i].strand
    print(vars(orf))
    # real stuff

    seq_len = end_c-start_c
    seq_lim_1 = start_c-30
    seq_lim_2 = end_c+30
    if len(gb_features) > 0:
        sorted_gb_starts = sorted(gb_starts)
        sorted_gb_ends = sorted(gb_ends)
        if start_c > int(sorted_gb_starts[0]) and end_c < int(sorted_gb_ends[int(len(sorted_gb_ends))-1]):
            seq_lim_1 = int(sorted_gb_starts[0])
            seq_lim_2 = int(sorted_gb_ends[int(len(sorted_gb_ends))-1])
        elif start_c > int(sorted_gb_starts[0]) and end_c > int(sorted_gb_ends[int(len(sorted_gb_ends))-1]):
            seq_lim_1 = int(sorted_gb_starts[0])
            seq_lim_2 = end_c
        elif start_c < int(sorted_gb_starts[0]) and end_c > int(sorted_gb_ends[int(len(sorted_gb_ends))-1]):
            seq_lim_1 = start_c
            seq_lim_2 = end_c
        elif start_c < int(sorted_gb_starts[0]) and end_c < int(sorted_gb_ends[int(len(sorted_gb_ends))-1]):
            seq_lim_1 = start_c
            seq_lim_2 = int(sorted_gb_ends[int(len(sorted_gb_ends))-1])

    #     seq_len = int(gb_ends[len(gb_ends)-1])-(int(gb_starts[0]))




    # print(gb_ends)
    # get seq limits
    feature_starts = []
    feature_ends = []
    for i in range(len(features)):
        feature_starts.append(features[i].start)
        feature_ends.append(features[i].end)
    # print(feature_starts)
    # print(feature_ends)
    lowest = feature_starts[0]
    highest = feature_ends[0]
    for i in range(len(feature_starts)):
        if int(feature_starts[i]) < lowest:
            lowest = feature_starts[i]
        if int(feature_ends[i]) > highest:
            highest = feature_ends[i]
    if len(gb_features) > 0:
        for i in range(len(gb_real_starts)):
            if int(gb_real_starts[i]) < lowest:
                lowest = int(gb_real_starts[i])
            if int(gb_real_ends[i]) > highest:
                highest = int(gb_real_ends[i])
    print(lowest, highest)
    # print(orf_for_tick)

    # GraphicFeature()
    orf_for_tick = [lowest]
    for aa in range(lowest, highest):
        if aa % 200 == 0:
            orf_for_tick.append(aa)
    if highest not in orf_for_tick:
        orf_for_tick.append(highest)
    record = GraphicRecord(sequence=orf.seq, sequence_length=seq_len, features=features)
    ax, _ = record.plot(figure_width=15)

    cooords = []
    cooords.append(orf_for_tick[0])
    cooords.append(orf_for_tick[len(orf_for_tick)-1])
    plt.xticks(np.arange(lowest, highest, step=200), orf_for_tick)

    plt.xlim([lowest-10, highest+10])
    plt.ylim([-5, spec_number+2+abs(int(orf.frame[3:]))])
    print(len(features))
    limit = len(features)
    # plt.ylim(-limit, limit)
    # plt.xlim([gb_starts[0], gb_ends[len(gb_ends)-1]])
    # plt.plot()



    plt.show()


class Browser(QWidget):
    def __init__(self, filetype, parent=None):
        super(Browser, self).__init__(parent)
        self.filetype = filetype
        layout = QFormLayout()
        self.btn = QPushButton("Choose from list")
        self.btn.clicked.connect(self.getItem)

        self.le = QLineEdit()
        layout.addRow(self.btn, self.le)

        self.setLayout(layout)
        self.setWindowTitle("Input Dialog demo")

        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        # qtRectangle = self.frameGeometry()

        # button = QPushButton('Find ORF', self)
        # button.setToolTip('Click to locate ORF along with its spec counts in the genome')
        # button.move(100, 70)
        # button.clicked.connect()
        # self.show()

    def getItem(self):
        orf_seqs, coords, spec_dic, items, orf_frames, accession = get_info(self.filetype)

        item, ok = QInputDialog.getItem(self, "select input dialog",
                                        "list of ORFs", items, 0, False)

        if ok and item:
            visualization(int(item), self.filetype)
    # @pyqtSlot()
    # def on_click(self):
    #     visualization()


def prog_browser(arg):
    app = QApplication(sys.argv)
    ex = Browser(arg)
    ex.show()
    sys.exit(app.exec_())

# visualization(115)
# prog_browser()


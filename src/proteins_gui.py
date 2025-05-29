# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
#
# This file is part of uProteInS.
#
# uProteInS is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# uProteInS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# uProteInS. If not, see <https://www.gnu.org/licenses/>.


import gui_args as gargs
import sys
import os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *

class Ui_ProteInS(object):
    def __init__(self):
        self.files = {}

    def setupUi(self, ProteInS):
        ProteInS.setObjectName("ProteInS")
        ProteInS.resize(800, 601)
        self.centralwidget = QtWidgets.QWidget(ProteInS)
        self.centralwidget.setObjectName("centralwidget")
        self.Modes = QtWidgets.QTabWidget(self.centralwidget)
        self.Modes.setGeometry(QtCore.QRect(0, 0, 801, 561))
        self.Modes.setObjectName("Modes")
        self.Assembly = QtWidgets.QWidget()
        self.Assembly.setObjectName("Assembly")
        self.read1 = QtWidgets.QLabel(self.Assembly)
        self.read1.setGeometry(QtCore.QRect(60, 50, 71, 16))
        self.read1.setObjectName("read1")

        self.read2 = QtWidgets.QLabel(self.Assembly)
        self.read2.setGeometry(QtCore.QRect(60, 90, 81, 16))
        self.read2.setObjectName("read2")

        # find read 1
        self.find_read1 = QtWidgets.QPushButton(self.Assembly)
        self.find_read1.setGeometry(QtCore.QRect(320, 50, 31, 21))
        self.find_read1.setObjectName("find_read1")
        self.find_read1.clicked.connect(lambda: self.find_files("rna_left"))

        # find read 2
        self.find_read2 = QtWidgets.QPushButton(self.Assembly)
        self.find_read2.setGeometry(QtCore.QRect(320, 90, 31, 21))
        self.find_read2.setObjectName("find_read2")
        self.find_read2.clicked.connect(lambda: self.find_files("rna_right"))


        self.specifics = QtWidgets.QGroupBox(self.Assembly)
        self.specifics.setGeometry(QtCore.QRect(440, 10, 341, 141))
        self.specifics.setObjectName("specifics")

        # cpu type
        self.rna_cpu_type = QtWidgets.QLineEdit(self.specifics)
        self.rna_cpu_type.setGeometry(QtCore.QRect(10, 61, 21, 20))
        self.rna_cpu_type.setObjectName("rna_cpu_type")

        self.rna_cpu = QtWidgets.QLabel(self.specifics)
        self.rna_cpu.setGeometry(QtCore.QRect(40, 57, 61, 31))
        self.rna_cpu.setObjectName("rna_cpu")
        self.k_mer = QtWidgets.QLabel(self.specifics)
        self.k_mer.setGeometry(QtCore.QRect(40, 19, 81, 31))
        self.k_mer.setObjectName("k_mer")
        self.k_mer_type = QtWidgets.QLineEdit(self.specifics)
        self.k_mer_type.setGeometry(QtCore.QRect(10, 23, 21, 20))
        self.k_mer_type.setObjectName("k_mer_type")

        # rna memory
        self.rna_memory = QtWidgets.QLabel(self.specifics)
        self.rna_memory.setGeometry(QtCore.QRect(50, 92, 81, 31))
        self.rna_memory.setObjectName("rna_memory")
        self.rna_memory_type = QtWidgets.QLineEdit(self.specifics)
        self.rna_memory_type.setGeometry(QtCore.QRect(10, 99, 31, 20))
        self.rna_memory_type.setObjectName("rna_memory_type")


        self.reads_box = QtWidgets.QGroupBox(self.Assembly)
        self.reads_box.setGeometry(QtCore.QRect(10, 10, 411, 221))
        self.reads_box.setObjectName("reads_box")
        self.paired_reads_box = QtWidgets.QGroupBox(self.reads_box)
        self.paired_reads_box.setGeometry(QtCore.QRect(20, 20, 371, 101))
        self.paired_reads_box.setObjectName("paired_reads_box")
        self.input_left_read = QtWidgets.QLabel(self.paired_reads_box)
        self.input_left_read.setGeometry(QtCore.QRect(100, 20, 47, 13))
        self.input_left_read.setObjectName("input_left_read")
        self.single_end_box = QtWidgets.QGroupBox(self.reads_box)
        self.single_end_box.setGeometry(QtCore.QRect(20, 130, 371, 80))
        self.single_end_box.setObjectName("single_end_box")

        # single-end reads find
        self.single_reads_find = QtWidgets.QPushButton(self.single_end_box)
        self.single_reads_find.setGeometry(QtCore.QRect(290, 40, 31, 21))
        self.single_reads_find.setObjectName("single_reads_find")
        self.single_reads_find.clicked.connect(lambda: self.find_files("rna_single"))


        self.single_reads = QtWidgets.QLabel(self.single_end_box)
        self.single_reads.setGeometry(QtCore.QRect(30, 40, 71, 16))
        self.single_reads.setObjectName("single_reads")
        self.rna_output_box = QtWidgets.QGroupBox(self.Assembly)
        self.rna_output_box.setGeometry(QtCore.QRect(10, 240, 411, 81))
        self.rna_output_box.setObjectName("rna_output_box")

        # rna output directory find
        self.rna_outdir_find = QtWidgets.QPushButton(self.rna_output_box)
        self.rna_outdir_find.setGeometry(QtCore.QRect(360, 30, 31, 21))
        self.rna_outdir_find.setObjectName("rna_outdir_find")
        self.rna_outdir_find.clicked.connect(lambda: self.find_directory("outdir", "YES"))


        self.rna_output_dir = QtWidgets.QLabel(self.rna_output_box)
        self.rna_output_dir.setGeometry(QtCore.QRect(19, 30, 91, 16))
        self.rna_output_dir.setObjectName("rna_output_dir")
        self.rna_run_proteins = QtWidgets.QToolButton(self.Assembly)
        self.rna_run_proteins.setGeometry(QtCore.QRect(440, 160, 341, 31))
        self.rna_run_proteins.setObjectName("rna_run_proteins")
        self.rna_run_proteins.clicked.connect(lambda: self.set_assembly_args())
        self.rna_output_box.raise_()
        self.reads_box.raise_()
        self.specifics.raise_()
        self.read1.raise_()
        self.read2.raise_()
        self.find_read1.raise_()
        self.find_read2.raise_()
        self.rna_run_proteins.raise_()
        self.Modes.addTab(self.Assembly, "")

        """ Database generator """
        self.Database = QtWidgets.QWidget()
        self.Database.setObjectName("Database")
        self.db_files = QtWidgets.QGroupBox(self.Database)
        self.db_files.setGeometry(QtCore.QRect(10, 10, 401, 101))
        self.db_files.setObjectName("db_files")

        ## Genome find ##
        self.db_genome_find = QtWidgets.QPushButton(self.db_files)
        self.db_genome_find.setGeometry(QtCore.QRect(340, 20, 31, 21))
        self.db_genome_find.setObjectName("db_genome_find")
        self.db_genome_find.clicked.connect(lambda: self.find_file("genome"))
        self.db_genome = QtWidgets.QLabel(self.db_files)
        self.db_genome.setGeometry(QtCore.QRect(20, 20, 81, 16))
        self.db_genome.setObjectName("db_genome")

        ## Proteome find ##
        self.db_proteome_find = QtWidgets.QPushButton(self.db_files)
        self.db_proteome_find.setGeometry(QtCore.QRect(340, 60, 31, 21))
        self.db_proteome_find.setObjectName("db_proteome_find")
        self.db_proteome_find.clicked.connect(lambda: self.find_file("proteome"))
        self.db_proteome = QtWidgets.QLabel(self.db_files)
        self.db_proteome.setGeometry(QtCore.QRect(20, 60, 91, 16))
        self.db_proteome.setObjectName("db_proteome")


        ## Database specifications ##
        self.db_specifics = QtWidgets.QGroupBox(self.Database)
        self.db_specifics.setGeometry(QtCore.QRect(430, 10, 351, 191))
        self.db_specifics.setObjectName("db_specifics")

        ## DB Start codons ##
        self.db_starts = QtWidgets.QLabel(self.db_specifics)
        self.db_starts.setGeometry(QtCore.QRect(130, 30, 71, 16))
        self.db_starts.setObjectName("db_starts")
        self.db_starts_type = QtWidgets.QLineEdit(self.db_specifics)
        self.db_starts_type.setGeometry(QtCore.QRect(10, 30, 111, 20))
        self.db_starts_type.setObjectName("db_starts_type")

        ## DB Stop codons ##
        self.db_ends = QtWidgets.QLabel(self.db_specifics)
        self.db_ends.setGeometry(QtCore.QRect(130, 64, 71, 16))
        self.db_ends.setObjectName("db_ends")
        self.db_ends_type = QtWidgets.QLineEdit(self.db_specifics)
        self.db_ends_type.setGeometry(QtCore.QRect(10, 64, 111, 20))
        self.db_ends_type.setObjectName("db_ends_type")

        ## DB Max ORF size ##
        self.db_maxsize = QtWidgets.QLabel(self.db_specifics)
        self.db_maxsize.setGeometry(QtCore.QRect(60, 100, 71, 16))
        self.db_maxsize.setObjectName("db_maxsize")
        self.db_maxsize_type = QtWidgets.QLineEdit(self.db_specifics)
        self.db_maxsize_type.setGeometry(QtCore.QRect(10, 100, 41, 20))
        self.db_maxsize_type.setObjectName("db_maxsize_type")

        ## DB Min ORF size ##
        self.db_minsize_type = QtWidgets.QLineEdit(self.db_specifics)
        self.db_minsize_type.setGeometry(QtCore.QRect(10, 130, 41, 20))
        self.db_minsize_type.setObjectName("db_minsize_type")
        self.db_minsize = QtWidgets.QLabel(self.db_specifics)
        self.db_minsize.setGeometry(QtCore.QRect(60, 130, 71, 16))
        self.db_minsize.setObjectName("db_minsize")

        ## DB Transcriptome check ##
        self.db_transcriptome = QtWidgets.QCheckBox(self.db_specifics)
        self.db_transcriptome.setGeometry(QtCore.QRect(10, 160, 91, 21))
        self.db_transcriptome.setObjectName("db_transcriptome")

        ## DB Run ProteInS ##
        self.db_run_proteins = QtWidgets.QToolButton(self.Database)
        self.db_run_proteins.setGeometry(QtCore.QRect(430, 210, 351, 31))
        self.db_run_proteins.setObjectName("db_run_proteins")
        self.db_run_proteins.clicked.connect(lambda: self.set_db_args())


        ## DB Output directory ##
        self.db_output_box = QtWidgets.QGroupBox(self.Database)
        self.db_output_box.setGeometry(QtCore.QRect(10, 120, 401, 81))
        self.db_output_box.setObjectName("db_output_box")
        self.db_outdir_find = QtWidgets.QPushButton(self.db_output_box)
        self.db_outdir_find.setGeometry(QtCore.QRect(340, 33, 31, 21))
        self.db_outdir_find.setObjectName("db_outdir_find")
        self.db_outdir = QtWidgets.QLabel(self.db_output_box)
        self.db_outdir.setGeometry(QtCore.QRect(10, 33, 101, 16))
        self.db_outdir.setObjectName("db_outdir")
        self.db_outdir_find.clicked.connect(lambda: self.find_directory("outdir", "YES"))


        self.Modes.addTab(self.Database, "")


        """ Peptide Search Mode """
        self.PepSearch = QtWidgets.QWidget()
        self.PepSearch.setObjectName("PepSearch")
        self.pep_specifications_area = QtWidgets.QScrollArea(self.PepSearch)
        self.pep_specifications_area.setGeometry(QtCore.QRect(410, 20, 371, 441))
        self.pep_specifications_area.setWidgetResizable(True)
        self.pep_specifications_area.setObjectName("pep_specifications_area")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 352, 545))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.formLayout = QtWidgets.QFormLayout(self.scrollAreaWidgetContents)
        self.formLayout.setObjectName("formLayout")

        ## Enzyme ##
        self.pep_enzyme_type = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.pep_enzyme_type.setMaximum(9)
        self.pep_enzyme_type.setObjectName("pep_enzyme_type")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.pep_enzyme_type)
        self.pep_enzyme = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_enzyme.setObjectName("pep_enzyme")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.pep_enzyme)

        ## Modifications ##
        self.pep_modifications_type = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.pep_modifications_type.setMinimum(0)
        self.pep_modifications_type.setMaximum(4)
        self.pep_modifications_type.setObjectName("pep_modifications_type")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.pep_modifications_type)
        self.pep_modifications = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_modifications.setObjectName("pep_modifications")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.pep_modifications)
        self.pep_fragmentation_type = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.pep_fragmentation_type.setMinimum(1)
        self.pep_fragmentation_type.setMaximum(4)
        self.pep_fragmentation_type.setObjectName("pep_fragmentation_type")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.pep_fragmentation_type)
        self.pep_fragmentation = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_fragmentation.setObjectName("pep_fragmentation")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.pep_fragmentation)
        self.pep_protocol_type = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.pep_protocol_type.setMaximum(5)
        self.pep_protocol_type.setObjectName("pep_protocol_type")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.pep_protocol_type)
        self.pep_protocol = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_protocol.setObjectName("pep_protocol")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.pep_protocol)
        self.pep_instrument_type = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.pep_instrument_type.setMaximum(3)
        self.pep_instrument_type.setObjectName("pep_instrument_type")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.pep_instrument_type)
        self.pep_instrument = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_instrument.setObjectName("pep_instrument")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.pep_instrument)
        self.spinBox = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.spinBox.setMaximum(2)
        self.spinBox.setObjectName("spinBox")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.spinBox)
        self.pep_instrument_2 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_instrument_2.setObjectName("pep_instrument_2")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.pep_instrument_2)
        self.pep_tolerance_type = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_tolerance_type.setObjectName("pep_tolerance_type")
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.pep_tolerance_type)
        self.pep_tolerance = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_tolerance.setObjectName("pep_tolerance")
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.pep_tolerance)
        self.pep_decoy_type = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_decoy_type.setObjectName("pep_decoy_type")
        self.formLayout.setWidget(10, QtWidgets.QFormLayout.LabelRole, self.pep_decoy_type)
        self.pep_decoy = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_decoy.setObjectName("pep_decoy")
        self.formLayout.setWidget(10, QtWidgets.QFormLayout.FieldRole, self.pep_decoy)
        self.pep_memory_type = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_memory_type.setObjectName("pep_memory_type")
        self.formLayout.setWidget(12, QtWidgets.QFormLayout.LabelRole, self.pep_memory_type)
        self.pep_memory = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory.setObjectName("pep_memory")
        self.formLayout.setWidget(12, QtWidgets.QFormLayout.FieldRole, self.pep_memory)
        self.pep_memory_type_2 = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_memory_type_2.setObjectName("pep_memory_type_2")
        self.formLayout.setWidget(14, QtWidgets.QFormLayout.LabelRole, self.pep_memory_type_2)
        self.pep_memory_2 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory_2.setObjectName("pep_memory_2")
        self.formLayout.setWidget(14, QtWidgets.QFormLayout.FieldRole, self.pep_memory_2)
        self.pep_memory_type_3 = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_memory_type_3.setObjectName("pep_memory_type_3")
        self.formLayout.setWidget(16, QtWidgets.QFormLayout.LabelRole, self.pep_memory_type_3)
        self.pep_memory_3 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory_3.setObjectName("pep_memory_3")
        self.formLayout.setWidget(16, QtWidgets.QFormLayout.FieldRole, self.pep_memory_3)
        self.pep_memory_type_4 = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_memory_type_4.setObjectName("pep_memory_type_4")
        self.formLayout.setWidget(18, QtWidgets.QFormLayout.LabelRole, self.pep_memory_type_4)
        self.pep_memory_4 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory_4.setObjectName("pep_memory_4")
        self.formLayout.setWidget(18, QtWidgets.QFormLayout.FieldRole, self.pep_memory_4)
        self.pep_memory_type_5 = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_memory_type_5.setObjectName("pep_memory_type_5")
        self.formLayout.setWidget(20, QtWidgets.QFormLayout.LabelRole, self.pep_memory_type_5)
        self.pep_memory_5 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory_5.setObjectName("pep_memory_5")
        self.formLayout.setWidget(20, QtWidgets.QFormLayout.FieldRole, self.pep_memory_5)
        self.pep_memory_type_6 = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_memory_type_6.setObjectName("pep_memory_type_6")
        self.formLayout.setWidget(22, QtWidgets.QFormLayout.LabelRole, self.pep_memory_type_6)
        self.pep_transcriptome = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.pep_transcriptome.setObjectName("pep_transcriptome")
        self.formLayout.setWidget(25, QtWidgets.QFormLayout.LabelRole, self.pep_transcriptome)
        self.pep_memory_6 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory_6.setObjectName("pep_memory_6")
        self.formLayout.setWidget(22, QtWidgets.QFormLayout.FieldRole, self.pep_memory_6)
        self.spinBox_3 = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.spinBox_3.setObjectName("spinBox_3")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.spinBox_3)
        self.pep_instrument_3 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_instrument_3.setObjectName("pep_instrument_3")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.pep_instrument_3)
        self.spinBox_4 = QtWidgets.QSpinBox(self.scrollAreaWidgetContents)
        self.spinBox_4.setObjectName("spinBox_4")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.spinBox_4)
        self.pep_instrument_4 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_instrument_4.setObjectName("pep_instrument_4")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.pep_instrument_4)
        self.pep_tolerance_type_2 = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_tolerance_type_2.setObjectName("pep_tolerance_type_2")
        self.formLayout.setWidget(23, QtWidgets.QFormLayout.LabelRole, self.pep_tolerance_type_2)
        self.pep_memory_7 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory_7.setObjectName("pep_memory_7")
        self.formLayout.setWidget(23, QtWidgets.QFormLayout.FieldRole, self.pep_memory_7)
        self.pep_tolerance_type_3 = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pep_tolerance_type_3.setObjectName("pep_tolerance_type_3")
        self.formLayout.setWidget(24, QtWidgets.QFormLayout.LabelRole, self.pep_tolerance_type_3)
        self.pep_memory_8 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pep_memory_8.setObjectName("pep_memory_8")
        self.formLayout.setWidget(24, QtWidgets.QFormLayout.FieldRole, self.pep_memory_8)
        self.pep_specifications_area.setWidget(self.scrollAreaWidgetContents)
        self.pep_specifications = QtWidgets.QLabel(self.PepSearch)
        self.pep_specifications.setGeometry(QtCore.QRect(410, 0, 91, 16))
        self.pep_specifications.setObjectName("pep_specifications")
        self.scrollArea = QtWidgets.QScrollArea(self.PepSearch)
        self.scrollArea.setGeometry(QtCore.QRect(10, 150, 381, 351))
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents_2 = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_2.setGeometry(QtCore.QRect(0, 0, 362, 456))
        self.scrollAreaWidgetContents_2.setObjectName("scrollAreaWidgetContents_2")
        self.gridLayout_24 = QtWidgets.QGridLayout(self.scrollAreaWidgetContents_2)
        self.gridLayout_24.setObjectName("gridLayout_24")
        self.mod1 = QtWidgets.QGroupBox(self.scrollAreaWidgetContents_2)
        self.mod1.setObjectName("mod1")
        self.gridLayout = QtWidgets.QGridLayout(self.mod1)
        self.gridLayout.setObjectName("gridLayout")
        self.position_type_mod1 = QtWidgets.QLineEdit(self.mod1)
        self.position_type_mod1.setObjectName("position_type_mod1")
        self.gridLayout.addWidget(self.position_type_mod1, 2, 2, 1, 1)
        self.type_type_mod1 = QtWidgets.QLineEdit(self.mod1)
        self.type_type_mod1.setObjectName("type_type_mod1")
        self.gridLayout.addWidget(self.type_type_mod1, 2, 0, 1, 1)
        self.aa_type_mod1 = QtWidgets.QLineEdit(self.mod1)
        self.aa_type_mod1.setObjectName("aa_type_mod1")
        self.gridLayout.addWidget(self.aa_type_mod1, 2, 5, 1, 1)
        self.mass_mod1 = QtWidgets.QLabel(self.mod1)
        self.mass_mod1.setObjectName("mass_mod1")
        self.gridLayout.addWidget(self.mass_mod1, 1, 6, 1, 1)
        self.aa_mod1 = QtWidgets.QLabel(self.mod1)
        self.aa_mod1.setObjectName("aa_mod1")
        self.gridLayout.addWidget(self.aa_mod1, 2, 6, 1, 1)
        self.position_mod1 = QtWidgets.QLabel(self.mod1)
        self.position_mod1.setObjectName("position_mod1")
        self.gridLayout.addWidget(self.position_mod1, 2, 4, 1, 1)
        self.name_mod1 = QtWidgets.QLabel(self.mod1)
        self.name_mod1.setObjectName("name_mod1")
        self.gridLayout.addWidget(self.name_mod1, 0, 6, 1, 1)
        self.type_mod1 = QtWidgets.QLabel(self.mod1)
        self.type_mod1.setObjectName("type_mod1")
        self.gridLayout.addWidget(self.type_mod1, 2, 1, 1, 1)
        self.mass_type_mod1 = QtWidgets.QLineEdit(self.mod1)
        self.mass_type_mod1.setObjectName("mass_type_mod1")
        self.gridLayout.addWidget(self.mass_type_mod1, 1, 0, 1, 5)
        self.name_type_mod1 = QtWidgets.QLineEdit(self.mod1)
        self.name_type_mod1.setObjectName("name_type_mod1")
        self.gridLayout.addWidget(self.name_type_mod1, 0, 0, 1, 5)
        self.gridLayout_24.addWidget(self.mod1, 0, 0, 1, 1)
        self.mod2 = QtWidgets.QGroupBox(self.scrollAreaWidgetContents_2)
        self.mod2.setObjectName("mod2")
        self.gridLayout_17 = QtWidgets.QGridLayout(self.mod2)
        self.gridLayout_17.setObjectName("gridLayout_17")
        self.type_type_mod2 = QtWidgets.QLineEdit(self.mod2)
        self.type_type_mod2.setObjectName("type_type_mod2")
        self.gridLayout_17.addWidget(self.type_type_mod2, 2, 0, 1, 1)
        self.position_type_mod2 = QtWidgets.QLineEdit(self.mod2)
        self.position_type_mod2.setObjectName("position_type_mod2")
        self.gridLayout_17.addWidget(self.position_type_mod2, 2, 2, 1, 1)
        self.name_type_mod2 = QtWidgets.QLineEdit(self.mod2)
        self.name_type_mod2.setObjectName("name_type_mod2")
        self.gridLayout_17.addWidget(self.name_type_mod2, 0, 0, 1, 5)
        self.aa_type_mod2 = QtWidgets.QLineEdit(self.mod2)
        self.aa_type_mod2.setObjectName("aa_type_mod2")
        self.gridLayout_17.addWidget(self.aa_type_mod2, 2, 5, 1, 1)
        self.aa_mod2 = QtWidgets.QLabel(self.mod2)
        self.aa_mod2.setObjectName("aa_mod2")
        self.gridLayout_17.addWidget(self.aa_mod2, 2, 6, 1, 1)
        self.mass_type_mod2 = QtWidgets.QLineEdit(self.mod2)
        self.mass_type_mod2.setObjectName("mass_type_mod2")
        self.gridLayout_17.addWidget(self.mass_type_mod2, 1, 0, 1, 5)
        self.type_mod2 = QtWidgets.QLabel(self.mod2)
        self.type_mod2.setObjectName("type_mod2")
        self.gridLayout_17.addWidget(self.type_mod2, 2, 1, 1, 1)
        self.mass_mod2 = QtWidgets.QLabel(self.mod2)
        self.mass_mod2.setObjectName("mass_mod2")
        self.gridLayout_17.addWidget(self.mass_mod2, 1, 6, 1, 1)
        self.name_mod2 = QtWidgets.QLabel(self.mod2)
        self.name_mod2.setObjectName("name_mod2")
        self.gridLayout_17.addWidget(self.name_mod2, 0, 6, 1, 1)
        self.position_mod2 = QtWidgets.QLabel(self.mod2)
        self.position_mod2.setObjectName("position_mod2")
        self.gridLayout_17.addWidget(self.position_mod2, 2, 4, 1, 1)
        self.gridLayout_24.addWidget(self.mod2, 1, 0, 1, 1)
        self.mod3 = QtWidgets.QGroupBox(self.scrollAreaWidgetContents_2)
        self.mod3.setObjectName("mod3")
        self.gridLayout_18 = QtWidgets.QGridLayout(self.mod3)
        self.gridLayout_18.setObjectName("gridLayout_18")
        self.mass_mod3 = QtWidgets.QLabel(self.mod3)
        self.mass_mod3.setObjectName("mass_mod3")
        self.gridLayout_18.addWidget(self.mass_mod3, 1, 6, 1, 1)
        self.aa_mod3 = QtWidgets.QLabel(self.mod3)
        self.aa_mod3.setObjectName("aa_mod3")
        self.gridLayout_18.addWidget(self.aa_mod3, 2, 6, 1, 1)
        self.type_type_mod3 = QtWidgets.QLineEdit(self.mod3)
        self.type_type_mod3.setObjectName("type_type_mod3")
        self.gridLayout_18.addWidget(self.type_type_mod3, 2, 0, 1, 1)
        self.type_mod3 = QtWidgets.QLabel(self.mod3)
        self.type_mod3.setObjectName("type_mod3")
        self.gridLayout_18.addWidget(self.type_mod3, 2, 1, 1, 1)
        self.name_type_mod3 = QtWidgets.QLineEdit(self.mod3)
        self.name_type_mod3.setObjectName("name_type_mod3")
        self.gridLayout_18.addWidget(self.name_type_mod3, 0, 0, 1, 5)
        self.mass_type_mod3 = QtWidgets.QLineEdit(self.mod3)
        self.mass_type_mod3.setObjectName("mass_type_mod3")
        self.gridLayout_18.addWidget(self.mass_type_mod3, 1, 0, 1, 5)
        self.aa_type_mod3 = QtWidgets.QLineEdit(self.mod3)
        self.aa_type_mod3.setObjectName("aa_type_mod3")
        self.gridLayout_18.addWidget(self.aa_type_mod3, 2, 5, 1, 1)
        self.name_mod3 = QtWidgets.QLabel(self.mod3)
        self.name_mod3.setObjectName("name_mod3")
        self.gridLayout_18.addWidget(self.name_mod3, 0, 6, 1, 1)
        self.position_type_mod3 = QtWidgets.QLineEdit(self.mod3)
        self.position_type_mod3.setObjectName("position_type_mod3")
        self.gridLayout_18.addWidget(self.position_type_mod3, 2, 2, 1, 1)
        self.position_mod3 = QtWidgets.QLabel(self.mod3)
        self.position_mod3.setObjectName("position_mod3")
        self.gridLayout_18.addWidget(self.position_mod3, 2, 4, 1, 1)
        self.gridLayout_24.addWidget(self.mod3, 2, 0, 1, 1)
        self.mod4 = QtWidgets.QGroupBox(self.scrollAreaWidgetContents_2)
        self.mod4.setObjectName("mod4")
        self.gridLayout_23 = QtWidgets.QGridLayout(self.mod4)
        self.gridLayout_23.setObjectName("gridLayout_23")
        self.name_mod4 = QtWidgets.QLabel(self.mod4)
        self.name_mod4.setObjectName("name_mod4")
        self.gridLayout_23.addWidget(self.name_mod4, 0, 6, 1, 1)
        self.position_type_mod4 = QtWidgets.QLineEdit(self.mod4)
        self.position_type_mod4.setObjectName("position_type_mod4")
        self.gridLayout_23.addWidget(self.position_type_mod4, 2, 2, 1, 1)
        self.mass_type_mod4 = QtWidgets.QLineEdit(self.mod4)
        self.mass_type_mod4.setObjectName("mass_type_mod4")
        self.gridLayout_23.addWidget(self.mass_type_mod4, 1, 0, 1, 5)
        self.type_type_mod4 = QtWidgets.QLineEdit(self.mod4)
        self.type_type_mod4.setObjectName("type_type_mod4")
        self.gridLayout_23.addWidget(self.type_type_mod4, 2, 0, 1, 1)
        self.aa_type_mod4 = QtWidgets.QLineEdit(self.mod4)
        self.aa_type_mod4.setObjectName("aa_type_mod4")
        self.gridLayout_23.addWidget(self.aa_type_mod4, 2, 5, 1, 1)
        self.position_mod4 = QtWidgets.QLabel(self.mod4)
        self.position_mod4.setObjectName("position_mod4")
        self.gridLayout_23.addWidget(self.position_mod4, 2, 4, 1, 1)
        self.mass_mod4 = QtWidgets.QLabel(self.mod4)
        self.mass_mod4.setObjectName("mass_mod4")
        self.gridLayout_23.addWidget(self.mass_mod4, 1, 6, 1, 1)
        self.aa_mod4 = QtWidgets.QLabel(self.mod4)
        self.aa_mod4.setObjectName("aa_mod4")
        self.gridLayout_23.addWidget(self.aa_mod4, 2, 6, 1, 1)
        self.type_mod4 = QtWidgets.QLabel(self.mod4)
        self.type_mod4.setObjectName("type_mod4")
        self.gridLayout_23.addWidget(self.type_mod4, 2, 1, 1, 1)
        self.name_type_mod4 = QtWidgets.QLineEdit(self.mod4)
        self.name_type_mod4.setObjectName("name_type_mod4")
        self.gridLayout_23.addWidget(self.name_type_mod4, 0, 0, 1, 5)
        self.gridLayout_24.addWidget(self.mod4, 3, 0, 1, 1)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents_2)
        self.Modifications_label = QtWidgets.QLabel(self.PepSearch)
        self.Modifications_label.setGeometry(QtCore.QRect(10, 130, 111, 16))
        self.Modifications_label.setObjectName("Modifications_label")
        self.groupBox = QtWidgets.QGroupBox(self.PepSearch)
        self.groupBox.setGeometry(QtCore.QRect(10, 14, 381, 111))
        self.groupBox.setObjectName("groupBox")

        ## Output directory ##
        self.pep_oudir_find = QtWidgets.QPushButton(self.groupBox)
        self.pep_oudir_find.setGeometry(QtCore.QRect(340, 70, 31, 21))
        self.pep_oudir_find.setObjectName("pep_oudir_find")
        self.pep_outdir = QtWidgets.QLabel(self.groupBox)
        self.pep_outdir.setGeometry(QtCore.QRect(10, 70, 101, 16))
        self.pep_outdir.setObjectName("pep_outdir")
        self.pep_oudir_find.clicked.connect(lambda: self.find_directory("outdir", "YES"))


        ## Mass Spec directory ##
        self.pep_oudir_find_2 = QtWidgets.QPushButton(self.groupBox)
        self.pep_oudir_find_2.setGeometry(QtCore.QRect(340, 30, 31, 21))
        self.pep_oudir_find_2.setObjectName("pep_oudir_find_2")
        self.pep_outdir_2 = QtWidgets.QLabel(self.groupBox)
        self.pep_outdir_2.setGeometry(QtCore.QRect(10, 30, 121, 16))
        self.pep_outdir_2.setObjectName("pep_outdir_2")
        self.pep_oudir_find_2.clicked.connect(lambda: self.find_directory("Mass_spec", "xablau"))


        ## Peptide Search - Run ProteInS ##
        self.pep_run_proteins = QtWidgets.QToolButton(self.PepSearch)
        self.pep_run_proteins.setGeometry(QtCore.QRect(410, 470, 371, 31))
        self.pep_run_proteins.setObjectName("pep_run_proteins")
        self.pep_run_proteins.clicked.connect(lambda: self.set_pep_search_args())


        self.Modes.addTab(self.PepSearch, "")

        """Post MS mode."""
        self.postms = QtWidgets.QWidget()
        self.postms.setObjectName("postms")
        self.Folders_pms = QtWidgets.QGroupBox(self.postms)
        self.Folders_pms.setGeometry(QtCore.QRect(10, 10, 381, 121))
        self.Folders_pms.setObjectName("Folders_pms")
        self.pms_outdir_find = QtWidgets.QPushButton(self.Folders_pms)
        self.pms_outdir_find.setGeometry(QtCore.QRect(340, 70, 31, 21))
        self.pms_outdir_find.setObjectName("pms_outdir_find")
        self.pms_outdir_find.clicked.connect(lambda: self.find_directory("outdir", "YES"))

        self.pms_outdir = QtWidgets.QLabel(self.Folders_pms)
        self.pms_outdir.setGeometry(QtCore.QRect(10, 70, 81, 16))
        self.pms_outdir.setObjectName("pms_outdir")
        self.pms_mass_spec_find = QtWidgets.QPushButton(self.Folders_pms)
        self.pms_mass_spec_find.setGeometry(QtCore.QRect(340, 30, 31, 21))
        self.pms_mass_spec_find.setObjectName("pms_mass_spec_find")
        self.pms_mass_spec_find.clicked.connect(lambda: self.find_directory("Mass_spec", "xablau"))

        self.pms_mass_spec = QtWidgets.QLabel(self.Folders_pms)
        self.pms_mass_spec.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.pms_mass_spec.setObjectName("pms_mass_spec")
        self.Files_pms = QtWidgets.QGroupBox(self.postms)
        self.Files_pms.setGeometry(QtCore.QRect(10, 140, 381, 151))
        self.Files_pms.setObjectName("Files_pms")
        self.proteome_pms_find = QtWidgets.QPushButton(self.Files_pms)
        self.proteome_pms_find.setGeometry(QtCore.QRect(340, 30, 31, 21))
        self.proteome_pms_find.setObjectName("proteome_pms_find")
        self.proteome_pms_find.clicked.connect(lambda: self.find_file("proteome"))

        self.proteome_pms = QtWidgets.QLabel(self.Files_pms)
        self.proteome_pms.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.proteome_pms.setObjectName("proteome_pms")
        self.genome_pms_find = QtWidgets.QPushButton(self.Files_pms)
        self.genome_pms_find.setGeometry(QtCore.QRect(340, 70, 31, 21))
        self.genome_pms_find.setObjectName("genome_pms_find")
        self.genome_pms_find.clicked.connect(lambda: self.find_file("genome"))

        self.genome_pms = QtWidgets.QLabel(self.Files_pms)
        self.genome_pms.setGeometry(QtCore.QRect(10, 70, 91, 16))
        self.genome_pms.setObjectName("genome_pms")
        self.genbank_pms_find = QtWidgets.QPushButton(self.Files_pms)
        self.genbank_pms_find.setGeometry(QtCore.QRect(340, 110, 31, 21))
        self.genbank_pms_find.setObjectName("genbank_pms_find")
        self.genbank_pms_find.clicked.connect(lambda: self.find_file("genbank"))

        self.genbank_pms = QtWidgets.QLabel(self.Files_pms)
        self.genbank_pms.setGeometry(QtCore.QRect(10, 110, 91, 16))
        self.genbank_pms.setObjectName("genbank_pms")
        self.specifications_pms = QtWidgets.QGroupBox(self.postms)
        self.specifications_pms.setGeometry(QtCore.QRect(410, 10, 361, 71))
        self.specifications_pms.setObjectName("specifications_pms")
        self.pms_transcriptome = QtWidgets.QCheckBox(self.specifications_pms)
        self.pms_transcriptome.setGeometry(QtCore.QRect(20, 30, 91, 17))
        self.pms_transcriptome.setObjectName("pms_transcriptome")
        self.pms_run_proteins = QtWidgets.QToolButton(self.postms)
        self.pms_run_proteins.setGeometry(QtCore.QRect(410, 90, 361, 31))
        self.pms_run_proteins.setObjectName("pms_run_proteins")
        self.pms_run_proteins.clicked.connect(lambda: self.set_postms_args())

        self.Modes.addTab(self.postms, "")

        """ Homologues mode """
        self.homologues = QtWidgets.QWidget()
        self.homologues.setObjectName("homologues")
        self.Homology_box = QtWidgets.QGroupBox(self.homologues)
        self.Homology_box.setGeometry(QtCore.QRect(10, 10, 391, 81))
        self.Homology_box.setObjectName("Homology_box")
        self.output_homo_find = QtWidgets.QPushButton(self.Homology_box)
        self.output_homo_find.setGeometry(QtCore.QRect(340, 30, 31, 21))
        self.output_homo_find.setObjectName("output_homo_find")
        self.output_homo = QtWidgets.QLabel(self.Homology_box)
        self.output_homo.setGeometry(QtCore.QRect(10, 30, 101, 16))
        self.output_homo.setObjectName("output_homo")
        self.specifications_homo = QtWidgets.QGroupBox(self.homologues)
        self.specifications_homo.setGeometry(QtCore.QRect(420, 10, 361, 121))
        self.specifications_homo.setObjectName("specifications_homo")
        self.organisms_homo = QtWidgets.QLabel(self.specifications_homo)
        self.organisms_homo.setGeometry(QtCore.QRect(280, 30, 71, 20))
        self.organisms_homo.setObjectName("organisms_homo")
        self.organisms_homo_type = QtWidgets.QLineEdit(self.specifications_homo)
        self.organisms_homo_type.setGeometry(QtCore.QRect(7, 30, 261, 20))
        self.organisms_homo_type.setObjectName("organisms_homo_type")
        self.table_homo = QtWidgets.QLabel(self.specifications_homo)
        self.table_homo.setGeometry(QtCore.QRect(50, 70, 195, 20))
        self.table_homo.setObjectName("table_homo")
        self.table_homo_type = QtWidgets.QSpinBox(self.specifications_homo)
        self.table_homo_type.setGeometry(QtCore.QRect(11, 70, 27, 20))
        self.table_homo_type.setMinimum(1)
        self.table_homo_type.setMaximum(33)
        self.table_homo_type.setObjectName("table_homo_type")
        self.homo_run_proteins = QtWidgets.QToolButton(self.homologues)
        self.homo_run_proteins.setGeometry(QtCore.QRect(420, 140, 361, 31))
        self.homo_run_proteins.setObjectName("homo_run_proteins")
        self.Modes.addTab(self.homologues, "")

        """ Visualization mode """
        self.visualization = QtWidgets.QWidget()
        self.visualization.setObjectName("visualization")
        self.folders_visu = QtWidgets.QGroupBox(self.visualization)
        self.folders_visu.setGeometry(QtCore.QRect(10, 10, 391, 71))
        self.folders_visu.setObjectName("folders_visu")
        self.output_visu_find = QtWidgets.QPushButton(self.folders_visu)
        self.output_visu_find.setGeometry(QtCore.QRect(340, 30, 31, 21))
        self.output_visu_find.setObjectName("output_visu_find")
        self.output_visu = QtWidgets.QLabel(self.folders_visu)
        self.output_visu.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.output_visu.setObjectName("output_visu")
        self.files_box_visu = QtWidgets.QGroupBox(self.visualization)
        self.files_box_visu.setGeometry(QtCore.QRect(10, 90, 391, 121))
        self.files_box_visu.setObjectName("files_box_visu")
        self.genome_visu_find = QtWidgets.QPushButton(self.files_box_visu)
        self.genome_visu_find.setGeometry(QtCore.QRect(340, 30, 31, 21))
        self.genome_visu_find.setObjectName("genome_visu_find")
        self.genome_visu = QtWidgets.QLabel(self.files_box_visu)
        self.genome_visu.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.genome_visu.setObjectName("genome_visu")
        self.genbank_visu = QtWidgets.QLabel(self.files_box_visu)
        self.genbank_visu.setGeometry(QtCore.QRect(10, 80, 91, 16))
        self.genbank_visu.setObjectName("genbank_visu")
        self.genbank_visu_find = QtWidgets.QPushButton(self.files_box_visu)
        self.genbank_visu_find.setGeometry(QtCore.QRect(340, 80, 31, 21))
        self.genbank_visu_find.setObjectName("genbank_visu_find")
        self.specifications_visu = QtWidgets.QGroupBox(self.visualization)
        self.specifications_visu.setGeometry(QtCore.QRect(420, 10, 361, 91))
        self.specifications_visu.setObjectName("specifications_visu")
        self.dataset_visu = QtWidgets.QLabel(self.specifications_visu)
        self.dataset_visu.setGeometry(QtCore.QRect(60, 40, 195, 20))
        self.dataset_visu.setObjectName("dataset_visu")
        self.dataset_visu_type = QtWidgets.QLineEdit(self.specifications_visu)
        self.dataset_visu_type.setGeometry(QtCore.QRect(10, 40, 41, 20))
        self.dataset_visu_type.setObjectName("dataset_visu_type")
        self.visu_run_proteins = QtWidgets.QToolButton(self.visualization)
        self.visu_run_proteins.setGeometry(QtCore.QRect(420, 110, 361, 31))
        self.visu_run_proteins.setObjectName("visu_run_proteins")
        self.Modes.addTab(self.visualization, "")


        ProteInS.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(ProteInS)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        ProteInS.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(ProteInS)
        self.statusbar.setObjectName("statusbar")
        ProteInS.setStatusBar(self.statusbar)
        self.actionOpen = QtWidgets.QAction(ProteInS)
        self.actionOpen.setObjectName("actionOpen")
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionOpen)
        self.menubar.addAction(self.menuFile.menuAction())


        self.retranslateUi(ProteInS)
        self.Modes.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(ProteInS)

    def find_file(self, file):
        """ Returns the path to the chosen file. """
        filePath = QtWidgets.QFileDialog.getOpenFileName()
        self.files[file] = filePath[0]

    def find_files(self, file):
        """ Returns the paths to the chosen files. """
        filePaths = QtWidgets.QFileDialog.getOpenFileNames()
        files = filePaths
        i = 1
        for_arg = ""
        for f in files[0]:
            if len(files) == i:
                for_arg += "{}".format(f)
            else:
                for_arg += "{},".format(f)
            i += 1
        self.files[file] = for_arg

    def find_directory(self, directory, output):
        if output == "YES":
            try:
                """Returns the path to the chosen directory. """
                dirPath = QtWidgets.QFileDialog.getExistingDirectory()
                print(dirPath)
                self.files[directory] = os.path.abspath("{}/ProteInS".format(dirPath))
                print(self.files[directory])
            except:
                print("Please inform the directory!")
        else:
            try:
                dirPath = QtWidgets.QFileDialog.getExistingDirectory()
                self.files[directory] = "{}".format(dirPath[0])

            except:
                print("Please inform the directory!")

    def get_text(self, gui_self, arg):
        """ Returns text of QLineEdit object - a text box. """
        text = gui_self.text()
        self.files[arg] = text

    def save_text(self, gui_self):
        text = gui_self.text()
        return text

    def test_check_box(self, box, arg):
        """ Checks value of QCheckBox object. """
        if box.isChecked() == True:
            self.files[arg] = "YES"

    def get_spin_value(self, box, arg):
        """ Returns the value of a QSpinBox - a numerical box containing a selection of numbers. """
        number = box.value()
        self.files[arg] = number

    def set_assembly_args(self):
        """ Should be assigned to a "Run ProteInS" button. """
        self.get_text(self.rna_cpu_type, "cpu")
        self.get_text(self.rna_memory_type, "memory")
        args = gargs.Arguments(self.files)
        args.assembly()
        print(sys.path)

    def set_db_args(self):
        """ Should be assigned to a "Run ProteInS" button. """
        self.get_text(self.db_starts_type, "starts")
        self.get_text(self.db_ends_type, "ends")
        self.get_text(self.db_maxsize_type, "maxsize")
        self.get_text(self.db_minsize_type, "minsize")
        self.test_check_box(self.db_transcriptome, "Transcriptome")
        args = gargs.Arguments(self.files)
        args.database()

    def organize_modifications(self):
        """ Should be called inside set_pep_search_args(). Organizes all inputs for peptide modifications in order to
         make it readable for the rest of the pipeline. """
        name1 = self.save_text(self.name_type_mod1)
        mass1 = self.save_text(self.mass_type_mod1)
        type1 = self.save_text(self.type_type_mod1)
        pos1 = self.save_text(self.position_type_mod1)
        aa1 = self.save_text(self.aa_type_mod1)
        if name1 != "":
            mod1 = "{},{},{},{},{}".format(name1, mass1, aa1, type1, pos1)
            self.files["modification1"] = mod1
        name2 = self.save_text(self.name_type_mod2)
        mass2 = self.save_text(self.mass_type_mod2)
        type2 = self.save_text(self.type_type_mod2)
        pos2 = self.save_text(self.position_type_mod2)
        aa2 = self.save_text(self.aa_type_mod2)
        if name2 != "":
            mod2 = "{},{},{},{},{}".format(name2, mass2, aa2, type2, pos2)
            self.files["modification2"] = mod2
        name3 = self.save_text(self.name_type_mod3)
        mass3 = self.save_text(self.mass_type_mod3)
        type3 = self.save_text(self.type_type_mod3)
        pos3 = self.save_text(self.position_type_mod3)
        aa3 = self.save_text(self.aa_type_mod3)
        if name3 != "":
            mod3 = "{},{},{},{},{}".format(name3, mass3, aa3, type3, pos3)
            self.files["modification3"] = mod3
        name4 = self.save_text(self.name_type_mod4)
        mass4 = self.save_text(self.mass_type_mod4)
        type4 = self.save_text(self.type_type_mod4)
        pos4 = self.save_text(self.position_type_mod4)
        aa4 = self.save_text(self.aa_type_mod4)
        if name4 != "":
            mod4 = "{},{},{},{},{}".format(name4, mass4, aa4, type4, pos4)
            self.files["modification4"] = mod4

    def set_pep_search_args(self):
        """ Should be assigned to a "Run ProteInS" button. """
        self.get_spin_value(self.pep_enzyme_type, "enzyme")
        self.get_spin_value(self.pep_modifications_type, "mods")
        self.get_spin_value(self.pep_fragmentation_type, "fragmentation")
        self.get_spin_value(self.pep_protocol_type, "protocol")
        self.get_spin_value(self.pep_instrument_type, "instrument")
        self.get_spin_value(self.spinBox, "termini")
        self.get_spin_value(self.spinBox_3, "missCleavages")
        self.get_spin_value(self.spinBox_4, "irrCleavages")
        self.get_text(self.pep_tolerance_type, "tolerance")
        self.get_text(self.pep_decoy_type, "decoy")
        self.get_text(self.pep_memory_type, "memory")
        self.get_text(self.pep_memory_type_2, "isotope")
        self.get_text(self.pep_memory_type_3, "length")
        self.get_text(self.pep_memory_type_4, "charge")
        self.get_text(self.pep_memory_type_5, "matches")
        self.get_text(self.pep_memory_type_6, "numMods")
        self.get_text(self.pep_tolerance_type_2, "ParentError")
        self.get_text(self.pep_tolerance_type_3, "Fdr")
        self.test_check_box(self.pep_transcriptome, "Transcriptome")
        self.organize_modifications()
        args = gargs.Arguments(self.files)
        args.peptide_search()

    def set_postms_args(self):
        self.test_check_box(self.pms_transcriptome, "Transcriptome")
        args = gargs.Arguments(self.files)
        print(args.dic)
        args.postms()

    def retranslateUi(self, ProteInS):
        _translate = QtCore.QCoreApplication.translate
        ProteInS.setWindowTitle(_translate("ProteInS", "ProteInS"))
        self.Modes.setToolTip(_translate("ProteInS", "<html><head/><body><p><br/></p></body></html>"))
        self.Modes.setWhatsThis(_translate("ProteInS", "<html><head/><body><p>Database</p></body></html>"))
        self.read1.setText(_translate("ProteInS", "Left read(s)"))
        self.read2.setText(_translate("ProteInS", "Right read(s)"))
        self.find_read1.setText(_translate("ProteInS", "..."))
        self.find_read2.setText(_translate("ProteInS", "..."))
        self.specifics.setTitle(_translate("ProteInS", "Specifications"))
        self.rna_cpu.setText(_translate("ProteInS", "CPU"))
        self.k_mer.setText(_translate("ProteInS", "K-mer length"))
        self.rna_memory.setText(_translate("ProteInS", "Max memory"))
        self.reads_box.setTitle(_translate("ProteInS", "Read files"))
        self.paired_reads_box.setTitle(_translate("ProteInS", "Paired-end reads"))
        self.input_left_read.setText(_translate("ProteInS", "."))
        self.single_end_box.setTitle(_translate("ProteInS", "Single-end reads"))
        self.single_reads_find.setText(_translate("ProteInS", "..."))
        self.single_reads.setText(_translate("ProteInS", "Read(s)"))
        self.rna_output_box.setTitle(_translate("ProteInS", "Output"))
        self.rna_outdir_find.setText(_translate("ProteInS", "..."))
        self.rna_output_dir.setText(_translate("ProteInS", "Output directory"))
        self.rna_run_proteins.setText(_translate("ProteInS", "Run ProteInS"))
        self.Modes.setTabText(self.Modes.indexOf(self.Assembly), _translate("ProteInS", "Transcriptome Assembly"))
        self.db_files.setTitle(_translate("ProteInS", "Files"))
        self.db_genome_find.setText(_translate("ProteInS", "..."))
        self.db_genome.setText(_translate("ProteInS", "Genome file"))
        self.db_proteome_find.setText(_translate("ProteInS", "..."))
        self.db_proteome.setText(_translate("ProteInS", "Proteome file"))
        self.db_specifics.setTitle(_translate("ProteInS", "Specifications"))
        self.db_starts.setText(_translate("ProteInS", "Start codons"))
        self.db_ends.setText(_translate("ProteInS", "Stop codons"))
        self.db_maxsize.setText(_translate("ProteInS", "Max ORF size"))
        self.db_minsize.setText(_translate("ProteInS", "Min ORF size"))
        self.db_transcriptome.setText(_translate("ProteInS", "Transcriptome"))
        self.db_run_proteins.setText(_translate("ProteInS", "Run ProteInS"))
        self.db_output_box.setTitle(_translate("ProteInS", "Output"))
        self.db_outdir_find.setText(_translate("ProteInS", "..."))
        self.db_outdir.setText(_translate("ProteInS", "Output directory"))
        self.Modes.setTabText(self.Modes.indexOf(self.Database), _translate("ProteInS", "Database Generator"))
        self.pep_enzyme_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>0: Unespecific Cleavage</p><p>1: Trypsin</p><p>2: Chymotrypsin</p><p>3: Lys-C</p><p>4: Lys-N</p><p>5: Glutamyl Endopeptidase</p><p>6: Arg-C</p><p>7: Asp-N</p><p>8: alphaLP</p><p>9: no cleavage</p><p><br/></p></body></html>"))
        self.pep_enzyme.setToolTip(_translate("ProteInS", "<html><head/><body><p>0: Unespecific Cleavage</p><p>1: Trypsin</p><p>2: Chymotrypsin</p><p>3: Lys-C</p><p>4: Lys-N</p><p>5: Glutamyl Endopeptidase</p><p>6: Arg-C</p><p>7: Asp-N</p><p>8: alphaLP</p><p>9: no cleavage</p><p><br/></p></body></html>"))
        self.pep_enzyme.setText(_translate("ProteInS", "Enzyme"))
        self.pep_modifications_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>Maximum number of modifications to be searched for. Max: 4</p></body></html>"))
        self.pep_modifications.setToolTip(_translate("ProteInS", "<html><head/><body><p>Maximum number of modifications to be searched for. Max: 4</p><p><br/></p></body></html>"))
        self.pep_modifications.setText(_translate("ProteInS", "Number of Modifications"))
        self.pep_fragmentation_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>1: CID</p><p>2: ETD</p><p>3: HCD</p><p>4: UVPD</p></body></html>"))
        self.pep_fragmentation.setToolTip(_translate("ProteInS", "<html><head/><body><p>1: CID</p><p>2: ETD</p><p>3: HCD</p><p>4: UVPD</p></body></html>"))
        self.pep_fragmentation.setText(_translate("ProteInS", "Fragmentation method"))
        self.pep_protocol_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>0: Automatic</p><p>1: Phosphorylation</p><p>2: iTRAQ</p><p>3: iTRAQPhospho</p><p>4: TMT</p><p>5: Standard</p></body></html>"))
        self.pep_protocol.setToolTip(_translate("ProteInS", "<html><head/><body><p>0: Automatic</p><p>1: Phosphorylation</p><p>2: iTRAQ</p><p>3: iTRAQPhospho</p><p>4: TMT</p><p>5: Standard</p></body></html>"))
        self.pep_protocol.setText(_translate("ProteInS", "Protocol"))
        self.pep_instrument_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>0: Low-res LCQ/LTQ</p><p>1: Orbitrap/FTICR/Lumos</p><p>2: TOF</p><p>3: Q-Exactive</p><p><br/></p></body></html>"))
        self.pep_instrument.setToolTip(_translate("ProteInS", "<html><head/><body><p>0: Low-res LCQ/LTQ</p><p>1: Orbitrap/FTICR/Lumos</p><p>2: TOF</p><p>3: Q-Exactive</p><p><br/></p></body></html>"))
        self.pep_instrument.setText(_translate("ProteInS", "Instrument"))
        self.spinBox.setToolTip(_translate("ProteInS", "<html><head/><body><p>Example for Trypsin:</p><p>0: non-tryptic</p><p>1: semi-tryptic</p><p>2: fully tryptic peptides only</p></body></html>"))
        self.pep_instrument_2.setToolTip(_translate("ProteInS", "<html><head/><body><p>Example for Trypsin:</p><p>0: non-tryptic</p><p>1: semi-tryptic</p><p>2: fully tryptic peptides only</p></body></html>"))
        self.pep_instrument_2.setText(_translate("ProteInS", "Number of Termini"))
        self.pep_tolerance_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>Precursor mass tolerance in ppm. Inform it such as: 50 ppm</p></body></html>"))
        self.pep_tolerance.setToolTip(_translate("ProteInS", "<html><head/><body><p>Precursor mass tolerance in ppm. Inform it such as: 50 ppm</p></body></html>"))
        self.pep_tolerance.setText(_translate("ProteInS", "Tolerance"))
        self.pep_decoy_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>Inform whether you want to search the decoy database or not. TRUE or FALSE.</p></body></html>"))
        self.pep_decoy.setToolTip(_translate("ProteInS", "<html><head/><body><p>Inform whether you want to search the decoy database or not. TRUE or FALSE.</p></body></html>"))
        self.pep_decoy.setText(_translate("ProteInS", "Decoy"))
        self.pep_memory.setText(_translate("ProteInS", "Memory"))
        self.pep_memory_type_2.setToolTip(_translate("ProteInS", "<html><head/><body><p>Inform it such as 0,1</p></body></html>"))
        self.pep_memory_2.setToolTip(_translate("ProteInS", "<html><head/><body><p>Inform it such as 0,1</p></body></html>"))
        self.pep_memory_2.setText(_translate("ProteInS", "Isotope error range"))
        self.pep_memory_type_3.setToolTip(_translate("ProteInS", "<html><head/><body><p>Peptide length range. I.e, for peptides between 6 and 40 amino acids long, type: 6,40</p></body></html>"))
        self.pep_memory_3.setToolTip(_translate("ProteInS", "<html><head/><body><p>Peptide length range. I.e, for peptides between 6 and 40 amino acids long, type: 6,40</p></body></html>"))
        self.pep_memory_3.setText(_translate("ProteInS", "Peptide length range"))
        self.pep_memory_type_4.setToolTip(_translate("ProteInS", "<html><head/><body><p>Minimum and maximum precursor ion charges. E.g. 2,3</p></body></html>"))
        self.pep_memory_4.setToolTip(_translate("ProteInS", "<html><head/><body><p>Minimum and maximum precursor ion charges. E.g. 2,3</p></body></html>"))
        self.pep_memory_4.setText(_translate("ProteInS", "Charges"))
        self.pep_memory_type_5.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of matches per spectrum to be reported.</p></body></html>"))
        self.pep_memory_5.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of matches per spectrum to be reported.</p></body></html>"))
        self.pep_memory_5.setText(_translate("ProteInS", "Matches"))
        self.pep_memory_type_6.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of dynamic/variable modifications per peptide to be reported.</p></body></html>"))
        self.pep_transcriptome.setText(_translate("ProteInS", "Transcriptome"))
        self.pep_memory_6.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of dynamic/variable modifications per peptide to be reported.</p></body></html>"))
        self.pep_memory_6.setText(_translate("ProteInS", "Number of mods per peptide"))
        self.spinBox_3.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of allowed miss cleavages.</p></body></html>"))
        self.pep_instrument_3.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of allowed miss cleavages.</p></body></html>"))
        self.pep_instrument_3.setText(_translate("ProteInS", "Miss cleavages"))
        self.spinBox_4.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of allowed irregular cleavages.</p></body></html>"))
        self.pep_instrument_4.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of allowed irregular cleavages.</p></body></html>"))
        self.pep_instrument_4.setText(_translate("ProteInS", "Irregular cleavages"))
        self.pep_tolerance_type_2.setToolTip(_translate("ProteInS", "<html><head/><body><p>PPM tolerance for the parent ion mass error. Inform it such as: 30</p></body></html>"))
        self.pep_memory_7.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of dynamic/variable modifications per peptide to be reported.</p></body></html>"))
        self.pep_memory_7.setText(_translate("ProteInS", "Parent Ion Mass error"))
        self.pep_tolerance_type_3.setToolTip(_translate("ProteInS", "<html><head/><body><p>False discovery rate. Default is 0,01.</p></body></html>"))
        self.pep_memory_8.setToolTip(_translate("ProteInS", "<html><head/><body><p>Number of dynamic/variable modifications per peptide to be reported.</p></body></html>"))
        self.pep_memory_8.setText(_translate("ProteInS", "False Discovery Rate (FDR)"))
        self.pep_specifications.setText(_translate("ProteInS", "Specifications"))
        self.mod1.setTitle(_translate("ProteInS", "Modification 1"))
        self.mass_mod1.setText(_translate("ProteInS", "Mass"))
        self.aa_mod1.setText(_translate("ProteInS", "Amino acid"))
        self.position_mod1.setText(_translate("ProteInS", "Position"))
        self.name_mod1.setText(_translate("ProteInS", "Name"))
        self.type_mod1.setText(_translate("ProteInS", "Type"))
        self.mod2.setTitle(_translate("ProteInS", "Modification 2"))
        self.aa_mod2.setText(_translate("ProteInS", "Amino acid"))
        self.type_mod2.setText(_translate("ProteInS", "Type"))
        self.mass_mod2.setText(_translate("ProteInS", "Mass"))
        self.name_mod2.setText(_translate("ProteInS", "Name"))
        self.position_mod2.setText(_translate("ProteInS", "Position"))
        self.mod3.setTitle(_translate("ProteInS", "Modification 3"))
        self.mass_mod3.setText(_translate("ProteInS", "Mass"))
        self.aa_mod3.setText(_translate("ProteInS", "Amino acid"))
        self.type_mod3.setText(_translate("ProteInS", "Type"))
        self.name_mod3.setText(_translate("ProteInS", "Name"))
        self.position_mod3.setText(_translate("ProteInS", "Position"))
        self.mod4.setTitle(_translate("ProteInS", "Modification 4"))
        self.name_mod4.setText(_translate("ProteInS", "Name"))
        self.position_mod4.setText(_translate("ProteInS", "Position"))
        self.mass_mod4.setText(_translate("ProteInS", "Mass"))
        self.aa_mod4.setText(_translate("ProteInS", "Amino acid"))
        self.type_mod4.setText(_translate("ProteInS", "Type"))
        self.Modifications_label.setText(_translate("ProteInS", "Modifications"))
        self.groupBox.setTitle(_translate("ProteInS", "Folders"))
        self.pep_oudir_find.setText(_translate("ProteInS", "..."))
        self.pep_outdir.setText(_translate("ProteInS", "Output directory"))
        self.pep_oudir_find_2.setText(_translate("ProteInS", "..."))
        self.pep_outdir_2.setText(_translate("ProteInS", "Mass Spec directory"))
        self.pep_run_proteins.setText(_translate("ProteInS", "Run ProteInS"))
        self.Modes.setTabText(self.Modes.indexOf(self.PepSearch), _translate("ProteInS", "Peptide Search"))
        self.Folders_pms.setTitle(_translate("ProteInS", "Folders"))
        self.pms_outdir_find.setText(_translate("ProteInS", "..."))
        self.pms_outdir.setText(_translate("ProteInS", "Output directory"))
        self.pms_mass_spec_find.setText(_translate("ProteInS", "..."))
        self.pms_mass_spec.setText(_translate("ProteInS", "Mass Spec Folder"))
        self.Files_pms.setTitle(_translate("ProteInS", "Files"))
        self.proteome_pms_find.setText(_translate("ProteInS", "..."))
        self.proteome_pms.setText(_translate("ProteInS", "Proteome file"))
        self.genome_pms_find.setText(_translate("ProteInS", "..."))
        self.genome_pms.setText(_translate("ProteInS", "Genome file"))
        self.genbank_pms_find.setText(_translate("ProteInS", "..."))
        self.genbank_pms.setText(_translate("ProteInS", "Genbank file"))
        self.specifications_pms.setTitle(_translate("ProteInS", "Specifications"))
        self.pms_transcriptome.setText(_translate("ProteInS", "Transcriptome"))
        self.pms_run_proteins.setText(_translate("ProteInS", "Run ProteInS"))
        self.Modes.setTabText(self.Modes.indexOf(self.postms), _translate("ProteInS", "Data Processing"))
        self.Homology_box.setTitle(_translate("ProteInS", "Folders"))
        self.output_homo_find.setText(_translate("ProteInS", "..."))
        self.output_homo.setText(_translate("ProteInS", "Output directory"))
        self.specifications_homo.setTitle(_translate("ProteInS", "Specifications"))
        self.organisms_homo.setText(_translate("ProteInS", "Organisms"))
        self.table_homo.setToolTip(_translate("ProteInS", "<html><head/><body><p><span style=\" font-family:\'Times New Roman\'; font-size:10pt; color:#000000;\">Inform the NCBI genetic code corresponding to the taxa of interest. This will determine the Start and Stop codons used during the search for orthologues.</span></p></body></html>"))
        self.table_homo.setText(_translate("ProteInS", "Genetic code"))

        self.table_homo_type.setToolTip(_translate("ProteInS", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:1; text-indent:0px;\"><span style=\" font-family:\'Times New Roman\'; font-size:10pt; color:#000000;\">Inform the NCBI genetic code corresponding to the taxa of interest. This will determine the Start and Stop codons used during the search for orthologues.</span></p>\n"
"<ul style=\"margin-top: 0px; margin-bottom: 0px; margin-left: 0px; margin-right: 0px; -qt-list-indent: 1;\"><li style=\" font-family:\'Times New Roman\'; font-size:medium; text-decoration: underline; color:#0000ff;\" align=\"justify\" style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1\"><span style=\" font-size:medium;\">1. The Standard Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">2. The Vertebrate Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG3\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">3. The Yeast Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG5\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">5. The Invertebrate Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG6\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">6. The Ciliate, Dasycladacean and Hexamita Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG9\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">9. The Echinoderm and Flatworm Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG10\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">10. The Euplotid Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">11. The Bacterial, Archaeal and Plant Plastid Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG12\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">12. The Alternative Yeast Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG13\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">13. The Ascidian Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG14\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">14. The Alternative Flatworm Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG16\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">16. Chlorophycean Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG21\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">21. Trematode Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG22\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">22. Scenedesmus obliquus Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG23\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">23. Thraustochytrium Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG24\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">24. Rhabdopleuridae Mitochondrial Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG25\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">25. Candidate Division SR1 and Gracilibacteria Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG26\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">26. Pachysolen tannophilus Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG27\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">27. Karyorelict Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG28\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">28. Condylostoma Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG29\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">29. Mesodinium Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG30\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">30. Peritrich Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG31\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">31. Blastocrithidia Nuclear Code</span></a></li>\n"
"<li style=\" font-family:\'Times New Roman\'; font-size:medium; color:#000000;\" align=\"justify\" style=\" margin-top:0px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><a href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG33\"><span style=\" font-size:medium; text-decoration: underline; color:#0000ff;\">33. Cephalodiscidae Mitochondrial UAA-Tyr Code</span></a></li></ul></body></html>"))
        self.homo_run_proteins.setText(_translate("ProteInS", "Run ProteInS"))
        self.Modes.setTabText(self.Modes.indexOf(self.homologues), _translate("ProteInS", "Homology"))
        self.folders_visu.setTitle(_translate("ProteInS", "Folders"))
        self.output_visu_find.setText(_translate("ProteInS", "..."))
        self.output_visu.setText(_translate("ProteInS", "Output directory"))
        self.files_box_visu.setTitle(_translate("ProteInS", "Files"))
        self.genome_visu_find.setText(_translate("ProteInS", "..."))
        self.genome_visu.setText(_translate("ProteInS", "Genome file"))
        self.genbank_visu.setText(_translate("ProteInS", "Genbank file"))
        self.genbank_visu_find.setText(_translate("ProteInS", "..."))
        self.specifications_visu.setTitle(_translate("ProteInS", "Specifications"))
        self.dataset_visu.setToolTip(_translate("ProteInS", "<html><head/><body><p>Choose which small ORF dataset you want to visualize. Type only one of the following, without the quotes: &quot;dna&quot;, &quot;rna&quot; or &quot;both&quot;.</p></body></html>"))
        self.dataset_visu.setText(_translate("ProteInS", "ORF Dataset"))
        self.dataset_visu_type.setToolTip(_translate("ProteInS", "<html><head/><body><p>Choose which small ORF dataset you want to visualize. Type only one of the following, without the quotes: &quot;dna&quot;, &quot;rna&quot; or &quot;both&quot;.</p></body></html>"))
        self.visu_run_proteins.setText(_translate("ProteInS", "Run ProteInS"))
        self.Modes.setTabText(self.Modes.indexOf(self.visualization), _translate("ProteInS", "Visualization"))
        self.menuFile.setTitle(_translate("ProteInS", "File"))
        self.actionOpen.setText(_translate("ProteInS", "Open"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ProteInS = QtWidgets.QMainWindow()
    ui = Ui_ProteInS()
    ui.setupUi(ProteInS)
    ProteInS.show()
    sys.exit(app.exec_())


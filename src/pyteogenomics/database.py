import sqlite3
import os

from Bio import SeqIO
from .orflib import ORFCollection
from .conversion import *
from .frame_translation import *


class DatabaseGenerator(object):
    """ Generates a database to manage all ORFs predicted using Genome or TranscriptomeTranslator classes. It can
    create either a SQL .db or a .xlsl to store the ORFs. the 'name' argument will determine the name of the file to be
    created. 'db_type' is either 'sql' or a 'df'. """
    def __init__(self, name="database", db_type='sql'):
        self.db = None
        self.__check_db_type(db_type, name)

    def __check_db_type(self, db_type, name):
        if db_type == 'sql':
            self.db = SQLDatabase(name)
        elif db_type == 'fasta':
            pass
        else:
            print("Inform a valid database type. It must be either 'sql',  'fasta' or 'df' (in case of a pandas "
                  "Data Frame. ")

    def add_orfs(self, orfs):
        """ Adds ORFs to the database. Provide the ORF object for this function. """
        self.db.add_orfs(orfs)
        return self

    def retrieve(self):
        data = self.db.retrieve()
        return data

    def to_fasta(self, filename="db_orfs.fasta"):
        """ Writes and entries and sequences inside the database to a fasta file. """
        data = self.retrieve()
        to_write = [f">{orf[1]}_{orf[4]}-{orf[5]}_{orf[6]}\n{orf[2]}\n" for orf in data]
        with open(filename, 'w') as fa:
            fa.writelines(to_write)


class FastaDatabase(object):
    def __init__(self):
        pass


class SQLDatabase(object):
    def __init__(self, name):
        self.name = f'{name}.db'
        self.__create_database()

    def __connect(self):
        conn = sqlite3.connect(self.name)
        return conn

    def __create_database(self):
        conn = self.__connect()
        try:
            conn.execute('''CREATE TABLE ORFOME(ID INT PRIMARY KEY     NOT NULL,
                                                NAME    TEXT    NOT NULL,
                                                SEQ     TEXT    NOT NULL,
                                                LENGTH  INT     NOT NULL,
                                                START   INT     NOT NULL,
                                                END     INT     NOT NULL,
                                                STRAND  TEXT    NOT NULL);''')
        except:
            pass

        return self

    def add_orfs(self, orfs):
        """ This method is only to be called within DatabaseGenerator class. """
        try:
            conn = self.__connect()
            cursor = conn.cursor()
            i = 1
            for orf in orfs:
                insert_query = f"""INSERT INTO ORFOME(ID, NAME, SEQ, LENGTH, START, END, STRAND)
                VALUES({i}, '{orf.name}', '{orf.seq}', {orf.__len__()}, {orf.start}, {orf.end}, '{orf.strand}')"""
                i += 1
                cursor.execute(insert_query)
            conn.commit()
            cursor.close()
        except sqlite3.Error as error:
            print("Failed to insert entries into db. Error: ", error)
        finally:
            if conn:
                conn.close()

    def retrieve(self):
        conn = self.__connect()
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM ORFOME")
        rows = cursor.fetchall()
        return rows

# seq = "ATGAGATGCGGC"
# tr = translate(seq)
# print(tr)

# genome_database = DatabaseGenerator(name="genome", db_type='sql')

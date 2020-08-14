import sqlite3
import os
from sqlite3 import Error

"""
    Old file for saving data into sqlit3 database
"""


class DB:
    def __init__(self, path):  # data_write
        self.path = path

    def database_file_exists(self):
        path_to_file = os.path.join(os.getcwd(), self.path, 'generated_params.db')
        if not os.path.exists(path_to_file):
            os.mknod(path_to_file)
        return path_to_file

    def get_path_to_file(self):
        return self.database_file_exists()

    def create_connection(self):
        """ create a database connection to SQLite database """
        conn = None
        try:
            conn = sqlite3.connect(self.get_path_to_file())
            return conn
        except Error as e:
            print(e)

    def create_tables(self):
        table_protracted_speciation_process = """ CREATE TABLE IF NOT EXISTS table_protracted_speciation_process (
                                            incipient_species_extinction_rate REAL PRIMARY KEY,
                                            speciation_initiation_from_orthospecies_rate REAL,
                                            speciation_initiation_from_incipient_species_rate REAL ,
                                            speciation_completion_rate REAL,
                                            orthospecies_extinction_rate REAL,
                                            aincipient_species_extinction_rate REAL
                                        ); """

        table_generate_sample = """ CREATE TABLE IF NOT EXISTS table_generate_sample (
                                            max_time REAL PRIMARY KEY,
                                            num_extant_orthospecies REAL,
                                            num_extant_lineages REAL,
                                            is_retry_on_total_extinction REAL,
                                            max_retries REAL
                                        ); """

        table_seq_gen = """ CREATE TABLE IF NOT EXISTS table_seq_gen (
                                            state_freqs BLOB PRIMARY KEY,
                                            general_rates BLOB
                                        ); """
        conn = self.create_connection()
        try:
            c = conn.cursor()
            c.execute(table_protracted_speciation_process)
            c.execute(table_generate_sample)
            c.execute(table_seq_gen)
        except Error as e:
            print(e)
        finally:
            if conn:
                conn.close()

    def write_to_db(self):
        conn = self.create_connection()
        try:
            c = conn.cursor()
            # testvalues
            c.execute('insert into table_protracted_speciation_process values (?,?,?,?,?,?)', [1, 2, 3, 4, 5, 6])
            c.execute('insert into table_generate_sample values (?,?,?,?,?)', [1, 2, 3, 4, 5])
            c.execute('insert into table_seq_gen values (?,?)', [1, 2])
        except Error as e:
            print(e)
        finally:
            if conn:
                conn.close()

    def read_db(self):
        conn = self.create_connection()
        cur = conn.cursor()
        try:
            with conn:
                cur.execute("SELECT * FROM table_protracted_speciation_process")
                print(cur.fetchall())
        except Error as e:
            print(e)

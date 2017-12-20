"Some basic methods for PQP parsing"

import sqlite3
import matplotlib.pyplot as plt


class Database:
  """
  Basic utility functions to deal with SQLite DBs
  """

  def __init__(self, db_file):
    "Class initialization"
    self.db_file = db_file
    self.conn, self.cur = self.create_connection()

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    self.conn.close()

  def create_connection(self):
    """
    Creates a database connection to the SQLite database
    specified by the db_file

    :param db_file: database file

    :return: Connection object or None
    """
    try:
      conn = sqlite3.connect(self.db_file)
      cur = conn.cursor()
      return conn, cur
    except sqlite3.Error as e:
      print(e)
    return None, None

  def execute_query(self, query, params=None):
    """
    Query all rows in the tasks table

    :param conn: the Connection object

    :return:
    """
    if params:
      self.cur.execute(query, params)
    else:
      self.cur.execute(query)
    rows = self.cur.fetchall()
    return rows

  def patch_decoys(self):
    """
    Patches decoy status of peptides and proteins based on decoy
    status of precursor.
    """
    # Update peptides
    query = """
    UPDATE peptide
    SET decoy = 1
    WHERE id IN
    (
    SELECT pep.id FROM precursor AS pre
    JOIN precursor_peptide_mapping AS prepep ON pre.id = prepep.precursor_id
    JOIN peptide AS pep ON prepep.peptide_id = pep.id
    WHERE pre.decoy = 1
    )"""
    self.execute_query(query)
    # Update proteins
    query = """
    UPDATE protein
    SET decoy = 1
    WHERE id IN
    (
    SELECT pro.id FROM peptide AS pep
    JOIN peptide_protein_mapping AS peppro on pep.id = peppro.peptide_id
    JOIN protein AS pro ON peppro.protein_id = pro.id
    WHERE pep.decoy = 1
    );
    """
    self.execute_query(query)
    # Commit changes
    self.conn.commit()


class Protein:
  """
  Statistics on proteins in the PQP database
  """

  def __init__(self, db):
    "Initialization"
    self.db = db
    self.proteins = self.query_proteins()

  def query_proteins(self, decoy=0):
    """Get a list of all proteins in PQP.

    Protein_accession consisting of multiple proteins are split into
    individual proteins and stored in a list of sets

    :param decoy: 0 targets only; 1 decoys only
    """
    query = """
    SELECT protein_accession FROM protein 
    WHERE decoy = ?"""
    params = (decoy, )
    rows = self.db.execute_query(query, params)
    return self.split_protein_names(rows)

  @staticmethod
  def split_protein_names(proteins):
    """
    Splits protein_accession into individual protein_accessions.
    E.g.: 2/sp|P55011|S12A2_HUMAN/sp|P37108|SRP14_HUMAN ->
    [{sp|P55011|S12A2_HUMAN, sp|P37108|SRP14_HUMAN}]

    :param rows: list of tuples of protein_accessions

    :return: List of individual protein_accession sets
    """
    return [{protein_accession[0].split('/')[1]} for protein_accession in proteins]

  def unique_proteins(self):
    """
    Creates a list of unique protein

    :param protein_set_list: a list of sets of protein_accessions (as
    output by split_protein_names

    :return: a set of unique proteins from a list of sets (output of
    split_protein_names)
    """
    return set([tuple(x) for x in self.proteins])

  def plot_proteins_per_group(self):
    """
    Plots distribution of protein counts in protein groups

    :return: matplotlib plot (display with plt.show())
    """
    pro_distribution = self.proteins_per_group()
    plt.hist(pro_distribution, bins=max(max(pro_distribution), 10))
    plt.title("Distribution of proteins per protein group")
    plt.xlabel("Number of proteins")
    plt.ylabel("Number of protein groups")
    return plt

  def proteins_per_group(self):
    """
    Calculates the number of proteins in each protein group.

    :return: a list with the number of proteins in each group
    """
    return [len(x) for x in self.proteins]

  def plot_peptides_per_protein(self):
    """
    Plots distribution of peptide counts per protein.

    :return: matplotlib plot (display with plt.show())
    """
    pep_distribution = [x[1] for x in self.peptides_per_protein()]
    plt.hist(pep_distribution, bins=max(pep_distribution))
    plt.title("Distribution of peptides per protein")
    plt.xlabel("Number of peptides")
    plt.ylabel("Number of proteins")
    return plt

  def peptides_per_protein(self, decoy=0):
    """
    Count distinct number of peptides per protein

    :param decoy: 0 targets only; 1 decoys only
    :return: list with number of distinct peptides per protein
    """
    query = """
    SELECT protein_accession, COUNT(DISTINCT(pep.modified_sequence)) AS num_peptides
    FROM protein AS pro
    JOIN peptide_protein_mapping AS ppm ON pro.ID = ppm.protein_ID
    JOIN peptide AS pep ON ppm.peptide_id = pep.ID
    WHERE pro.decoy = ?
    GROUP BY protein_accession;
    """
    params = (decoy, )
    return self.db.execute_query(query, params)

  def count_target_proteins(self, decoy=0):
    """
    Counts the number of protein groups

    :param decoy: 0 targets only; 1 decoys only
    :return: the number of protein groups in the selected category
    """
    if decoy:
      query = """
      SELECT COUNT(id) FROM protein 
      WHERE decoy = ?
      """
      params = (decoy, )
    return self.db.execute_query(query, params)


class Peptide:
  """
  Statistics on peptides in the PQP database
  """

  def __init__(self, db):
    "Initialization"
    self.db = db

  def unique_peptides(self, decoy=0):
    """
    Gets a list of all unique non-decoy peptides in the PQP (modifications are preserved)

    :param decoy: 0 targets only; 1 decoys only
    :return: List of unique peptides in PQP
    """
    query = """
    SELECT DISTINCT(modified_sequence) FROM peptide
    WHERE decoy = ?
    """
    params = (decoy, )
    return self.db.execute_query(query, params)

  def get_proteotypic_peptides(self):
    """
    Gets a list of proteotypic peptides

    :return: List of proteotypic peptides
    """
    query = """
    SELECT * FROM
    (
    SELECT COUNT(DISTINCT(protein_accession)) AS num_pro_with_pep FROM peptide AS pep
    JOIN peptide_protein_mapping AS peppro on pep.id = peppro.peptide_id
    JOIN protein AS pro ON peppro.protein_id = pro.id
    WHERE pep.decoy = 0
    GROUP BY pep.id
    ) AS pep_pro_tbl
    WHERE pep_pro_tbl.num_pro_with_pep = 1
    """
    return self.db.execute_query(query)

  def plot_peptide_promiscuity(self):
    """
    Plots distribution of how number of proteins associated with peptides in the PQP library

    :return: matplotlib plot (display with plt.show())
    """
    pep_promiscuity = [x[0] for x in self.peptide_promiscuity()]
    plt.hist(pep_promiscuity, log=True)
    plt.title("Distribution of peptides per protein")
    plt.xlabel("Number of proteins")
    plt.ylabel("log(Number of peptides)")
    return plt

  def peptide_promiscuity(self):
    """
    Gets the distribution of number of proteins associated with
    peptides in the PQP library

    :return: list of the number of protein groups a peptide associates with
    """
    query = """
    SELECT COUNT(DISTINCT(protein_accession)) AS num_pro_with_pep FROM peptide AS pep
    JOIN peptide_protein_mapping AS peppro on pep.id = peppro.peptide_id
    JOIN protein AS pro ON peppro.protein_id = pro.id
    WHERE pep.decoy = 0
    GROUP BY pep.id
    """
    return self.db.execute_query(query)

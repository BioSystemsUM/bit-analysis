import gzip
import json
import os
import time
import requests
import shutil
import urllib.request as urllib_request
from contextlib import closing
from typing import List, Union

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

from settings import ENTREZ_E_MAIL, ENTREZ_API_KEY


def _get(url) -> requests.Response:
    try:

        response = requests.get(url)

        try:
            response.raise_for_status()

        except requests.exceptions.HTTPError as exception:
            response = None

    except requests.exceptions.RequestException as exception:
        response = None

    if response is None:
        response = requests.Response()
        response.status_code = 404

    return response


class Organism:

    def __init__(self,
                 bigg_id=None,
                 name=None):

        self._bigg_id = bigg_id
        self._name = name
        self._taxonomy_record = {}
        self._genome_record = {}
        self._assembly_record = {}

        self._status = False

    # ------------------------------
    # Polymorphism
    # ------------------------------
    @classmethod
    def from_dict(cls, model_dict):
        bigg_id = model_dict.get('bigg_id')
        organism = model_dict.get('organism')

        return cls(bigg_id=bigg_id,
                   name=organism)

    # ------------------------------
    # Built-in
    # ------------------------------
    def __str__(self):
        return f'{self.bigg_id}_{self.name}'

    def __repr__(self):
        return f'Organism: {self.bigg_id}'

    # ------------------------------
    # Static attributes
    # ------------------------------
    @property
    def bigg_id(self):
        return self._bigg_id

    @property
    def name(self):
        return self._name

    @property
    def status(self):
        return self._status

    @property
    def taxonomy_record(self):
        return self._taxonomy_record

    @taxonomy_record.setter
    def taxonomy_record(self, value):

        if value:
            value = dict(value)

        self._taxonomy_record = value

    @property
    def genome_record(self):
        return self._genome_record

    @genome_record.setter
    def genome_record(self, value):

        if value:
            value = dict(value)

        self._genome_record = value

    @property
    def assembly_record(self):
        return self._assembly_record

    @assembly_record.setter
    def assembly_record(self, value):

        if value:
            value = dict(value)

        self._assembly_record = value

    # ------------------------------
    # Dynamic attributes
    # ------------------------------
    @property
    def taxonomy_ids(self) -> list:
        return [str(identifier) for identifier in self._taxonomy_record.get('IdList', [])]

    @property
    def taxonomy_id(self) -> Union[str, None]:

        id_list = self.taxonomy_ids

        if not id_list:
            return

        elif len(id_list) > 1:
            return

        return str(id_list[0])

    @property
    def assembly_id(self):
        return self._genome_record.get('AssemblyID')

    @property
    def assembly_accession(self):
        return str(self._assembly_record.get('AssemblyAccession', ''))

    @property
    def assembly_name(self):
        return str(self._assembly_record.get('AssemblyName', ''))

    @property
    def genbank_ftp(self):
        return str(self._assembly_record.get('FtpPath_GenBank', '/'))

    @property
    def refseq_ftp(self):
        return str(self._assembly_record.get('FtpPath_RefSeq', '/'))

    @property
    def rna_ftp(self):
        file_name = self.genbank_ftp.split('/')[-1]
        file_name = f'{file_name}_rna_from_genomic.fna.gz'
        return f'{self.genbank_ftp}/{file_name}'

    @property
    def rna_file_name(self):
        file_name = self.genbank_ftp.split('/')[-1]
        return f'{self.name}_{self.taxonomy_id}_{file_name}.fna.gz'

    def has_valid_taxonomy(self):

        if self.taxonomy_id:
            return True

        return False

    def has_valid_genome(self):
        if self.assembly_id:
            return True

        return False

    def has_valid_assembly(self):
        if self.assembly_accession:
            return True

        return False

    def has_valid_ftp(self):

        url = self.genbank_ftp.replace('ftp://', 'https://')

        ftp_response = _get(url)

        if ftp_response.status_code == 200:
            return True

        return False

    def is_valid(self):

        tax = self.has_valid_taxonomy()
        genome = self.has_valid_genome()
        assembly = self.has_valid_assembly()
        ftp = self.has_valid_ftp()

        if tax and genome and assembly and ftp:
            return True

        return False

    def get_rna_from_ftp(self, workdir: str):

        if self.has_valid_ftp():

            try:

                with closing(urllib_request.urlopen(self.rna_ftp)) as request:

                    file_path = os.path.join(workdir, self.rna_file_name)

                    with open(file_path, 'wb') as file:
                        shutil.copyfileobj(request, file)

                self._status = True

            except BaseException:

                pass

        return self._status

    def to_dataframe(self):

        description = {'bigg_id': self.bigg_id,
                       'name': self.name,
                       'taxonomy_id': self.taxonomy_id,
                       'taxonomy_ids': self.taxonomy_ids,
                       'assembly_id': self.assembly_id,
                       'assembly_accession': self.assembly_accession,
                       'assembly_name': self.assembly_name,
                       'refseq': self.refseq_ftp,
                       'rna_ftp': self.rna_ftp,
                       'rna_file': self.rna_file_name,
                       'status': self.status}

        return pd.DataFrame.from_dict(description)


def read_bigg_models(workdir: str) -> dict:
    file_path = os.path.join(workdir, 'bigg_models.json')

    with open(file_path, 'r') as file:
        return json.load(file)


def get_bigg_models(url: str = 'http://bigg.ucsd.edu/api/v2/models',
                    workdir: str = None,
                    by_request: bool = True,
                    verbose: bool = False) -> dict:
    res = {}

    if by_request:

        if verbose:
            print('### Get BiGG models using the REST API ###')

        response = _get(url)

        if response.status_code == 200:
            res = response.json()

            if 'results' in res:
                res['results'] = tuple(res['results'])

        if workdir:
            file_path = os.path.join(workdir, 'bigg_models.json')

            with open(file_path, 'w') as file:
                json.dump(res, file)

    else:

        if verbose:
            print('### Reading local BiGG models ###')

        if not workdir:
            raise ValueError('If by_request is false, please provide a valid workdir')

        res = read_bigg_models(workdir)

    return res


def get_organisms(models_dict: dict,
                  organisms_filter: List[str] = None,
                  verbose: bool = False) -> List[Organism]:
    if verbose:
        print('### Parsing BiGG organisms ###')

    res = []

    models_dict = models_dict.get('results', [])

    if organisms_filter:

        filtered = []
        for model_dict in models_dict:
            if model_dict.get('bigg_id', '') in organisms_filter:
                filtered.append(model_dict)

    else:
        filtered = models_dict

    for model_dict in filtered:
        res.append(Organism.from_dict(model_dict))

    return res


def get_organisms_tax_id(organisms: List[Organism],
                         interval: float = 0.2,
                         verbose: bool = False,
                         workdir: str = None):
    if verbose:
        print('### Get organisms taxonomy identifiers ###')

    for organism in organisms:

        if verbose:
            print(f'Taxonomy identifier for {organism.bigg_id}')

        # ncbi api allows 3 requests per second, 10 if you have an api key
        time.sleep(interval)

        try:
            handle = Entrez.esearch(db='taxonomy', term=organism.name, email=ENTREZ_E_MAIL, api_key=ENTREZ_API_KEY)
            record = Entrez.read(handle)
            handle.close()

            organism.taxonomy_record = record

            if verbose:
                print(f'Taxonomy identifier for {organism.bigg_id} was found: {organism.taxonomy_id}')

        except IOError:

            if verbose:
                print(f'Could not find taxonomy identifier for {organism.bigg_id}')


def get_organisms_assembly_id(organisms: List[Organism],
                              interval: float = 0.2,
                              verbose: bool = False):
    if verbose:
        print('### Get organisms assembly identifiers ###')

    for organism in organisms:

        if verbose:
            print(f'Genome record for {organism.taxonomy_id}')

        if organism.has_valid_taxonomy():

            time.sleep(interval)

            term = f'txid{organism.taxonomy_id}[Organism:noexp]'

            try:
                handle = Entrez.esearch(db='genome', term=term, email=ENTREZ_E_MAIL, api_key=ENTREZ_API_KEY)
                record = Entrez.read(handle)
                handle.close()

                ids = record.get('IdList', [])

                if len(ids) == 1:

                    genome_id = ids[0]

                    handle2 = Entrez.esummary(db='genome', id=str(genome_id),
                                              email=ENTREZ_E_MAIL, api_key=ENTREZ_API_KEY)
                    record2 = Entrez.read(handle2)
                    handle2.close()

                    organism.genome_record = record2[0]

                    if verbose:
                        print(f'Assembly identifier for {organism.bigg_id} was found: {organism.assembly_id}')

                else:
                    if verbose:
                        print(f'Could not find genome record for {ids}')

            except IOError:
                if verbose:
                    print(f'Could not find genome identifier for {term}')


def get_organisms_assembly_record(organisms: List[Organism],
                                  interval: float = 0.2,
                                  verbose: bool = False):
    if verbose:
        print('### Get organisms assembly record ###')

    for organism in organisms:

        if verbose:
            print(f'Assembly record for {organism.bigg_id}')

        if organism.has_valid_genome():

            time.sleep(interval)

            try:
                handle = Entrez.esummary(db='assembly', id=organism.assembly_id,
                                         email=ENTREZ_E_MAIL, api_key=ENTREZ_API_KEY)
                record = Entrez.read(handle)
                handle.close()

                if record:

                    assembly_records = record['DocumentSummarySet']['DocumentSummary']

                    if assembly_records:
                        organism.assembly_record = assembly_records[0]

                if verbose:
                    print(f'Assembly record for {organism.bigg_id} was found')

            except IOError:
                if verbose:
                    print(f'Could not find assembly record for {organism.bigg_id}')


def get_organisms_rna(organisms: List[Organism],
                      workdir: str,
                      interval: float = 0.3,
                      verbose: bool = False):
    if verbose:
        print('### Get organisms rna ###')

    gz_dir = os.path.join(workdir, 'gz_fnas')
    if not os.path.exists(gz_dir):
        os.makedirs(gz_dir)

    ftps = []

    for organism in organisms:

        if organism.rna_ftp in ftps:
            continue

        time.sleep(interval)

        status = organism.get_rna_from_ftp(workdir=gz_dir)

        ftps.append(organism.rna_ftp)

        if verbose:

            if status:

                print(f'RNA for {organism.bigg_id} was obtained')

            else:
                print(f'Failed to obtain RNA for {organism.bigg_id}')


def create_report(organisms: List[Organism],
                  workdir: str,
                  verbose: bool = False):

    if verbose:
        print('### Create report ###')

    dfs = []

    for organism in organisms:
        df = organism.to_dataframe()

        dfs.append(df)

    final_df = pd.concat(dfs)

    file_path = os.path.join(workdir, 'report.xlsx')
    final_df.to_excel(file_path, index=False)


def unpacking_gz(workdir: str,
                 verbose: bool = False):

    if verbose:
        print('### Unpacking gunzip ###')

    gz_dir = os.path.join(workdir, 'gz_fnas')

    if os.path.exists(gz_dir):

        files = os.listdir(gz_dir)

        fna_dir = os.path.join(workdir, 'fnas')
        if not os.path.exists(fna_dir):
            os.makedirs(fna_dir)

        for file in files:

            if file.endswith('.fna.gz'):
                file_path = os.path.join(gz_dir, file)

                with gzip.open(file_path, 'rb') as f:
                    fna_path = os.path.join(fna_dir, file.replace('.gz', ''))

                    with open(fna_path, 'wb') as f_out:
                        shutil.copyfileobj(f, f_out)


def parse_file_name_to_description(file_name):
    # f'{self.name}_{self.taxonomy_id}_{file_name}.fna'

    file_name = file_name.replace('.fna', '')
    names = file_name.split('_')

    name = names[0]

    name = name.replace(' ', '_')

    return name, file_name


def to_compile(workdir: str,
               verbose: bool = False):

    if verbose:
        print('### Compiling organisms rna sequences ###')

    fna_dir = os.path.join(workdir, 'fnas')

    if os.path.exists(fna_dir):
        files = os.listdir(fna_dir)

        rrna_dir = os.path.join(workdir, 'rrna_fna')
        if not os.path.exists(rrna_dir):
            os.makedirs(rrna_dir)

        file_out = os.path.join(rrna_dir, 'rrnas.fna')

        with open(file_out, 'w') as f_out:

            for file in files:

                sequence = Seq('')

                file_path = os.path.join(fna_dir, file)

                for seq_record in SeqIO.parse(file_path, 'fasta'):
                    if '18S' in seq_record.description.upper():
                        sequence += seq_record.seq

                    elif '16S' in seq_record.description.upper():
                        sequence += seq_record.seq

                if len(sequence) > 0:

                    name, description = parse_file_name_to_description(file)

                    record = SeqRecord(sequence,
                                       id=name,
                                       description=description)

                    SeqIO.write(record, f_out, 'fasta')


def main(workdir: str = None,
         rnas: bool = True,
         unpack: bool = True,
         compiling: bool = True,
         verbose: bool = False,
         **kwargs):

    if rnas:

        url = kwargs.get('url', 'http://bigg.ucsd.edu/api/v2/models')
        by_request = kwargs.get('by_request', True)

        models = get_bigg_models(url=url,
                                 workdir=workdir,
                                 by_request=by_request,
                                 verbose=verbose)

        organisms_filter = kwargs.get('organisms_filter', None)
        bigg_organisms = get_organisms(models,
                                       organisms_filter=organisms_filter,
                                       verbose=verbose)

        interval = kwargs.get('interval', 0.2)
        get_organisms_tax_id(bigg_organisms,
                             interval=interval,
                             verbose=verbose,
                             workdir=workdir)

        get_organisms_assembly_id(bigg_organisms,
                                  interval=interval,
                                  verbose=verbose)

        get_organisms_assembly_record(bigg_organisms,
                                      interval=interval,
                                      verbose=verbose)

        if not workdir:
            workdir = os.path.join(os.getcwd(), 'organisms_rna')

            if not os.path.exists(workdir):
                os.makedirs(workdir)

        interval = kwargs.get('interval', 0.2)
        interval += 0.1
        get_organisms_rna(bigg_organisms,
                          workdir=workdir,
                          interval=interval,
                          verbose=verbose)

        create_report(bigg_organisms, workdir=workdir, verbose=verbose)

    if unpack:

        if not workdir:
            workdir = os.path.join(os.getcwd(), 'organisms_rna')

            if not os.path.exists(workdir):
                os.makedirs(workdir)

        unpacking_gz(workdir, verbose=verbose)

    if compiling:

        if not workdir:
            workdir = os.path.join(os.getcwd(), 'organisms_rna')

            if not os.path.exists(workdir):
                os.makedirs(workdir)

        to_compile(workdir, verbose=verbose)


if __name__ == '__main__':

    directory = os.path.join(os.getcwd(), 'organisms_rna')

    if not os.path.exists(directory):
        os.makedirs(directory)

    main(workdir=directory, rnas=False,unpack=False, verbose=True)

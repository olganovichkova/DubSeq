import os
import random
import math
import json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import poisson
from scipy.optimize import nnls
from statistics import median
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge


class BarseqLayoutItem:
    def __init__(self, itnum, item_type, experiment_condition):
        self.__itnum = itnum
        self.__item_type = item_type
        self.__experiment_condition = experiment_condition

    @property
    def itnum(self):
        return self.__itnum

    @property
    def item_type(self):
        return self.__item_type

    @property
    def experiment_condition(self):
        return self.__experiment_condition


class BarseqLayout:
    def __init__(self, layout_file_name):
        self.__layout_file_name = layout_file_name
        self.__df = None
        self.__load()

        # Check if the layout has time zero items
        if len(self.time_zero_items) == 0:
            raise ValueError(
                'No time zero experiments were found in the layout')

    def save(self, fname):
        with open(fname, 'w') as f:
            f.write('\t'.join((
                'itnum',
                'type',
                'name'
            )))
            f.write('\n')

            for item in self.all_items:
                f.write(
                    '\t'.join(
                        str(x) for x in [
                            item.itnum,
                            item.item_type,
                            item.experiment_condition
                        ]
                    )
                )
                f.write('\n')

    def __load(self):
        self.__df = pd.read_csv(self.__layout_file_name, sep='\t')

    @property
    def layout_file_name(self):
        return self.__layout_file_name

    @property
    def time_zero_items(self):
        df = self.__df
        return self.__to_items(df[df.type == 'Time0'])

    @property
    def lb_items(self):
        df = self.__df
        return self.__to_items(df[df.type == 'LB'])

    @property
    def stress_items(self):
        df = self.__df
        return self.__to_items(df[df.type == 'stress'])

    @property
    def non_time_zero_items(self):
        df = self.__df
        return self.__to_items(df[df.type != 'Time0'])

    @property
    def all_items(self):
        return self.__to_items(self.__df)

    def __to_items(self, df):
        items = []
        for _, row in df.iterrows():
            # print(row['name'])
            # items.append(BarseqLayoutItem(row.itnum, row.type, row.name))
            items.append(BarseqLayoutItem(row.itnum, row.type, row['name']))
        return items

    @property
    def experiment_types(self):
        return self.__df.type.unique()


class BpagItem:
    __slots__ = ['barcode_up', 'barcode_dn',
                 'bpair_read_count', 'up_read_count', 'dn_read_count',
                 'contig_id',  'pos_from', 'pos_to']

    def __init__(self, barcode_up, barcode_dn, bpair_read_count,
                 up_read_count, dn_read_count,
                 contig_id, pos_from, pos_to):
        self.barcode_up = barcode_up
        self.barcode_dn = barcode_dn
        self.bpair_read_count = bpair_read_count
        self.up_read_count = up_read_count
        self.dn_read_count = dn_read_count
        self.contig_id = contig_id
        self.pos_from = pos_from
        self.pos_to = pos_to


class BpagSet:
    def __init__(self, blag_file_name):
        self.__blag_file_name = blag_file_name
        self.__items = []
        self.__up_barcode_2_item = {}
        self.__load()

    @property
    def blag_file_name(self):
        return self.__blag_file_name

    @property
    def size(self):
        return len(self.__items)

    def get_item(self, index):
        return self.__items[index]

    def find_up_item(self, barcode_up):
        return self.__up_barcode_2_item.get(barcode_up)

    def __load(self):
        df = pd.read_csv(self.__blag_file_name, sep='\t')
        for _, row in df.iterrows():
            if row.recommended == '+':
                item = BpagItem(
                    row.barcode_up,
                    row.barcode_dn,
                    row.bpair_read_count,
                    row.up_read_count,
                    row.dn_read_count,
                    row.up_contig_id,
                    row.pos_from,
                    row.pos_to
                )
                self.__items.append(item)
                self.__up_barcode_2_item[item.barcode_up] = item


class TimeZeroItem:
    def __init__(self, barcode, time0_experiments_count):
        self.__barcode = barcode
        self.__read_counts = [0] * time0_experiments_count

    @property
    def barcode(self):
        return self.__barcode

    @property
    def total_read_count(self):
        return sum(self.__read_counts)

    @property
    def max_read_count(self):
        return max(self.__read_counts)

    def set_read_count(self, experiment_index, count):
        self.__read_counts[experiment_index] = count


class TimeZeroSet:
    def __init__(self, bpag_set, barseq_layout, barseq_dir):
        self.__time0_itnums = []
        for item in barseq_layout.time_zero_items:
            self.__time0_itnums.append(item.itnum)

        self.__barcode_2_item = {}
        self.__items = []
        self.__load(bpag_set, barseq_dir)

    def filter_items(self, good_item_method):
        for i in range(self.size)[::-1]:
            if not good_item_method(self.__items[i]):
                barcode = self.__items[i].barcode
                del self.__barcode_2_item[barcode]
                del self.__items[i]

    @property
    def size(self):
        return len(self.__items)

    def __load(self, bpag_set, barseq_dir):
        for experiment_index, itnum in enumerate(self.__time0_itnums):
            bstat_fname = self.__get_bstat_file(itnum, barseq_dir)
            if not bstat_fname:
                raise ValueError(
                    'Can not find bstat file for itnum %s in %s directory' % (itnum, barseq_dir))

            df = pd.read_csv(bstat_fname, sep='\t')
            for _, row in df.iterrows():

                if row.recommnended != '+':
                    continue

                if not bpag_set.find_up_item(row.barcode):
                    continue

                self.__register_read_count(
                    row.barcode, experiment_index, int(row.reads_count))

    def __get_bstat_file(self, itnum, barseq_dir):
        file_path = None
        for file_name in os.listdir(barseq_dir):
            if not file_name.endswith('.bstat.tsv'):
                continue
            if '_' + itnum + '_' in file_name:
                file_path = os.path.join(barseq_dir, file_name)
                break

        return file_path

    @property
    def experiment_count(self):
        return len(self.__time0_itnums)

    def __register_read_count(self, barcode, exp_index, read_count):
        t0_item = self.__barcode_2_item.get(barcode)
        if not t0_item:
            t0_item = TimeZeroItem(barcode, self.experiment_count)
            self.__barcode_2_item[barcode] = t0_item

        if t0_item:
            t0_item.set_read_count(exp_index, read_count)


class Fitness:

    # GFF_FILE = None
    # BPAG_FILE = None

    SCORE_TYPE_MEAN = 0
    SCORE_TYPE_NNLS = 1
    SCORE_TYPE_C_NNLS = 2
    SCORE_TYPE_RIDGE = 3
    SCORE_TYPE_LASSO = 4
    SCORE_TYPE_ELASTIC_NET = 5
    SCORE_TYPE_NAMES = ['mean', 'nnls', 'cnnls', 'ridge', 'lasso', 'enet']

    MIN_TIME0_READ_COUNT = 10
    CONDITIONS = {}
    BARCODE_2_INDEX = {}
    BARCODE_COUNTS = []
    BARCODE_INDICES = []
    BARCODE_REPLICATES = []

    GENES = []
    GENOME_SEGMENTS = []

# Ecoli
# ITNUM_INDEX = 3

# scoreType = SCORE_TYPE_RIDGE
# GENE_SCORES_DIR = ROOT_DIR + '/narrative_data/data/barseq/gene_scores/ridge_1e-3'
    RIDGE_PARAM_ALPHA = 0.001


# scoreType = SCORE_TYPE_LASSO
# GENE_SCORES_DIR = ROOT_DIR + '/narrative_data/data/barseq/gene_scores/lasso_5e-5'
    LASSO_PARAM_ALPHA = 0.00005

# scoreType = SCORE_TYPE_ELASTIC_NET
    ELASTIC_NET_PARAM_A = 0.0005
    ELASTIC_NET_PARAM_B = 0.0001
# GENE_SCORES_DIR = ROOT_DIR + '/narrative_data/data/barseq/gene_scores/enet_4e-5_1e-7'


# CONDITIONS_CONFIG = {
#     'FEBA_133': {
#         'barseqDir': ROOT_DIR + '/narrative_data/data/barseq/FEBA_BS_133/',
#         'conditionFile': ROOT_DIR + '/narrative_data/data/barseq/FEBA_133_conditions.txt',
#         'df': None
#     },
#     'FEBA_134': {
#         'barseqDir': ROOT_DIR + '/narrative_data/data/barseq/FEBA_BS_134/',
#         'conditionFile': ROOT_DIR + '/narrative_data/data/barseq/FEBA_134_conditions.txt',
#         'df': None
#     },
#     'FEBA_136': {
#         'barseqDir': ROOT_DIR + '/narrative_data/data/barseq/FEBA_BS_136/',
#         'conditionFile': ROOT_DIR + '/narrative_data/data/barseq/FEBA_136_conditions.txt',
#         'df': None
#     }
# }

    #######################
    # Init
    #######################

    # @staticmethod
    # def init():
    #     for expId in Fitness.CONDITIONS_CONFIG:
    #         exp = Fitness.CONDITIONS_CONFIG[expId]
    #         exp['df'] = pd.read_csv(exp['conditionFile'],
    #                                 quotechar='"', delimiter='\t')

    @staticmethod
    def init(barseq_layout, barseq_dir, bpag_fname, genes_gff_fname=None):
        Fitness.initConditions(barseq_layout, barseq_dir)

        t0Indeces = Fitness.getTimeZeroIndeces()

        # Load data
        Fitness.loadBPAG(bpag_fname)
        Fitness.loadCounts()
        Fitness.buildREF_TIME0(t0Indeces)
        Fitness.cleanBARCODE_COUNTS(t0Indeces)

        if genes_gff_fname:
            Fitness.loadGenes(genes_gff_fname)
            Fitness.associateGenesWithBarcodes()
            # cleanGENES()
            Fitness.buildGENOME_SEGMENTS()

    @staticmethod
    def initConditions(barseq_layout, barseq_dir):

        for index, item in enumerate(barseq_layout.all_items):
            Fitness.CONDITIONS[item.itnum] = {
                "index": index,
                "type": item.item_type,
                "desc": item.experiment_condition}

        for file_name in os.listdir(barseq_dir):
            if not file_name.endswith('.bstat.tsv'):
                continue

            vals = file_name.split("_")
            for val in vals:
                if val in Fitness.CONDITIONS:
                    itnum = val
                    Fitness.CONDITIONS[itnum]['file'] = os.path.join(
                        barseq_dir, file_name)
                    break

        # @staticmethod
        # def initConditions(expId):
        #     Fitness.CONDITIONS.clear()

        #     df = Fitness.CONDITIONS_CONFIG[expId]['df']
        #     for index, row in df.iterrows():
        #         itnum = row['itnum']
        #         itype = row['type']
        #         itdesc = row['name']
        #         Fitness.CONDITIONS[itnum] = {
        #             "index": index, "type": itype, "desc": itdesc}

        #     barseqDir = Fitness.CONDITIONS_CONFIG[expId]['barseqDir']
        #     for file in os.listdir(barseqDir):
        #         if not file.endswith('.b1'):
        #             continue
        #         vals = file.split("_")
        #         if len(vals) > ITNUM_INDEX:
        #             itnum = vals[ITNUM_INDEX]
        #             if itnum in Fitness.CONDITIONS:
        #                 Fitness.CONDITIONS[itnum]['file'] = barseqDir + file

        #######################
        # Gettters
        #######################

    @staticmethod
    def getTimeZeroIndeces():
        indeces = []
        for itnum in Fitness.CONDITIONS:
            it = Fitness.CONDITIONS[itnum]
            if it['type'] == 'Time0':
                indeces.append(it['index'])
        return indeces

    @staticmethod
    def getSample(itIndex):
        sample = []
        for bIndex in Fitness.BARCODE_INDICES:
            row = Fitness.BARCODE_COUNTS[bIndex]
            count = row['counts'][itIndex]
            sample.append(count)
        return sample

    @staticmethod
    def getRefTime0Sample():
        sample = []
        for bIndex in Fitness.BARCODE_INDICES:
            row = Fitness.BARCODE_COUNTS[bIndex]
            count = row['time0']
            sample.append(count)
        return sample

    @staticmethod
    def getTotalCount(itIndex):
        total = 0
        for bIndex in Fitness.BARCODE_INDICES:
            row = Fitness.BARCODE_COUNTS[bIndex]
            total += row['counts'][itIndex]
        return total

    @staticmethod
    def getItNum(itIndex):
        for itnum in Fitness.CONDITIONS:
            it = Fitness.CONDITIONS[itnum]
            if it['index'] == itIndex:
                return itnum
        return None

    #######################
    # Exporters
    #######################

    @staticmethod
    def save_gscore_base(fname):
        with open(fname, 'w') as f:
            f.write(
                '\t'.join((
                    'gene_index',
                    'covering_fragment_count',
                    'name',
                    'locus_tag',
                    'gene_type',
                    'contig_id',
                    'pos_from',
                    'pos_to',
                    'strand',
                    'product',
                    'note',
                    'description',
                    'barcodes'
                ))
            )
            f.write('\n')
            for gene_index, gene in enumerate(Fitness.GENES):
                vals = [
                    gene_index,
                    len(gene['barcodeIndeces']),
                    gene['name'],
                    gene['locusTag'],
                    gene['geneType'],
                    gene['contigId'],
                    gene['posFrom'],
                    gene['posTo'],
                    gene['strand'],
                    gene['product'],
                    gene['note'],
                    gene['description'],
                    ','.join(Fitness.BARCODE_COUNTS[i]['barcode']
                             for i in gene['barcodeIndeces'])
                ]
                f.write('\t'.join(str(x) for x in vals))
                f.write('\n')

    @staticmethod
    def save_fscore_base(fname):
        with open(fname, 'w') as f:
            f.write('\t'.join((
                'barcode',
                'contig_id',
                'pos_from',
                'pos_to',
                't0_count',
                't0_reads_avg',
                't0_reads_total',
                't0_reads',
                't0_itnums'
            )))
            f.write('\n')

            time_zero_indeces = Fitness.getTimeZeroIndeces()
            for item in Fitness.BARCODE_COUNTS:

                time_zero_vals = [item['counts'][i] for i in time_zero_indeces]
                time_zero_itnums = [Fitness.getItNum(
                    i) for i in time_zero_indeces]

                vals = [
                    item['barcode'],
                    item['contigId'],
                    item['posFrom'],
                    item['posTo'],
                    len(time_zero_vals),
                    sum(time_zero_vals) / float(len(time_zero_vals)),
                    sum(time_zero_vals),
                    ','.join(str(x) for x in time_zero_vals),
                    ','.join(str(x) for x in time_zero_itnums)
                ]

                f.write('\t'.join(str(x) for x in vals))
                f.write('\n')

    @staticmethod
    def save_fscores(score_fname, fs, ss, ts):
        with open(score_fname, 'w') as f:
            f.write('%s\n' % '\t'.join(
                ['barcode', 'score', 'stress_read_count', 't0_total_read_count']))
            for index, score in enumerate(fs):
                barcode = Fitness.BARCODE_COUNTS[index]['barcode']
                f.write('%s\n' % '\t'.join(str(x)
                                           for x in [barcode, score, ss[index], ts[index]]))

    @staticmethod
    def save_gscores(score_fname, score_types, gss):
        with open(score_fname, 'w') as f:
            column_names = ['index', 'gene_name', 'locus_tag']
            for score_type in score_types:
                column_names.append(
                    'score_' + Fitness.SCORE_TYPE_NAMES[score_type])

            f.write('\t'.join(column_names) + '\n')
            for index, gene in enumerate(Fitness.GENES):

                vals = [index, gene['name'], gene['locusTag']]
                for gs in gss:
                    vals.append(gs[index])

                f.write('\t'.join(str(x) for x in vals) + '\n')
                # '%s\t%s\t%s\n' % (index, name, score))

    # @staticmethod
    # def exportGeneScores(fileName, geneScores):
    #     with open(fileName, 'w') as f:
    #         for index in range(0, len(Fitness.GENES)):
    #             score = geneScores[index]
    #             gene = Fitness.GENES[index]
    #             f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
    #                 score, gene['locusTag'], gene['name'], gene['strand'], gene['posFrom'], gene['posTo'], gene['product']))

    # @staticmethod
    # def exportRegionData(fileName, sampleCounts, fitnessScores):
    #     with open(fileName, 'w') as f:
    #         for index in range(0, len(Fitness.BARCODE_COUNTS)):
    #             row = Fitness.BARCODE_COUNTS[index]
    #             f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % (row['pairBarcodeUp'], row['pairBarcodeDn'], row['posFrom'], row['posTo'], 0,
    #                                                                     row['barcode'], sampleCounts[index], Fitness.REF_TIME0[index], 0, fitnessScores[index]))

    # def exportProjectDescriptor(name, datasetName, scoreMethod):

    #     cnds = []

    #     for i in range(0, len(CONDITIONS_LIST) / 3):
    #         itNum = CONDITIONS_LIST[i * 3]
    #         itDesc = CONDITIONS_LIST[i * 3 + 2]
    #         cnds.append({
    #             'itNumber': itNum,
    #             'name': itNum + ": " + itDesc
    #         })

    #     descriptor = {
    #         'name': name,
    #         'datasetName': datasetName,
    #         'scoreMethod': scoreMethod,
    #         'conditions': cnds
    #     }
    #     return json.dumps(descriptor)

    #######################
    # Loaders
    #######################

    @staticmethod
    def loadBPAG(bpag_fname):
        print("loadBPAG:...",  end='', flush=True)
        del Fitness.BARCODE_COUNTS[:]
        df = pd.read_csv(bpag_fname, sep='\t')

        for _, row in df.iterrows():
            if row.recommended != '+':
                continue
            Fitness.BARCODE_COUNTS.append({
                'pairBarcodeUp': row.barcode_up,
                'pairBarcodeDn': row.barcode_dn,

                "barcode": row.barcode_up,
                "contigId": row.up_contig_id,
                "posFrom": row.pos_from,
                "posTo": row.pos_end,
                "time0": 0,
                "counts": [0] * len(Fitness.CONDITIONS)
            })

            # with open(bpag_fname, 'r') as f:
            #     for line in f:
            #         vals = line.split("\t")
            #         barcodes = vals[0].split("=")
            #         barcode = barcodes[1]
            #         Fitness.BARCODE_COUNTS.append({
            #             'pairBarcodeUp': barcodes[0],
            #             'pairBarcodeDn': barcodes[1],

            #             "barcode": barcode,
            #             "contigId": vals[3],
            #             "posFrom": int(vals[4]),
            #             "posTo": int(vals[5]),
            #             "time0": 0,
            #             "counts": [0] * len(Fitness.CONDITIONS)
            #         })
        Fitness.updateBARCODE_INDICES()
        print("Done!")

    @staticmethod
    def updateBARCODE_INDICES():
        Fitness.BARCODE_2_INDEX.clear()
        del Fitness.BARCODE_INDICES[:]
        del Fitness.BARCODE_REPLICATES[:]
        for index, br in enumerate(Fitness.BARCODE_COUNTS):
            Fitness.BARCODE_2_INDEX[br['barcode']] = index
            Fitness.BARCODE_INDICES.append(index)
            Fitness.BARCODE_REPLICATES.append(1)
    #     print "\t from updateBARCODE_INDICES", len(BARCODE_2_INDEX), len(BARCODE_INDICES), len(BARCODE_REPLICATES)

    @staticmethod
    def loadCounts():
        print("loadCounts:...", end='', flush=True)
        for itnum in Fitness.CONDITIONS:
            it = Fitness.CONDITIONS[itnum]
            print('Doing file: ', it['file'])
            df = pd.read_csv(it['file'], sep='\t')
            for _, row in df.iterrows():
                if row.sim_recommended != '+':
                    continue
                barcode = row.barcode
                count = row.reads_count
                if barcode in Fitness.BARCODE_2_INDEX:
                    barcodeIndex = Fitness.BARCODE_2_INDEX[barcode]
                    itIndex = it['index']
                    Fitness.BARCODE_COUNTS[barcodeIndex]['counts'][itIndex] = count

            # with open(it['file'], 'r') as f:
            #     for line in f:
            #         vals = line.split('\t')
            #         barcode = vals[0]
            #         count = int(vals[1])
            #         if barcode in Fitness.BARCODE_2_INDEX:
            #             barcodeIndex = Fitness.BARCODE_2_INDEX[barcode]
            #             itIndex = it['index']
            #             Fitness.BARCODE_COUNTS[barcodeIndex]['counts'][itIndex] = count
        print("Done!")

    @staticmethod
    def cleanBARCODE_COUNTS(time0Indeces):
        len0 = len(Fitness.BARCODE_COUNTS)
        delCount = 0
        for i in range(len0)[::-1]:
            counts = Fitness.BARCODE_COUNTS[i]['counts']
            hasData = False
            for t0Index in time0Indeces:
                if counts[t0Index] >= Fitness.MIN_TIME0_READ_COUNT:
                    hasData = True
                    break
            if not hasData:
                del Fitness.BARCODE_COUNTS[i]
                delCount += 1

        Fitness.updateBARCODE_INDICES()
        # print "\tLen before = ", len0
        # print "\tLen after = ", len(BARCODE_COUNTS)
        # print "\tDel count = ", delCount
        # print "cleanBARCODE_COUNTS: Done!"

    @staticmethod
    def loadGenes(genes_gff_fname):
        del Fitness.GENES[:]

        # First read all features that has "Parent" property and hash them
        id2features = {}
        with open(genes_gff_fname, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                vals = line.split('\t')

                f_contig = vals[0]
                f_pos_from = int(vals[3])
                f_pos_to = int(vals[4])
                f_strand = vals[6]
                f_description = vals[8].strip()

                f_parent = None
                f_name = ""
                f_product = ""
                f_note = ""
                f_pseudo = False
                for dval in f_description.split(";"):
                    if dval.startswith("Parent="):
                        f_parent = dval[len("Parent="):].strip()
                    elif dval.startswith("gene="):
                        f_name = dval[len("gene="):].strip()
                    elif dval.startswith("product="):
                        f_product = dval[len("product="):].strip()
                    elif dval.startswith("Note="):
                        note = dval[len("Note="):].strip()
                    elif 'pseudo=true' in dval:
                        f_pseudo = True

                if f_parent:
                    features = id2features.get(f_parent)
                    if not features:
                        features = []
                        id2features[f_parent] = features
                    features.append({
                        'gene_type': vals[2],
                        'gene_name': f_name,
                        'contig': f_contig,
                        'pos_from': f_pos_from,
                        'pos_to': f_pos_to,
                        'strand': f_strand,
                        'pseudo': f_pseudo,
                        'product': f_product,
                        'note': f_note,
                        'description': f_description
                    })

        # Now read all "gene" features and collect of children
        with open(genes_gff_fname, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                vals = line.split('\t')
                if vals[2] == 'gene':
                    gene_contig = vals[0]
                    gene_pos_from = int(vals[3])
                    gene_pos_to = int(vals[4])
                    gene_strand = vals[6]
                    gene_description = vals[8].strip()

                    gene_locus_tag = None
                    gene_id = None
                    for term in vals[8].split(';'):
                        (key, value) = term.split('=')
                        if key == 'locus_tag':
                            gene_locus_tag = value.strip()
                        elif key == 'ID':
                            gene_id = value.strip()

                    if not gene_id:
                        continue

                    features = id2features.get(gene_id)
                    if not features:
                        continue

                    # build features related to this gene and locations are correct
                    gene_features = []
                    for f in features:
                        if f['contig'] != gene_contig:
                            continue
                        if f['strand'] != gene_strand:
                            continue
                        if f['pos_from'] < gene_pos_from:
                            continue
                        if f['pos_to'] > gene_pos_to:
                            continue
                        gene_features.append(f)

                    if len(gene_features) == 0:
                        continue

                    # if there are more than one feature, check that the type of feature is the same
                    gene_types = {}
                    for f in gene_features:
                        gene_types[f['gene_type']] = 1
                    if len(gene_types) > 1:
                        raise ValueError(
                            "More than one gene type for a given gene: " + gene_id)

                    f = gene_features[0]

                    Fitness.GENES.append({
                        'contigId': f['contig'],
                        'geneType': f['gene_type'],
                        'posFrom': f['pos_from'],
                        'posTo': f['pos_to'],
                        'strand': f['strand'],
                        'name': f['gene_name'],
                        'product': f['product'],
                        'locusTag': gene_locus_tag,
                        'note': f['note'],
                        'description': f['description'],
                        'index': 0,
                        'barcodeIndeces': []
                    })

        Fitness.GENES.sort(key=lambda x: x['posFrom'], reverse=False)
        print('Load genes: Done!')

    @staticmethod
    def _loadGenes(genes_gff_fname):
        del Fitness.GENES[:]

        gene_found = False
        gene_contig = None
        gene_pos_from = None
        gene_pos_to = None
        gene_strand = None
        gene_locus_tag = None

        with open(genes_gff_fname, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                vals = line.split('\t')
                gene_type = vals[2]

                if gene_type == 'gene':
                    gene_contig = vals[0]
                    gene_pos_from = int(vals[3])
                    gene_pos_to = int(vals[4])
                    gene_strand = vals[6]
                    gene_found = True
                    gene_locus_tag = None
                    for term in vals[8].split(';'):
                        (key, value) = term.split('=')
                        if key == 'locus_tag':
                            gene_locus_tag = value.strip()

                if gene_type == 'CDS' or 'rna' in gene_type.lower():
                    if gene_found:
                        f_contig = vals[0]
                        f_pos_from = int(vals[3])
                        f_pos_to = int(vals[4])
                        f_strand = vals[6]
                        f_description = vals[8].strip()

                        if f_contig == gene_contig and f_pos_from == gene_pos_from and f_pos_to == gene_pos_to and f_strand == gene_strand:

                            f_locus_tag = gene_locus_tag
                            f_name = ""
                            f_product = ""
                            f_note = ""
                            for dval in f_description.split(";"):
                                if dval.startswith("gene="):
                                    f_name = dval[len("gene="):].strip()
                                elif dval.startswith("product="):
                                    f_product = dval[len("product="):].strip()
                                elif dval.startswith("Note="):
                                    f_note = dval[len("Note="):].strip()

                            Fitness.GENES.append({
                                'contigId': f_contig,
                                'geneType': gene_type,
                                'posFrom': f_pos_from,
                                'posTo': f_pos_to,
                                'strand': f_strand,
                                'name': f_name,
                                'product': f_product,
                                'locusTag': f_locus_tag,
                                'note': f_note,
                                'description': f_description,
                                'index': 0,
                                'barcodeIndeces': []
                            })

                        gene_found = False
        print("loadGenes: Done!")

        # with open(genes_gff_fname, 'r') as f:
        #     for line in f:
        #         if line.startswith("#"):
        #             continue
        #         vals = line.split("\t")
        #         contigid = vals[0]
        #         geneType = vals[2]
        #         posFrom = int(vals[3])
        #         posTo = int(vals[4])
        #         strand = vals[6]
        #         description = vals[8].strip()

        #         name = ""
        #         product = ""
        #         locusTag = ""
        #         note = ""
        #         for dval in description.split(";"):
        #             if dval.startswith("gene="):
        #                 name = dval[len("gene="):].strip()
        #             elif dval.startswith("product="):
        #                 product = dval[len("product="):].strip()
        #             elif dval.startswith("Name="):
        #                 locusTag = dval[len("Name="):].strip()
        #             elif dval.startswith("Note="):
        #                 note = dval[len("Note="):].strip()

        #         if geneType == 'CDS' or "RNA" in geneType:
        #             Fitness.GENES.append({
        #                 'contigId': contigid,
        #                 'geneType': geneType,
        #                 'posFrom': posFrom,
        #                 'posTo': posTo,
        #                 'strand': strand,
        #                 'name': name,
        #                 'product': product,
        #                 'locusTag': locusTag,
        #                 'note': note,
        #                 'description': description,
        #                 'index': 0,
        #                 'barcodeIndeces': []
        #             })

    @staticmethod
    def associateGenesWithBarcodes():
        print("associateGenesWithBarcodes:...", end='', flush=True)
        for bIndex, bcode in enumerate(Fitness.BARCODE_COUNTS):
            for gene in Fitness.GENES:
                if bcode['contigId'] == gene['contigId']:
                    if bcode['posFrom'] <= gene['posFrom'] and bcode['posTo'] >= gene['posTo']:
                        gene['barcodeIndeces'].append(bIndex)

        print("Done!")

    # Delete genes that do not have associated barcodes
    @staticmethod
    def cleanGENES():
        len0 = len(Fitness.GENES)
        delCount = 0
        for i in range(len0)[::-1]:
            gene = Fitness.GENES[i]
            if len(gene['barcodeIndeces']) == 0:
                del Fitness.GENES[i]
                delCount += 1

        # print "\tLen before = ", len0
        # print "\tLen after = ", len(GENES)
        # print "\tDel count = ", delCount
        # print "cleanGENES: Done!"

    #######################
    # Core Builders
    #######################
    @staticmethod
    def buildGENOME_SEGMENTS():
        print("buildGENOME_SEGMENTS:... ", end='', flush=True)
        del Fitness.GENOME_SEGMENTS[:]

        prevGene = None
        geneIndeces = []
        for gIndex, gene in enumerate(Fitness.GENES):
            if prevGene is not None:
                hasSameBarcode = False
                for bIndex1 in prevGene['barcodeIndeces']:
                    for bIndex2 in gene['barcodeIndeces']:
                        if bIndex1 == bIndex2:
                            hasSameBarcode = True
                            break
                if not hasSameBarcode:
                    # new segment
                    if len(geneIndeces) > 0:
                        Fitness.GENOME_SEGMENTS.append({
                            'geneIndeces': geneIndeces
                        })
                        geneIndeces = []
            geneIndeces.append(gIndex)
            prevGene = gene

        # final check
        if len(geneIndeces) > 0:
            Fitness.GENOME_SEGMENTS.append({
                'geneIndeces': geneIndeces
            })
        print('Done!')

    @staticmethod
    def buildREF_TIME0(time0Indeces):
        print("buildREF_TIME0:...", end='', flush=True)
        for bIndex in Fitness.BARCODE_INDICES:
            row = Fitness.BARCODE_COUNTS[bIndex]
            counts = row['counts']
            total = 0
            for t0Index in time0Indeces:
                total += counts[t0Index]
            Fitness.BARCODE_COUNTS[bIndex]['time0'] = total
        print("Done!")

    # No needs to adjust by total...
    @staticmethod
    def buildFitnessScore(sample, sampleT0):
        scores = []
        stotal = sum(sample)
        stotalT0 = sum(sampleT0)

        for index, t in enumerate(sampleT0):
            s = sample[index]
            score = (s + 1.0) / (t + 1.0) * stotalT0 / stotal
            scores.append(score)
        # normalize by median
        scoreMedian = median(scores)
        for index, val in enumerate(scores):
            scores[index] = math.log(val * 1.0 / scoreMedian, 2)

        return scores

    @staticmethod
    def buildGeneScores(fscores, scoreType):
        geneScores = [0] * len(Fitness.GENES)

        for genomeSegment in Fitness.GENOME_SEGMENTS:

            # v - array of fragment scores
            v = []
            # A - 2d array of presence/absence
            A = []

            geneIndeces = list(genomeSegment['geneIndeces'])

            # define total list of barcodeIndeces
            bIndeces = []
            bIndecesHash = {}
            for gIndex in geneIndeces:
                for bIndex in Fitness.GENES[gIndex]['barcodeIndeces']:
                    bIndecesHash[str(bIndex)] = bIndex

            for key in bIndecesHash:
                bIndex = bIndecesHash[key]
                for i in range(Fitness.BARCODE_REPLICATES[bIndex]):
                    bIndeces.append(bIndex)

            # Review gene indeces (remove thouse which are not coverd by bIndeces)
            gCount = len(geneIndeces)
            for i in range(gCount)[::-1]:
                gIndex = geneIndeces[i]
                if len(set(Fitness.GENES[gIndex]['barcodeIndeces']).intersection(bIndeces)) == 0:
                    del geneIndeces[i]

            if len(geneIndeces) > 0 and len(bIndeces) > 0:
                # build a subset of barcode (fragemtn) scores corresponding to bIndeces
                for bIndex in bIndeces:
                    v.append(fscores[bIndex])

                # build matrix of presence/absence
                for gIndex in geneIndeces:
                    row = [0] * len(bIndeces)
                    for vIndex, bIndex1 in enumerate(bIndeces):
                        for bIndex2 in Fitness.GENES[gIndex]['barcodeIndeces']:
                            if bIndex1 == bIndex2:
                                row[vIndex] = 1
                    A.append(row)

                # convert to numpy array
                v = np.array(v)
                A = np.array(A)
                A = A.T

                scores = []
                if scoreType == Fitness.SCORE_TYPE_MEAN:
                    for gIndex in range(A.shape[1]):
                        score = 0
                        n = 0
                        for bIndex in range(A.shape[0]):
                            if A[bIndex, gIndex] == 1:
                                score += v[bIndex]
                                n += 1
                        if n > 0:
                            score /= n
                        scores.append(score)
                elif scoreType == Fitness.SCORE_TYPE_NNLS:
                    x = nnls(A, v)
                    scores = x[0]
                elif scoreType == Fitness.SCORE_TYPE_C_NNLS:
                    x = nnls(A, v)
                    scoresDirect = x[0]

                    x = nnls(A, v * (-1))
                    scoresReverse = x[0]
                    for i, scoreDirect in enumerate(scoresDirect):
                        scoreReverse = -scoresReverse[i]
                        score = 0
                        if scoreDirect != 0 and scoreReverse == 0:
                            score = scoreDirect
                        elif scoreDirect == 0 and scoreReverse != 0:
                            score = scoreReverse
                        scores.append(score)

                elif scoreType == Fitness.SCORE_TYPE_ELASTIC_NET:
                    # alpha = a + b and l1_ratio = a / (a + b)
                    # a * L1 + b * L2

                    alpha = Fitness.ELASTIC_NET_PARAM_A + Fitness.ELASTIC_NET_PARAM_B
                    l1_ratio = Fitness.ELASTIC_NET_PARAM_A / alpha
                    estimator = ElasticNet(
                        alpha=alpha, l1_ratio=l1_ratio, normalize=True)
                    estimator.fit(A, v)
                    scores = estimator.coef_

                elif scoreType == Fitness.SCORE_TYPE_LASSO:
                    estimator = Lasso(
                        alpha=Fitness.LASSO_PARAM_ALPHA, normalize=True)
                    estimator.fit(A, v)
                    scores = estimator.coef_

                elif scoreType == Fitness.SCORE_TYPE_RIDGE:
                    estimator = Ridge(
                        alpha=Fitness.RIDGE_PARAM_ALPHA, normalize=True, solver='lsqr')
                    estimator.fit(A, v)
                    scores = estimator.coef_

                for i, geneScore in enumerate(scores):
                    gIndex = geneIndeces[i]
                    geneScores[gIndex] = geneScore

    #     print 'buildGeneScores: Done!'
        return geneScores

    @staticmethod
    def buildPoissonNoisedSample(sample):
        pSample = []
        for count in sample:
            pcount = poisson.rvs(count, size=1)[0] if count > 0 else 0
            pSample.append(pcount)
        return pSample

    @staticmethod
    def bootstrapSampleReadCounts(sample):
        bootSample = [0] * len(sample)
        total = sum(sample)
        indices = list(range(len(sample)))
        probs = [0] * len(sample)
        for i, val in enumerate(sample):
            probs[i] = (val) * 1.0 / total

        bIndeces = np.random.choice(indices, total, replace=True, p=probs)
        for bIndex in bIndeces:
            bootSample[bIndex] += 1

        return bootSample

    @staticmethod
    def bootstrapBARCODE_INDICES():
        barcodesNumber = len(Fitness.BARCODE_COUNTS)
        indices = list(range(barcodesNumber))
        probs = [0] * barcodesNumber
        for i, val in enumerate(indices):
            probs[i] = 1.0 / barcodesNumber

        BARCODE_INDICES = np.random.choice(
            indices, barcodesNumber, replace=True, p=probs)
        for i in range(len(Fitness.BARCODE_REPLICATES)):
            Fitness.BARCODE_REPLICATES[i] = 0

        for bIndex in BARCODE_INDICES:
            Fitness.BARCODE_REPLICATES[bIndex] += 1

    @staticmethod
    def buildNoisedGeneScores(nCycles, sample, scoreType, doBootstrapIndeces, doBootstrapReadCounts, doPoissonNoise, flNoiseT0):

        geneScores = []

        # init geneScores array
        for i in range(len(Fitness.GENES)):
            geneScores.append([0] * nCycles)

        for cycleIndex in range(nCycles):
            if cycleIndex % 10 == 0:
                print("\t", cycleIndex, ": ", end='', flush=True)
            print('.', end='', flush=True)
            if (cycleIndex + 1) % 10 == 0:
                print('')

            sampleStress = sample
            sampleT0 = Fitness.getRefTime0Sample()

            # 1. bootstrap barcode indeces if needed
            Fitness.updateBARCODE_INDICES()
            if doBootstrapIndeces:
                Fitness.bootstrapBARCODE_INDICES()
                sampleT0 = Fitness.getRefTime0Sample()

            # 2. bootstrap read counts if needed
            if doBootstrapReadCounts:
                sampleStress = Fitness.bootstrapSampleReadCounts(sampleStress)
                if flNoiseT0:
                    sampleT0 = Fitness.bootstrapSampleReadCounts(sampleT0)

            # 3. add Poisson noise if neede
            if doPoissonNoise:
                sampleStress = Fitness.buildPoissonNoisedSample(sampleStress)
                if flNoiseT0:
                    sampleT0 = Fitness.buildPoissonNoisedSample(sampleT0)

            fscores = Fitness.buildFitnessScore(sampleStress, sampleT0)
            gscores = Fitness.buildGeneScores(fscores, scoreType)
            for gIndex, geneScore in enumerate(gscores):
                geneScores[gIndex][cycleIndex] = geneScore

        # sort scores
        for gs in geneScores:
            gs.sort()

        print('buildGrandNoisedGeneScores: Done!')
        return geneScores


# @staticmethod
# def processAllConditions(scoreType):
#     # init bracode indices to avoid previous bootstrap
#     Fitness.updateBARCODE_INDICES()

#     for index in range(0, len(Fitness.CONDITIONS)):
#         ss = Fitness.getSample(index)
#         fs = Fitness.buildFitnessScore(ss, Fitness.getRefTime0Sample())
#         gs = Fitness.buildGeneScores(fs, scoreType)

#         itNumber = Fitness.CONDITIONS_LIST[index * 3]
#         fileName = Fitness.BARSEQ_DIR + \
#             "scores_carbon_v0/nnls/gene_scores." + itNumber + ".txt"
#         Fitness.exportGeneScores(fileName, gs)

#         fileName = Fitness.BARSEQ_DIR + "scores_carbon_v0/region_data." + itNumber + ".txt"
#         Fitness.exportRegionData(fileName, ss, fs)

#         print(index, max(gs),
#               Fitness.CONDITIONS_LIST[index * 3],  Fitness.CONDITIONS_LIST[index * 3 + 2])

    @staticmethod
    def updateGENE_scores(gs, gscoresBI, gscoresPN, sample):
        for index, gene in enumerate(Fitness.GENES):
            gene['score'] = gs[index]
            gene['Mb'] = np.mean(gscoresBI[index])
            gene['Vb'] = np.var(gscoresBI[index])
            gene['Mp'] = np.mean(gscoresPN[index])
            gene['Vp'] = np.var(gscoresPN[index])
            gene['V'] = max(gene['Vb'], gene['Vp'])
            gene['n'] = len(gene['barcodeIndeces'])
            gene['barcodeCounts'] = []
            for bIndex in gene['barcodeIndeces']:
                gene['barcodeCounts'].append(sample[bIndex])

    @staticmethod
    def buildStat():
        # Effective variance estimated from median
        Veffm = 0

        # Effective variance estimated from weighted averaged
        Veffa = 0
        barcodeCount = 0
        geneTotal = 0
        nTotal = 0
        VList = []

        for index, gene in enumerate(Fitness.GENES):
            if gene['n'] > 3 and gene['score'] > 0:
                VList.append(gene['V'])
                nTotal += gene['n']
                Veffa += gene['V'] * gene['n']
                geneTotal += 1
                barcodeCount += np.average(gene['barcodeCounts'])

        # from average
        Veffa /= nTotal
        barcodeCount /= geneTotal

        # from median
        VList.sort()
        Veffm = np.median(VList)

        # calculate tscore and pvalue
        alpha = 0.1
        for index, gene in enumerate(Fitness.GENES):
            if gene['n'] > 0:
                gene['avg_barcodeCount'] = barcodeCount
                gene['Veffa'] = Veffa
                gene['Vm'] = ((gene['n'] - 1) * gene['V'] + Veffa) / gene['n']
                gene['tscore'] = gene['score'] / math.sqrt(gene['Vm'] + alpha)
                gene['pvalue'] = stats.t.sf(gene['tscore'], gene['n'] - 1)
            else:
                gene['Vm'] = 0
                gene['tscore'] = 0
                gene['pvalue'] = 1

    @staticmethod
    def cleanGeneScores():
        for gene in Fitness.GENES:
            gene['score'] = 0
            gene['Mb'] = 0
            gene['Vb'] = 0
            gene['Mp'] = 0
            gene['Vp'] = 0
            gene['V'] = 0
            gene['n'] = 0
            gene['Vm'] = 0
            gene['tscore'] = 0
            gene['pvalue'] = 0
            gene['barcodeCounts'] = []

    # @staticmethod
    # def initData(expId):
    #     Fitness.initConditions(expId)
    #     t0Indeces = Fitness.getTimeZeroIndeces()

    #     # Load data
    #     Fitness.loadBPAG()
    #     Fitness.loadCounts()
    #     Fitness.buildREF_TIME0(t0Indeces)
    #     Fitness.cleanBARCODE_COUNTS(t0Indeces)

    #     Fitness.loadGenes()
    #     Fitness.associateGenesWithBarcodes()
    #     # cleanGENES()
    #     Fitness.buildGENOME_SEGMENTS()

    # @staticmethod
    # def doTimeZeros(expId):
    #     indeces = Fitness.getTimeZeroIndeces()

    #     for i in indeces:
    #         t0StressIndex = i
    #         t0Indeces = []
    #         for j in indeces:
    #             if i != j:
    #                 t0Indeces.append(j)
    #         print(t0StressIndex, t0Indeces)
    #         Fitness.buildREF_TIME0(t0Indeces)
    #         doSample(t0StressIndex, Fitness.GENE_SCORES_DIR +
    #                  '/time0_' + expId + '_')

    # @staticmethod
    # def doAll(expId, doNoise=True):
    #     for index in range(0, max(45, len(Fitness.CONDITIONS))):
    #         doSample(index, Fitness.GENE_SCORES_DIR +
    #                  '/sample_' + expId + '_', doNoise)

    # @staticmethod
    # def doSample(sampleIndex, filePrefix, scoreType, doNoise=True):
    #     nCycles = 100
    #     Fitness.cleanGeneScores()
    #     sample = Fitness.getSample(sampleIndex)

    #     fs = Fitness.buildFitnessScore(sample, Fitness.getRefTime0Sample())
    #     gs = Fitness.buildGeneScores(fs, scoreType)

    #     if doNoise:
    #         gscoresBI = Fitness.buildNoisedGeneScores(nCycles, sample,  doBootstrapIndeces=True,
    #                                                   doBootstrapReadCounts=False, doPoissonNoise=False, flNoiseT0=False)
    #         gscoresPN = Fitness.buildNoisedGeneScores(nCycles, sample,  doBootstrapIndeces=False,
    #                                                   doBootstrapReadCounts=False, doPoissonNoise=True, flNoiseT0=False)
    #     else:
    #         gscoresBI = [[0] * nCycles] * len(gs)
    #         gscoresPN = [[0] * nCycles] * len(gs)

    #     updateGENE_scores(gs, gscoresBI, gscoresPN, sample)
    #     buildStat()

    #     fileName = filePrefix + Fitness.getItNum(sampleIndex) + '.json'
    #     with open(fileName, "w") as f:
    #         json.dump({"genes": Fitness.GENES}, f)

    # def main0():
    #     initConditions('FEBA_133')
    #     t0Indeces = getTimeZeroIndeces()

    #     # Load data
    #     loadBPAG()
    #     loadCounts()
    #     buildREF_TIME0(t0Indeces)
    #     cleanBARCODE_COUNTS(t0Indeces)

    # def main():
    #     print('Lets do it!')
    #     doNoise = False
    #     init()

    #     for expId in CONDITIONS_CONFIG:
    #         print('DOING ' + expId)
    #         initData(expId)
    #         doAll(expId, doNoise)
    #         # doTimeZeros(expId)

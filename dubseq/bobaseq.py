import json
import os
import subprocess
import argparse
import sys
from .core import util

class Context:

    bobaseq_data_dir = None
    bobaseq_pipeline_path = None
    cfg_template_filename = '/home/olga/Development/dubseq/DubSeq/dubseq/cfg/_bobaseq_template_cfg.json'
    cfg_json = None

    @staticmethod
    def build_context(args):
        Context.bobaseq_data_dir = args.bobaseq_data_dir

        with open(Context.cfg_template_filename) as f:
            Context.cfg_json = json.load(f)

        Context.bobaseq_pipeline_path = Context.cfg_json['bobaseq_pipeline_path']

def parse_args():
    parser = argparse.ArgumentParser(
        description='''
        This program executes the Boba-seq pipeline.
        
        Example to run the Boba-seq pipeline:

        python3 -m dubseq.bobaseq -d /path/to/bobaseq_data_dir
        ''',
        formatter_class=util.RawDescriptionArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-d', '--data_dir', 
                        dest='bobaseq_data_dir', 
                        help='input path to bobaseq data directory', 
                        type=str, 
                        required=True)    
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()

def check_args(args):
    print('-----', args.bobaseq_data_dir )
    if(args.bobaseq_data_dir is None):
        print('''ERROR: No bobaseq data direcotry provided.''')
        return False
    
    return True


class Bobaseq:

    def generate_bpag_file(self):

        #get lib_name
        cfg_file_name = os.path.join(Context.bobaseq_data_dir, "out/bobaseq_cfg.json")
        with open(cfg_file_name) as f:
            cfg_json = json.load(f)
        lib_name = cfg_json['lib_names'][0]

        #get the path to the main output file with coordinates and barcodes
        bobaseq_barcodes_filename = os.path.join(Context.bobaseq_data_dir, 
            f"out/bobaseq_pipeline/{lib_name}/05-BC_and_genes_dfs/BC2best_pos.tsv")

        # generate bpag
        bpag_filename = os.path.join(Context.bobaseq_data_dir, "out/bpag.tsv")
        with open(bpag_filename, 'w') as fout:
            fout.write('barcode_up\tbarcode_dn\tbpair_count\tpos_from\tpos_end\tregion_len\trecommended\tup_read_count\tup_contig_id\tup_pos\tup_strand\tdn_read_count\tdn_contig_id\tdn_pos\tdn_strand\n')
            with open(bobaseq_barcodes_filename) as fin:
                # skip header
                fin.readline()

                #process all lines
                for line in fin:
                    line = line.strip()
                    vals = line.split('\t')
                    barcode_up = vals[0]
                    barcode_dn = '----------'
                    bpair_count = vals[7]
                    pos_from = vals[1]
                    pos_end = vals[2]
                    region_len = vals[6]
                    recommended = '----------'
                    up_read_count = vals[7]
                    up_contig_id = vals[3]
                    up_pos = 0
                    up_strand = vals[5]
                    dn_read_count = 0
                    dn_contig_id = '----------'
                    dn_pos = 0
                    dn_strand = '-'

                    fout.write(f'{barcode_up}\t{barcode_dn}\t{bpair_count}\t{pos_from}\t{pos_end}\t{region_len}\t{recommended}\t{up_read_count}\t{up_contig_id}\t{up_pos}\t{up_strand}\t{dn_read_count}\t{dn_contig_id}\t{dn_pos}\t{dn_strand}\n')

    def run_bobaseq_pipeline(self, start_step="1"):
        cfg_file = os.path.join(Context.bobaseq_data_dir, "out", "bobaseq_cfg.json")
        inp_dir = os.path.join(Context.bobaseq_data_dir, "in")
        op_dir = os.path.join(Context.bobaseq_data_dir, "out", "bobaseq_pipeline")
        subprocess.call(["python3", Context.bobaseq_pipeline_path, cfg_file, inp_dir, op_dir, start_step])


    def create_cfg_file(self):

        # retrieve the file names for the cfg file
        lib_names = self.find_files_by_ext('.fq', Context.bobaseq_data_dir, 'in/split_files')
        lib_genome_filenames = self.find_files_by_ext('.fna', Context.bobaseq_data_dir, 'genome')
        lib_genome_gffs = self.find_files_by_ext('.gff', Context.bobaseq_data_dir, 'genome')


        Context.cfg_json['lib_names'] = [lib_names[0].split(".")[0]]
        Context.cfg_json['lib_genome_filenames'] = lib_genome_filenames
        Context.cfg_json['lib_genome_gffs'] = lib_genome_gffs
        Context.cfg_json['lib_genome_dir'] =  os.path.join(Context.bobaseq_data_dir, "genome")

        # save the bobaseq config file
        bobaseq_cfg_filename = os.path.join(Context.bobaseq_data_dir, "out", "bobaseq_cfg.json")
        with open(bobaseq_cfg_filename, 'w') as f:
            json.dump(Context.cfg_json, f, indent=4, sort_keys=False)

    def find_files_by_ext(self, ext, base_dir, sub_dir):
        files = []
        file_dir = os.path.join(base_dir, sub_dir)
        for item in os.listdir(file_dir):
            if item.endswith(ext):
                files.append(item)
        return files
    

def main():
    bs = Bobaseq()
    bs.create_cfg_file()
    bs.run_bobaseq_pipeline()
    bs.generate_bpag_file()

# bs = Bobaseq()
# bs.create_cfg_file("/home/olga/Development/dubseq/DubSeq/data/bobaseq_demux_test/N4/bobaseq")
# bs.run_bobaseq_pipeline("/home/olga/Development/boba-seq/Boba-seq/src/run_steps.py", "/home/olga/Development/dubseq/DubSeq/data/bobaseq_demux_test/N4/bobaseq")
# bs.generate_bpag_file("/home/olga/Development/dubseq/DubSeq/data/bobaseq_demux_test/N4/bobaseq")


if __name__ == "__main__":
    args = parse_args()
    if check_args(args):
        Context.build_context(args)
        main()
import json
import os
import subprocess

#find_cfg_file finds the path to the specific name of file by searching in specific directory
def find_cfg_file(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)    

def find_files_by_ext(ext, base_dir, sub_dir):
    files = []
    file_dir = os.path.join(base_dir, sub_dir)
    for item in os.listdir(file_dir):
        if item.endswith(ext):
            files.append(item)
    return files


        
class Bobaseq:

    def generate_bpag_file(self, bobaseq_data_dir):

        #get the output dir names
        # out_dir_name = find_files_by_ext(".fq", bobaseq_data_dir, "in/split_files")
        # out_dir_name = [out_dir_name[0].split(".")[0]]

        cfg_file_name = os.path.join(bobaseq_data_dir, "out/bobaseq_cfg.json")
        with open(cfg_file_name) as f:
            cfg_json = json.load(f)
        lib_name = cfg_json['lib_names'][0]


        #get the path to the main output file with coordinates and barcodes
        bobaseq_barcodes_filename = os.path.join(bobaseq_data_dir, f"out/bobaseq_pipeline/{lib_name}/05-BC_and_genes_dfs/BC2best_pos.tsv")

        # #read file line by line and save to an array(whitespaces removed)
        # with open(bobaseq_barcodes_filename) as f:
        #     lines = [line.rstrip() for line in f]

        bpag_filename = os.path.join(bobaseq_data_dir, "out/bpag.tsv")




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
                    




    def run_bobaseq_pipeline(self, bobaseq_pipeline_dir, bobaseq_data_dir, start_step="1"):
        cfg_file = os.path.join(bobaseq_data_dir, "out", "bobaseq_cfg.json")
        inp_dir = os.path.join(bobaseq_data_dir, "in")
        op_dir = os.path.join(bobaseq_data_dir, "out", "bobaseq_pipeline")
        subprocess.call(["python3", bobaseq_pipeline_dir, cfg_file, inp_dir, op_dir, start_step])


    def create_cfg_file(self, bobaseq_data_dir):

        # Reads the JSON from the bobaseq template config file
        cfg_template_filename = find_cfg_file('bobaseq_template_cfg.json', './')
        with open(cfg_template_filename) as f:
            cfg_json = json.load(f)

        # Retrieve the file names for the cfg file
        lib_names = find_files_by_ext('.fq', bobaseq_data_dir, 'in/split_files')
        # lib_names = find_files_by_ext('.fq', bobaseq_data_dir, 'in')
        lib_genome_filenames = find_files_by_ext('.fna', bobaseq_data_dir, 'genome')
        lib_genome_gffs = find_files_by_ext('.gff', bobaseq_data_dir, 'genome')

       
        lib_names = [lib_names[0].split(".")[0]]


        cfg_json['lib_names'] = lib_names
        cfg_json['lib_genome_filenames'] = lib_genome_filenames
        cfg_json['lib_genome_gffs'] = lib_genome_gffs
        cfg_json['lib_genome_dir'] =  os.path.join(bobaseq_data_dir, "genome")

        # Save the bobaseq config file
        bobaseq_cfg_filename = os.path.join(bobaseq_data_dir, "out", "bobaseq_cfg.json")
        with open(bobaseq_cfg_filename, 'w') as f: 
            json.dump(cfg_json, f, indent=4, sort_keys=False)

    def run_history(self):
        subprocess.call(["ls", "-l"])

bs = Bobaseq()
bs.create_cfg_file("./../data/bobaseq")
bs.run_bobaseq_pipeline("/home/olga/Development/boba-seq/Boba-seq/src/run_steps.py", "/home/olga/Development/dubseq/DubSeq/data/bobaseq")
bs.generate_bpag_file("./../data/bobaseq")
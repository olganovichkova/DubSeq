import argparse
import sys
from .core import util


class Context:

    gff_fname_in = None
    gff_fname_out = None
    check_gff = False
    optimize_gff = False


    @staticmethod
    def build_context(args):
        Context.gff_fname_in = args.gff_fname_in
        Context.gff_fname_out = args.gff_fname_out
        Context.check_gff = args.check_gff
        Context.optimize_gff = args.optimize_gff


def parse_args():

    parser = argparse.ArgumentParser(
        description='''
        The is a tool to simplify parsing user input commands with Python Arguments Parser.

        Examples to run the bpag program:

        python -m dubseq.gffutil -p /path/to/bpseq.tsv 
               -u  /path/to/bagseq_up.tsv  -d /path/to/bagseq_dn.tsv
               -o /output/dir                                    
        ''',
        formatter_class=util.RawDescriptionArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        dest='gff_fname_in',
                        help='input gff file name for checking or trasformation',
                        type=str,
                        required=True
                        )

    parser.add_argument('-o', '--output',
                        dest='gff_fname_out',
                        help='output filename... ',
                        type=str,
                        required=False
                        )

    parser.add_argument('--check',
                        dest='check_gff',
                        help='add this argument to check if input gff file is correct',
                        action='store_true')

    parser.add_argument('--optimize',
                        dest='optimize_gff',
                        help='add this argument to fix and optimize the gff file',
                        action='store_true')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()

def check_args(args):
    if not args.gff_fname_in:
        print('ERROR: No input file provided.')
        return False
    if args.optimize_gff:
        if args.gff_fname_out:
            return True
        else:
            print('''ERROR: No output file provided. An output file is required if the "--optimize" argument is added.''')
            return False

    if not args.check_gff and not args.optimize_gff:
        print('ERROR: To run, please add either the "--optimize" argument or the "--check" argument.')
        return False
            

def main():
    if Context.check_gff:
        check_gff()
    elif Context.optimize_gff:
        optimize_gff()


def check_gff():
    has_empty_line = False
    has_gene_feature = False

    with open(Context.gff_fname_in, 'r') as f_in:
        for line in f_in:
            line = line.strip()

            if line.startswith('#'):
                continue

            if line == '':
                has_empty_line = True
                continue
            
            vals = line.split('\t')
            if len(vals) < 9:
                raise Exception(f'The number of columns should be at least 9. The problem line: {line}')
                
            f_seqid, f_source, f_type, f_start, f_end, f_score, f_strand, f_phase, f_attributes = vals
            if f_type == 'gene':
                has_gene_feature = True

    
    # print report 
    if has_empty_line:
        print('ERROR: the gff file has empty lines. Please remove them')
    if has_gene_feature == False:
        print('''ERROR: the gff file does not have the required "gene" features. Please use the "--optimize" argument to automatically inject "gene" features in the gff file''')

    if has_empty_line == False and has_gene_feature == True:
        print('The gff file looks valid for processing by dubseq/bobaseq modules')





    

def optimize_gff():
    
    with open(Context.gff_fname_out, 'w') as f_out:    
        with open(Context.gff_fname_in, 'r') as f_in:
            for line in f_in:
                line = line.strip()

                if line == '':
                    continue

                if line.startswith('#'):
                    f_out.write(line)
                    f_out.write('\n')
                    continue

                vals = line.split('\t')
                if len(vals) < 9:
                    raise Exception(f'The number of columns should be at least 9. The problem line: {line}')
                
                f_seqid, f_source, f_type, f_start, f_end, f_score, f_strand, f_phase, f_attributes = vals
                
                # process region type
                if f_type == 'region':
                    f_out.write(line)
                    f_out.write('\n')
                    continue

                # process all other types
                atr_vals = {}
                for attr_val in f_attributes.split(';'):
                    key, val = attr_val.split('=')
                    atr_vals[key] = val


                # generate GENE line
                gene_id = 'gene-' + atr_vals['ID']
                gene_locus_tag = atr_vals['locus_tag'] if 'locus_tag' in atr_vals else 'locus_tag-' + atr_vals['ID']
                gene_line = '\t'.join( [f_seqid, f_source, 'gene', f_start, f_end, f_score, f_strand, f_phase, f'ID={gene_id};locus_tag={gene_locus_tag}'])
                f_out.write(gene_line)                
                f_out.write('\n')

                # mofify and write the original line
                atr_vals['Parent'] = gene_id
                f_attributes = ';'.join( [ f'{key}={val}'  for key, val in atr_vals.items() ] )
                f_line = '\t'.join( [f_seqid, f_source, f_type, f_start, f_end, f_score, f_strand, f_phase, f_attributes])
                f_out.write(f_line)
                f_out.write('\n')

if __name__ == "__main__":
    args = parse_args()
    if check_args(args):
        Context.build_context(args)
        main()
    
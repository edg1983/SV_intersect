'''
SV Intersect script for quick analysis of SV in a trio
Allow annotation with genes, gnomAD AF and additional external samples

Author: Edoardo Giacopuzzi
Email: edoardo.giacopuzzi@well.ox.ac.uk
'''

import re
import sys
import os
from shutil import which,rmtree,move
import argparse
import subprocess

def get_stdout(command, shell=False):
    program = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in iter(program.stdout.readline, b''):
        yield line.decode("utf-8")

def run(command, shell=False):
    program = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = program.communicate()
    return program.returncode, stdout, stderr

def step(name,exitcode,err):
    if exitcode == 0:
        print (name, "completed")
    else:
        error_message="{process} failed with exitcode {code}. Error was:\n{error}".format(process=name,code=exitcode,error=err.decode("utf-8"))
        sys.exit(error_message)

def checkFile(file):
    if not os.path.isfile(file):
        sys.exit("CRITICAL! " + file + " do not exists!")

def getFiles(filelist,mode="dict"):
    if mode == "dict":
        files = {}
        with open(filelist, "r") as f:
            line = f.readline()
            while line:
                line = line.rstrip("\n")
                tokens = line.split('\t')
                files[line[0]] = line[1]
                line = line.rstrip("\n")
    elif mode == "list":
        files = []
        with open(filelist, "r") as f:
            line = f.readline()
            while line:
                line = line.rstrip("\n")
                files.append(line)
                line = line.rstrip("\n")
    return files

def annotate_gnomad(input_file,gnomad_file,output_file,overlap=0.5):
    print("## gnomAD annotation...")
    print ("\tAnnotate", input_file, "with gnomAD AF from", gnomad_file, "- overlap", overlap)
    command = ['bedtools','intersect','-wb','-f',overlap,'-a',input_file,'-b',gnomad_file,'>',output_file]
    exitcode, out, err = run(" ".join(str(c) for c in command), True)
    return exitcode, err

def annotate_external(input_file, bed_files, overlap=0.5):
    print("## External cohort annotation...")
    print("\tAnnotate", input_file, "with AF from", len(bed_files), "additional files")
    print("\toverlap", overlap)
    command = ['bedtools','intersect','-u','-f',overlap,'-a',input_file,'-b']
    command.extend(bed_files)
    external_counts = {}
    for bedline in get_stdout(" ".join(str(c) for c in command), True):
        bedline = bedline.rstrip('\n')
        bedline = bedline.split('\t')
        var_id = "_".join(bedline[0:3])
        if var_id not in external_counts.keys():
            external_counts[var_id] = 1
        else:
            external_counts[var_id] += 1
    return external_counts

def annotate_genes(input_file,genes_file,output_file):
    print("## genes annotation...")
    print ("\tAnnotate", input_file, "with genes definitions from", genes_file)
    command = ['bedtools','intersect','-wo','-a',input_file,'-b',genes_file,'>',output_file]
    exitcode, out, err = run(" ".join(command), True)
    return exitcode, err

def makeBed(input_list, outdir, mode="list"):
    print("### Make bed files ###")
    os.makedirs(outdir, exist_ok=True)
    n=0
    if mode == "list": outfiles = []
    if mode == "dict": outfiles = {}

    with open(input_list) as filelist:
        line = filelist.readline()
        while line:
            n += 1
            line = line.rstrip('\n')

            if mode == "list":
                file = line
                checkFile(file)
                filebase = os.path.basename(file)
                prefix,ext = os.path.splitext(filebase)
                outfile = outdir + '/' + prefix + '.bed'
                outfiles.append(outfile)
            elif mode == "dict":
                tokens = line.split("\t")
                file = tokens[1]
                checkFile(file)
                outfile = outdir + '/' + tokens[0] + '.bed'
                outfiles[tokens[0]] = outfile

            print ("\tProcessing file", n, end="\r")
            command = 'bcftools query -f "%CHROM\\t%POS\\t%INFO/END\\t%FILTER\\t%SVTYPE\\t[%GT]\\n" -i \'GT="alt"\' ' + file + ' | grep PASS | grep -v BND | bedtools sort -i stdin > ' + outfile
            exitcode, out, err = run(command, True)
            if exitcode != 0:
                error_message = "CRITICAL! Error in bed conversion:\n" + err.decode('utr-8')
                sys.exit(error_message)
            line = filelist.readline()

    print("\n\tConversion finished!")
    print("\t", n, "files saved in", outdir)
    return outfiles

def make_anno(var_id,do_gnomad,do_additional,do_genes):
    annotations = []
    if do_gnomad == 1:
        try:
            annotations.append(gnomad_af[var_id]['AC'])
            annotations.append(gnomad_af[var_id]['SVTYPE'])
            annotations.extend(gnomad_af[var_id]['AFs'])
        except:
            annotations.extend(['NONE',0,0,0,0,0,0,0])
    if do_additional == 1:
        try:
            annotations.append(external_af[var_id])
            annotations.append(int(external_af[var_id])/tot_ext_samples)
        except:
            annotations.extend([0,0])
    if do_genes == 1:
        try:
            annotations.append(";".join(genes_anno[var_id]))
        except:
            annotations.append('NO_GENE')
    return annotations

def makeheader(mode,do_gnomad,do_additional,do_genes):
    header = []
    if mode == "denovo":
        header.extend(['CHROM','START','END','FILTER','SVTYPE','GT'])
    if mode == "intersect":
        header.extend(['CHROM','START','END','FILTER','SVTYPE','GT','FATHER_SVTYPE','FATHER_GT','MOTHER_SVTYPE','MOTHER_GT'])

    if do_gnomad == 1:
        header.extend(['gnomAD_tot_AC','gnomAD_SVTYPE','gnomAD_AF','gnomAD_POPMAX_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','gnomAD_EUR_AF'])
    if do_additional == 1:
        header.extend(['EXT_COUNT','EXT_PCT'])
    if do_genes == 1:
        header.append('GENES')
    return header

def analysis(args):
    print("### Trio analysis ###")
    os.makedirs(args.tmpdir, exist_ok=True)

    print("\tTrio files loaded")
    trio_files = makeBed(args.input_trio, args.tmpdir, "dict")
    for key, value in trio_files.items():
        print(key,":",value)

    do_gnomad = 0
    do_additional = 0
    do_genes = 0

    if args.gnomad is not None:
        checkFile(args.gnomad)
        do_gnomad = 1
        gnomad_bed = args.tmpdir + "/tmp_gnomad.tsv"
        (exitcode, err) = annotate_gnomad(trio_files['PROBAND'],args.gnomad,gnomad_bed, args.overlap_others)
        step("\tgnomAD annotation", exitcode, err)
        print("\tFormatting annotations...")
        with open(gnomad_bed) as f:
            line = f.readline()
            while line:
                line = line.rstrip('\n')
                tokens = line.split('\t')
                if args.gnomad_same_svtype == True:
                    if tokens[4] != tokens[9]:
                        line = f.readline()
                        continue

                var_id = '_'.join(tokens[0:3])
                if var_id not in gnomad_af.keys():
                    gnomad_af[var_id] = {}
                    gnomad_af[var_id]['GLOBAL_AF'] = tokens[11]
                    gnomad_af[var_id]['SVTYPE'] = tokens[9]
                    gnomad_af[var_id]['AC'] = tokens[10]
                    gnomad_af[var_id]['AFs'] = tokens[11:]
                else:
                    gnomad_af[var_id]['AC'] += tokens[10]
                    if tokens[11] > gnomad_af[var_id]['GLOBAL_AF']:
                        gnomad_af[var_id]['GLOBAL_AF'] = tokens[11]
                        gnomad_af[var_id]['SVTYPE'] = tokens[9]
                        gnomad_af[var_id]['AFs'] = tokens[11:]
                line = f.readline()

    if args.additional_vcfs is not None or args.additional_beds is not None:
        print("## Additional cohort provided")
        do_additional = 1
        if args.additional_vcfs is not None:
            checkFile(args.additional_vcfs)
            additional_beds = makeBed(args.additional_vcfs, args.tmpdir)
        elif args.additional_beds is not None:
            abbitional_beds = getFiles(args.additional_beds)
            checkFile(args.additional_beds)

        tot_ext_samples = len(additional_beds)
        external_af = annotate_external(trio_files['PROBAND'],additional_beds,args.overlap_others)
        step("\texternal samples annotation", 0, "NONE")

    if args.genes is not None:
        checkFile(args.genes)
        do_genes = 1
        genes_bed = args.tmpdir + "/tmp_genes.tsv"
        (exitcode, err) = annotate_genes(trio_files['PROBAND'],args.genes,genes_bed)
        step("\tGenes annotation", exitcode, err)
        print("\tFormatting annotations...")
        with open(genes_bed, "r") as f:
            line = f.readline()
            while line:
                line = line.rstrip('\n')
                tokens = line.split('\t')
                var_id = '_'.join(tokens[0:3])
                if var_id not in genes_anno.keys(): genes_anno[var_id] = []
                genes_anno[var_id].append('|'.join(tokens[9:]))
                line = f.readline()

    #De novo SV (report non overlapping vars)
    print("## Calculating candidate de novo")
    denovo_bed = args.tmpdir + "/denovo_SV.bed"
    command = ['bedtools','intersect','-v','-f',args.overlap_trio,'-a',trio_files['PROBAND'],'-b',trio_files['MOTHER'],trio_files['FATHER']]
    command.extend(['>',denovo_bed])
    exitcode, out, err = run(" ".join(command), True)
    step("\tdenovo identification", exitcode, err)

    #Intersect with parents (report var and parent overlaps 1 per line, col7:parentID)
    print("## Calculating proband-parents intersect")
    intersect_bed = args.tmpdir + "/intersect_SV.bed"
    command = ['bedtools','intersect','-filenames','-wb','-f',args.overlap_trio,'-a',trio_files['PROBAND'],'-b',trio_files['MOTHER'],trio_files['FATHER']]
    command.extend(['>',intersect_bed])
    exitcode, out, err = run(" ".join(command), True)
    step("\tproband-parents intersection", exitcode, err)

    #Finally, prepare annotated tables: denovo, recessive, intersect
    print("## Performing annotation and creating tables")
    denovo_tab = open(args.outdir + "/denovo_SV.tsv", "w+")
    header = makeheader("denovo",do_gnomad,do_additional,do_genes)
    denovo_tab.write("\t".join(header) + "\n")
    with open(denovo_bed) as f:
        line = f.readline()
        while line:
            line = line.rstrip('\n')
            tokens = line.split('\t')
            var_id = '_'.join(tokens[0:3])
            annotations = make_anno(var_id,do_gnomad,do_additional,do_genes)
            denovo_tab.write(line + "\t" + "\t".join(str(a) for a in annotations) + "\n")
            line=f.readline()
    denovo_tab.close()
    print("\tDenovo table written to: denovo_SV.tsv")

    recessive_tab = open(args.outdir + "/recessive_SV.tsv", "w+")
    intersect_tab = open(args.outdir + "/intersect_SV.tsv", "w+")
    header = makeheader("intersect",do_gnomad,do_additional,do_genes)
    recessive_tab.write("\t".join(header) + "\n")
    intersect_tab.write("\t".join(header) + "\n")
    variant = {}
    with open(intersect_bed) as f:
        line = f.readline()
        while line:
            line = line.rstrip('\n')
            tokens = line.split('\t')
            var_id = '_'.join(tokens[0:3])
            annotations = make_anno(var_id,do_gnomad,do_additional,do_genes)
            parent = tokens[6]
            m = re.search(r'(MOTHER|FATHER)',parent)
            parent = m.group(1)
            variant[var_id] = {}
            variant[var_id]['VAR'] = tokens[0:6]
            variant[var_id]['ANNO'] = annotations
            variant[var_id][parent] = {}
            variant[var_id][parent]['GT'] = tokens[12]
            variant[var_id][parent]['SVTYPE'] = tokens[11]
            line=f.readline()

    for var_id in variant.keys():
        annotation = "\t".join(str(a) for a in variant[var_id]['ANNO'])
        var = "\t".join(str(a) for a in variant[var_id]['VAR'])
        if 'FATHER' not in variant[var_id].keys():
            variant[var_id]['FATHER'] = {}
            variant[var_id]['FATHER']['SVTYPE'] = "NONE"
            variant[var_id]['FATHER']['GT'] = "MISS"
        if 'MOTHER' not in variant[var_id].keys():
            variant[var_id]['MOTHER'] = {}
            variant[var_id]['MOTHER']['SVTYPE'] = "NONE"
            variant[var_id]['MOTHER']['GT'] = "MISS"
        if variant[var_id]['FATHER'] == "0/1" and variant[var_id]['MOTHER'] == "0/1":
            recessive_tab.write("\t".join([var,
                variant[var_id]['FATHER']['SVTYPE'],
                variant[var_id]['FATHER']['GT'],
                variant[var_id]['MOTHER']['SVTYPE'],
                variant[var_id]['MOTHER']['GT'],
                annotation]) + "\n")
        intersect_tab.write("\t".join([var,
            variant[var_id]['FATHER']['SVTYPE'],
            variant[var_id]['FATHER']['GT'],
            variant[var_id]['MOTHER']['SVTYPE'],
            variant[var_id]['MOTHER']['GT'],
            annotation]) + "\n")
    recessive_tab.close()
    intersect_tab.close()

    print("\tRecessive table written to recessive_SV.tsv")
    print("\tIntersection table written to intersect_SV.tsv")

    if args.keep_temp == False:
        print("## Removing temp folder")
        rmtree(args.tmpdir)

    print("All done!")

def makeCohortBed(input_list, outdir, overlap):
    print("### Make cohort bed file ###")
    os.makedirs(outdir, exist_ok=True)
    n=0
    del_file = outdir + '/cases_DEL'
    dup_file = outdir + '/cases_DUP'
    inv_file = outdir + '/cases_INV'
    out_del = open(del_file, "w+")
    out_dup = open(dup_file, "w+")
    out_inv = open(inv_file, "w+")
    with open(input_list) as filelist:
        line = filelist.readline()
        while line:
            n += 1
            line = line.rstrip('\n')

            file = line
            checkFile(file)
            filebase = os.path.basename(file)
            prefix,ext = os.path.splitext(filebase)
            out = {}

            print ("\tProcessing file", n, end="\r")
            command = 'bcftools query -f "%CHROM\\t%POS\\t%INFO/END\\t%FILTER\\t%SVTYPE\\n" -i \'GT="alt"\' ' + file + ' | grep PASS | grep -v BND | bedtools sort -i stdin '
            for bedline in get_stdout(command, True):
                bedline = bedline.rstrip("\n")
                bedline = bedline.split("\t")
                if bedline[4] == "DEL": out_del.write("\t".join(bedline[0:3]) + "\t" + prefix + "\n")
                if bedline[4] == "DUP": out_dup.write("\t".join(bedline[0:3]) + "\t" + prefix + "\n")
                if bedline[4] == "INV": out_inv.write("\t".join(bedline[0:3]) + "\t" + prefix + "\n")
            line = filelist.readline()

    print("\n\tConversion finished!")
    out_del.close()
    out_dup.close()
    out_inv.close()

    print("\tMerging regions with overlap threhsold", overlap)
    outfiles = []
    for file in [del_file,dup_file,inv_file]:
        out_file = file + ".merged"
        outfiles.append(out_file)
        out_merge = open(out_file, "w+")
        command = 'bedtools sort -i {input} | bedmap --count --echo-map-range --echo-map-id-uniq --fraction-both {overlap} --delim "\\t" - | awk \'$1>0\' | sort -u | cut -f2-'.format(overlap=overlap, input=file)
        for bedline in get_stdout(command, True):
                bedline = bedline.rstrip("\n")
                bedline = bedline.split("\t")
                samples = bedline[3].split(';')
                out_merge.write("\t".join(bedline[0:3])+ "\t" + str(len(samples)) + "\t" + bedline[3] + "\n")
        out_merge.close()

    print("\tMerging completed!")
    return n, outfiles

def updateRegions(regions, new_anno, name):
    for var_id in regions.keys():
        if var_id in new_anno.keys():
            regions[var_id][name] = new_anno[var_id]
        else:
            regions[var_id][name] = 0
    return regions

def cohort(args):
    print("### Cohort analysis ###")
    os.makedirs(args.tmpdir, exist_ok=True)

    do_gnomad = 0
    do_genes = 0

    checkFile(args.cases_vcfs)
    checkFile(args.controls_vcfs)
    controls_beds = makeBed(args.controls_vcfs, args.tmpdir)
    (tot_cases, merged_regions) = makeCohortBed(args.cases_vcfs, args.tmpdir, args.overlap_cases)

    print("## Analyzing regions")
    print("\tOverlap threshold:", args.overlap_controls)
    candidate_files = []
    tot_controls = len(controls_beds)
    for region_file in merged_regions:
        m = re.search(r'(DEL|DUP|INV)', region_file)
        svtype = m.group(1)
        print("\t# Processing", svtype)
        print("\tIntersect with controls")
        regions = {}
        with open(region_file) as f:
            line = f.readline()
            while line:
                line = line.rstrip("\n")
                line = line.split("\t")
                var_id = "_".join(line[0:3])
                regions[var_id] = {}
                regions[var_id]['region'] = "\t".join(line)
                line = f.readline()

        if args.controls_same_svtype == True:
            use_controls_beds = []
            for bed in controls_beds:
                svtype_bed = bed + "_" + svtype
                command = ['grep',svtype,bed,'>',svtype_bed]
                exitcode, out, err = run(" ".join(command), True)
                use_controls_beds.append(svtype_bed)
        else:
            use_controls_beds = controls_beds.copy()
        controls_count = annotate_external(region_file, use_controls_beds, args.overlap_controls)
        regions = updateRegions(regions, controls_count, 'controls_c')
        header = ['CHROM','START','STOP','CASE_COUNT','CASES','TOT_CASES','CONTROLS_COUNT','TOT_CONTROLS']
        new_file = region_file + ".candidates"

        if args.gnomad is not None:
            checkFile(args.gnomad)
            do_gnomad = 1
            gnomad_bed = args.tmpdir + "/tmp_gnomad.tsv"
            (exitcode, err) = annotate_gnomad(region_file,args.gnomad,gnomad_bed,args.overlap_gnomad)
            step("\tgnomAD annotation", exitcode, err)
            print("\tFormatting annotations...")
            with open(gnomad_bed) as f:
                line = f.readline()
                while line:
                    line = line.rstrip('\n')
                    tokens = line.split('\t')
                    if args.gnomad_same_svtype == True:
                        if tokens[8] != svtype:
                            line = f.readline()
                            continue

                    var_id = '_'.join(tokens[0:3])
                    if var_id not in gnomad_af.keys():
                        gnomad_af[var_id] = {}
                        gnomad_af[var_id]['GLOBAL_AF'] = tokens[11]
                        gnomad_af[var_id]['SVTYPE'] = tokens[9]
                        gnomad_af[var_id]['AC'] = tokens[10]
                        gnomad_af[var_id]['AFs'] = tokens[11:]
                    else:
                        gnomad_af[var_id]['AC'] += tokens[10]
                        if tokens[11] > gnomad_af[var_id]['GLOBAL_AF']:
                            gnomad_af[var_id]['GLOBAL_AF'] = tokens[11]
                            gnomad_af[var_id]['SVTYPE'] = tokens[9]
                            gnomad_af[var_id]['AFs'] = tokens[11:]
                    line = f.readline()
            regions = updateRegions(regions, gnomad_af, 'gnomad')
            header.extend(['gnomAD_tot_AC','gnomAD_SVTYPE','gnomAD_AF','gnomAD_POPMAX_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','gnomAD_EUR_AF'])
            new_file = new_file + ".gnomad"
            print("gnomad annotation done")

        if args.genes is not None:
            checkFile(args.genes)
            do_genes = 1
            genes_bed = args.tmpdir + "/tmp_genes.tsv"
            (exitcode, err) = annotate_genes(region_file,args.genes,genes_bed)
            step("\tGenes annotation", exitcode, err)
            print("\tFormatting annotations...")
            with open(genes_bed, "r") as f:
                line = f.readline()
                while line:
                    line = line.rstrip('\n')
                    tokens = line.split('\t')
                    var_id = '_'.join(tokens[0:3])
                    if var_id not in genes_anno.keys(): genes_anno[var_id] = []
                    genes_anno[var_id].append('|'.join(tokens[8:]))
                    line = f.readline()
            regions = updateRegions(regions,genes_anno,'genes')
            header.append('GENES')
            new_file = new_file + ".genes"

        filebase = os.path.basename(new_file)
        out_file = args.outdir + '/' + filebase + ".bed"
        print("\tSave annotated file:", out_file)
        out_tab = open(out_file, "w+")
        out_tab.write("\t".join(header) + "\n")
        for var_id in regions.keys():
            outline = [regions[var_id]['region'],tot_cases,regions[var_id]['controls_c'],tot_controls]
            if do_gnomad == 1:
                if regions[var_id]['gnomad'] == 0:
                    outline.extend(['NONE',0,0,0,0,0,0,0])
                else:
                    outline.extend([regions[var_id]['gnomad']['SVTYPE'],regions[var_id]['gnomad']['AC']])
                    outline.extend(regions[var_id]['gnomad']['AFs'])
            if do_genes == 1:
                if regions[var_id]['genes'] == 0:
                    outline.append("NO_GENES")
                else:
                    outline.append(";".join(genes_anno[var_id]))
            out_tab.write("\t".join(str(a) for a in outline) + "\n")
    print("Analysis completed!")

    if args.keep_temp == False:
        print("## Removing temp folder")
        rmtree(args.tmpdir)

ap = argparse.ArgumentParser(description="Analysis of SVs")
#ap.add_argument('mode', help="Running mode", choices=['analysis','makebed'], action="store")
subargs = ap.add_subparsers(help="Arguments for each mode", dest="run_mode", required=True)

sub_analysis = subargs.add_parser('trio', help="analysis mode arguments")
sub_analysis.add_argument('-i','--input_trio',help="List of SV vcf for the pedigree. Proband in first line", action="store", required=True)
sub_analysis.add_argument('-o','--outdir',help="Output folder", action="store", required=True)
sub_analysis.add_argument('-g','--genes',help="Gene definitions table", action="store", required=True)
sub_analysis.add_argument('-k','--gnomad',help="gnomad AF for annotation", action="store", required=False)
sub_analysis.add_argument('-s','--gnomad_same_svtype',help="Annotate gnomAD AF only it the SVTYPE match", action="store_true", required=False)
sub_analysis.add_argument('-f','--overlap_trio',help="Fraction of overlap when comparing trio samples", action="store", default=0.5, required=True)
sub_analysis.add_argument('-r','--overlap_others',help="Fraction of overlap when comparing additional samples", action="store", default=0.5, required=False)
sub_analysis.add_argument('-p','--keep_temp',help="Keep temp files", action="store_true", required=False)
sub_analysis.add_argument('-t','--tmpdir',help="Directory for temporary files", action="store", default="intersect_tmp", required=False)
group = sub_analysis.add_mutually_exclusive_group(required=True)
group.add_argument('-v','--additional_vcfs',help="List of SV vcf files for additional individuals to filter", action="store")
group.add_argument('-b','--additional_beds',help="List of bed files for additional individuals to filter", action="store")

sub_makebed = subargs.add_parser('makebed', help="makebed mode arguments")
sub_makebed.add_argument('-i','--input',help="List of SV vcf to convert", action="store", required=True)
sub_makebed.add_argument('-o','--outdir',help="Folder to save converted bed", action="store", required=True)

sub_cohort = subargs.add_parser('cohort', help="cohort mode arguments")
sub_cohort.add_argument('-c','--cases_vcfs',help="List of SV vcf for cases", action="store", required=True)
sub_cohort.add_argument('-n','--controls_vcfs',help="Folder of SV vcf for controls", action="store", required=True)
sub_cohort.add_argument('-f','--overlap_cases',help="Fraction of overlap to merge regions within cases", action="store", default=0.5, required=False)
sub_cohort.add_argument('-r','--overlap_controls',help="Fraction of overlap for filtering a control region", action="store", default=0.5, required=False)
sub_cohort.add_argument('-a','--overlap_gnomad',help="Fraction of overlap to annotate gnomAD regions", action="store", default=0.5, required=False)
sub_cohort.add_argument('-o','--outdir',help="Folder to save result files", action="store", required=True)
sub_cohort.add_argument('-k','--gnomad',help="gnomad AF for annotation", action="store", required=False)
sub_cohort.add_argument('-s','--gnomad_same_svtype',help="Annotate gnomAD AF only it the SVTYPE match", action="store_true", required=False)
sub_cohort.add_argument('-v','--controls_same_svtype',help="Count overlapping regions from controls only it the SVTYPE match", action="store_true", required=False)
sub_cohort.add_argument('-g','--genes',help="Gene definitions table", action="store", required=False)
sub_cohort.add_argument('-p','--keep_temp',help="Keep temp files", action="store_true", required=False)
sub_cohort.add_argument('-t','--tmpdir',help="Directory for temporary files", action="store_true", default="intersect_tmp", required=False)


args = ap.parse_args()

#check modules and python version
print("Python version", sys.version_info[0], "detected")
if sys.version_info[0] < 3:
    sys.exit("Python v3 required. Unload python and then load module python/3.7")

#global variables
genes_anno = {}
external_af = {}
gnomad_af = {}
tot_ext_samples = 0

if which("bedtools") is None: run('module load bedtools', True)
if which("bcftools") is None: run('module load bcftools', True)

#Make output folder
os.makedirs(args.outdir, exist_ok=True)

if args.run_mode == "makebed":
    checkFile(args.input)
    makeBed(args.input, args.outdir)
if args.run_mode == "trio":
    checkFile(args.input_trio)
    analysis(args)
if args.run_mode == "cohort":
    cohort(args)

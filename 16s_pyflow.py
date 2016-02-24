#!/usr/bin/env python2.7
# coding=utf-8

"""
	16s_pyflow.py

	Date:   16/02/2016
	Usage:  For usage instructions run with option --help
	Author: Madis Rumming <mrumming@cebitec.uni-bielefeld.de>
"""



__author__  = "Madis Rumming <mrumming@cebitec.uni-bielefeld.de>"
__copyright__ = "Copyright 2016, Computational Metagenomics, Faculty of Technology, Bielefeld University"

__version__ = "1.0"
__maintainer__ = "Madis Rumming"
__email__ = "mrumming@cebitec.uni-bielefeld.de"
__status__ = "Production"


import argparse
import os.path
import sys
from string import maketrans

# Extending PYTHONPATH for pyflow
sys.path.append("/vol/cmg/share/virtualenvironments/pyflows/lib/python2.7/site-packages/pyflow/")
from pyflow import WorkflowRunner



# Global paths to the binaries
flash = "/vol/cmg/bin/flash"
fastqc = "/vol/cmg/bin/fastqc"



def parse_arguments():
	parser = argparse.ArgumentParser("Performs analysis of CeBiTec sequenced 16s rRNA sequenced samples incl. read merging, simple QC, OTU clustering via QIIME 1.9 (open reference based), taxonomical assignment and basiv overview plots. It starts with scanning all subdirectories (each subdir is one sample with R001 and R002 forward/backward mate pair reads) from designated input direcotry and takes the directory name in a normalized form as sample identifier. For details of name conversion see description of --qiime-metadata option.")
	parser.add_argument("-i", "--input-dir", dest='input_dir', help="Path to home directory of sample subdirectories with .fastq(.gz). Default: %s" % (os.path.abspath('.')), default=os.path.abspath('.'), required=False, type=str)
	parser.add_argument("-o", "--output-dir", dest='output_dir', help="Path to final output directory to write results to. Default: %s" % (os.path.abspath('.')), default=os.path.abspath('.'), required=False, type=str)
	
	
	parser.add_argument("--flash-path", dest='flash', help="Path to flash binary to use. Default: %s" % (flash), default=flash, required=False, type=str)
	parser.add_argument("--flash-threshold", dest='flash_th', help="Threshold of max errors within overlap region. Default: %f" % (0.01), default=0.01, type=float, required=False)
	parser.add_argument("--flash-min-overlap", dest='flash_min', help="Minimal overlap length. Default: %iBP" % (10), default=10, type=int, required=False)
	parser.add_argument("--flash-max-overlap", dest='flash_max', help="Maximal overlap length. Default: %iBP" % (65), default=65, type=int, required=False)
	
	parser.add_argument("--fastqc-path", dest='fastqc', help="Path to fastqc binary to use. Default: %s" % (fastqc), default=fastqc, required=False, type=str)
	
	parser.add_argument("--qiime-steps", dest='qiime_steps', help="Step 4 realted settings for OTU clustering. Default: %s" % ('suppress_step4'), choices=['with_step4', 'suppress_step4', 'both'], default='suppress_step4', required=False)
	parser.add_argument("--qiime-settings", dest='qiime_settings', help="Path to qiime parameter file (optional).", required=False, default=None, type=str)
	parser.add_argument("--qiime-metadata", dest='qiime_metadata', help="Tab-separated metadata for all samples with QIIME-compatible normalized SampleIDs. 'a-z', 'A-Z' and '.' are allowed characters. Sample names as found in the input directory as subdirectories are normalized through replacing the following characters with a '.': _-+%%<BLANK WHITESPACE>;:,/ A header with speaking category names is required, starting with '#SampleID'. Allowed characters for the header are: 'a-z', 'A-Z', '0-9' and '_'. Example header: '#SampleID\\tCondition\\tMedication'", required=False, default=None, type=str)
	
	
	
	p_group = parser.add_argument_group("Stop after", "Stop computation after the choosen step.")
	group = p_group.add_mutually_exclusive_group(required=True)
	group.add_argument("--readmerging", dest='stop_at', action='store_const', const=0)
	group.add_argument("--qc", dest='stop_at', action='store_const', const=1)
	group.add_argument("--demultiplexing", dest='stop_at', action='store_const', const=2)
	group.add_argument("--clustering", dest='stop_at', action='store_const', const=3)
   	group.add_argument("--assign-taxnominy", dest='stop_at', action='store_const', const=4)
	group.add_argument("--alpha-diversity", dest='stop_at', action='store_const', const=5)
	group.add_argument("--beta-diversity", dest='stop_at', action='store_const', const=6)
	group.add_argument("--complete", help="Perform the whole pipeline including plot generation.", dest='stop_at', action='store_const', const=7)
	
	parser.add_argument("--cores", dest='nCores', help="Amount of CPUs to use for parallel jobs. Default: %i" % (48), default=48, required=False, type=int)
	parser.add_argument("--is-continued", dest='continued', help="Enables continuing an erroneous or paused workflow. MUST use the same dataDirRoot as before.", default=False, required=False, action='store_true')
	parser.add_argument("--is-dry-run", dest='dry_run', help="Check workflow without execution.", default=False, action='store_true', required=False)
	
	
	
	args = parser.parse_args()

	return(args)







class RRNa16sWorkflow(WorkflowRunner):
	
	def __init__(self, stop_at, output_dir, samples, flash, flash_th, flash_min, flash_max, fastqc, qiime_steps, qiime_settings, qiime_metadata, nCores):
		self.stop_at = stop_at
		self.output_dir = output_dir
		self.samples = samples
		self.flash = flash
		self.flash_th = flash_th
		self.flash_min = flash_min
		self.flash_max = flash_max
		self.fastqc = fastqc
		self.qiime_steps = qiime_steps
		self.qiime_settings = qiime_settings
		self.qiime_metadata = qiime_metadata
		self.nCores = nCores
	
	def workflow(self):
		os.mkdir(os.path.join(self.output_dir, "flash"))
		print("KEYS: %s" % (", ".join(self.samples.keys())))
		flash_tasks = []
		for key in self.samples.keys():
			print("TaskID: %s for Inputs %s and %s" % (key, self.samples[key][0], self.samples[key][1]))
			flash_tasks.append(key)
			
			cmd = "%s" % (self.flash)\
				+ " %s" % (self.samples[key][0])\
				+ " %s" % (self.samples[key][1])\
				+ " -m %i" % (self.flash_min)\
				+ " -M %i" % (self.flash_max)\
				+ " -d %s" % (os.path.join(self.output_dir, "flash"))\
				+ " -o %s" % (key)\
				+ " -x %f" % (self.flash_th)\
				+ " 2>&1 | tee %s" % (os.path.join(self.output_dir, "flash", key+".log") )
			
			self.addTask(label=key, command=cmd, nCores=1, memMb=8192)

# Perform QC
		if self.stop_at > 0:
			os.mkdir(os.path.join(self.output_dir, "fastqc"))
			cmd = "%s" % (self.fastqc)\
				+ " -o %s" % (os.path.join(self.output_dir, "fastqc"))\
				+ " %s" % (os.path.join(self.output_dir, "flash", "*.extendedFrags.fastq"))
			self.addTask(label="fastqc", command=cmd, isForceLocal=True, dependencies=flash_tasks)

# Perform demultiplexing

		cmd_qiime_base = "source activate qiime1_9_1 && "
		tr_table = maketrans("_-+% ;:,/",".........")

		if self.stop_at > 1:
			os.mkdir(os.path.join(self.output_dir, "qiime"))
			
			metadata_ = {}
			metadata_header = None
			
			if self.qiime_metadata:
				metadataIn = open(os.path.abspath(self.qiime_metadata))
				header = metadataIn.readline()
				header = header.strip().split('\t')
				metadata_header = header[1:]
				for line in metadataIn:
					line = line.strip().split('\t')
					metadata_[line[0]] = line[1:]
				metadataIn.close()
			
			
			
			
			mapping_file = open(os.path.join(self.output_dir, "qiime", "combined_mapping.txt"), "w")
			if metadata_header:
				mapping_file.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tLocalisation\tUnmodifiedIdentifier\t%s\tDescription\n" % ("\t".join(metadata_header)))
				for key in self.samples.keys():
					mapping_file.write("%s\t\t\t\t%s\t%s\tNoNe\n" % (key.translate(tr_table), key, "\t".join(metadata_[key.translate(tr_table)])))
			else:
				mapping_file.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tLocalisation\tUnmodifiedIdentifier\tDescription\n")
				for key in self.samples.keys():
					mapping_file.write("%s\t\t\t\t%s\tNoNe\n" % (key.translate(tr_table), key))

			mapping_file.close()
			
			alphaparam_file = open(os.path.join(self.output_dir, "qiime", "alpha_parameters.txt"), "w")
			alphaparam_file.write("alpha_diversity:metrics chao1,PD_whole_tree,observed_otus,shannon\n")
			alphaparam_file.close()
			
			keys_ = []
			files_ = []
			for key in self.samples.keys():
				keys_.append(key.translate(tr_table))
				files_.append(os.path.join(self.output_dir, "flash", "%s.extendedFrags.fastq" % (key)))
			
			cmd = cmd_qiime_base\
				+ "split_libraries_fastq.py"\
				+ " -i %s" % (",".join(files_))\
				+ " --sample_id %s" % (",".join(keys_))\
				+ " -o %s" % (os.path.join(self.output_dir, "qiime", "splitted_library"))\
				+ " -q 19"\
				+ " --barcode_type 'not-barcoded'"\
				+ " -m %s" % (os.path.join(self.output_dir, "qiime", "combined_mapping.txt"))
			
			self.addTask(label="split_qiime", command=cmd, dependencies=flash_tasks, env={"PATH": "/vol/qiime-1.9/anaconda-2.1x64/bin:/usr/bin:$PATH"}, isForceLocal=True)
	

# Perform OTU clustering and taxonomical classification

		qiime_tasks = []
		qiime_output_dirs = []

		if self.stop_at > 2:
			
			# 'with_step4', 'suppress_step4', 'both'
			
			nCores = self.nCores
			if self.qiime_steps == 'both':
				nCores = self.nCores/2
			
			cmd = cmd_qiime_base\
				+ "pick_open_reference_otus.py"\
				+ " -i %s" % (os.path.join(self.output_dir, "qiime", "splitted_library", "seqs.fna"))\
				+ " -a"\
				+ " -O %i" % (nCores)\
				+ " -m usearch61"\
				+ " -v"
			
			if self.qiime_settings:
				cmd += " -p %s" % (self.qiime_settings)
			
			
			
			if self.qiime_steps == 'with_step4' or self.qiime_steps == 'both':
				if self.stop_at > 3:
					qiime_tasks.append("openref_qiime_tax_w4")
					qiime_output_dirs.append("open_ref_otus_tax_w4")
					cmd_ = cmd\
						+ " -o %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_tax_w4"))
					self.addTask(label="openref_qiime_tax_w4", command=cmd_, dependencies="split_qiime", env={"PATH": "/vol/qiime-1.9/anaconda-2.1x64/bin:/usr/bin:/bin:/usr(local/bin:$PATH"}, nCores=nCores, memMb=4096)
			
					cmd_add = cmd_qiime_base\
						+ "biom add-metadata"\
						+ " -i %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_tax_w4", "otu_table_mc2_w_tax_no_pynast_failures.biom"))\
						+ " -o %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_tax_w4", "otu_table_mc2_w_tax_no_pynast_failures_w_metadata.biom"))\
						+ " -m %s" % (os.path.join(self.output_dir, "qiime", "combined_mapping.txt"))
					qiime_tasks.append("add_metadata_w4")
					self.addTask(label="add_metadata_w4", command=cmd_add, dependencies="openref_qiime_tax_w4", env={"PATH": "/vol/qiime-1.9/anaconda-2.1x64/bin:/usr/bin:/bin:/usr(local/bin:$PATH"}, isForceLocal=True)
				else:
					# Perform OTU clustering ONLY
					qiime_tasks.append("openref_qiime_w4")
					cmd_ = cmd\
						+ " -o %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_w4"))\
						+ " --suppress_taxonomy_assignment"
					self.addTask(label="openref_qiime_w4", command=cmd_, dependencies="split_qiime", env={"PATH": "/vol/qiime-1.9/anaconda-2.1x64/bin:/usr/bin:/bin:/usr(local/bin:$PATH"}, nCores=nCores, memMb=4096)


			if self.qiime_steps == 'suppress_step4' or self.qiime_steps == 'both':
				if self.stop_at > 3:
					qiime_tasks.append("openref_qiime_tax_no4")
					qiime_output_dirs.append("open_ref_otus_tax_no4")
					cmd_ = cmd\
						+ " -o %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_tax_no4"))\
						+ " --suppress_step4"

					self.addTask(label="openref_qiime_tax_no4", command=cmd_, dependencies="split_qiime", env={"PATH": "/bin:/usr/local/bin:/usr/bin:/vol/qiime-1.9/anaconda-2.1x64/bin:$PATH"}, nCores=nCores, memMb=4096)
			
					cmd_add = cmd_qiime_base\
						+ "biom add-metadata"\
						+ " -i %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_tax_no4", "otu_table_mc2_w_tax_no_pynast_failures.biom"))\
						+ " -o %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_tax_no4", "otu_table_mc2_w_tax_no_pynast_failures_w_metadata.biom"))\
						+ " -m %s" % (os.path.join(self.output_dir, "qiime", "combined_mapping.txt"))
					qiime_tasks.append("add_metadata_no4")
					self.addTask(label="add_metadata_no4", command=cmd_add, dependencies="openref_qiime_tax_no4", env={"PATH": "/vol/qiime-1.9/anaconda-2.1x64/bin:/usr/bin:/bin:/usr(local/bin:$PATH"}, isForceLocal=True)
			
				else:
					# Perform OTU clustering ONLY
					qiime_tasks.append("openref_qiime_no4")
					cmd_ = cmd\
						+ " -o %s" % (os.path.join(self.output_dir, "qiime", "open_ref_otus_no4"))\
						+ " --suppress_step4"\
						+ " --suppress_taxonomy_assignment"
					
					
					self.addTask(label="openref_qiime_no4", command=cmd_, dependencies="split_qiime", env={"PATH": "/bin:/usr/local/bin:/usr/bin:/vol/qiime-1.9/anaconda-2.1x64/bin:$PATH"}, nCores=nCores, memMb=4096)

# Compute alpha diversity measures and plots
		if self.stop_at > 4:
			nCores = self.nCores
			if self.qiime_steps == 'both':
				nCores = self.nCores/2
			
			for entry in qiime_output_dirs:
				cmd = cmd_qiime_base\
					+ "alpha_rarefaction.py"\
					+ " -i %s" % (os.path.join(self.output_dir, "qiime", entry, "otu_table_mc2_w_tax_no_pynast_failures_w_metadata.biom"))\
					+ " -o %s" % (os.path.join(self.output_dir, "qiime", entry, "alpha_diversity"))\
					+ " -t %s" % (os.path.join(self.output_dir, "qiime", entry, "rep_set.tre"))\
					+ " -m %s" % (os.path.join(self.output_dir, "qiime", "combined_mapping.txt"))\
					+ " -a"\
					+ " -O %i" % (nCores)\
					+ " -p %s" % (os.path.join(self.output_dir, "qiime", "alpha_parameters.txt"))
				self.addTask(label="alphadiversity_"+entry, command=cmd, dependencies=qiime_tasks, env={"PATH": "/bin:/usr/local/bin:/usr/bin:/vol/qiime-1.9/anaconda-2.1x64/bin:$PATH"}, memMb=1024, nCores=nCores)

# Compute beta diversity measures and plots
		if self.stop_at > 5:
			nCores = self.nCores
			if self.qiime_steps == 'both':
				nCores = self.nCores/2
			
			for entry in qiime_output_dirs:
				cmd = cmd_qiime_base\
					+ "beta_diversity_through_plots.py"\
					+ " -i %s" % (os.path.join(self.output_dir, "qiime", entry, "otu_table_mc2_w_tax_no_pynast_failures_w_metadata.biom"))\
					+ " -m %s" % (os.path.join(self.output_dir, "qiime", "combined_mapping.txt"))\
					+ " -o %s" % (os.path.join(self.output_dir, "qiime", entry, "beta_diversity"))\
					+ " -t %s" % (os.path.join(self.output_dir, "qiime", entry, "rep_set.tre"))\
					+ " -a"\
					+ " -O %i" % (nCores)
				self.addTask(label="betadiversity_"+entry, command=cmd, dependencies=qiime_tasks, env={"PATH": "/bin:/usr/local/bin:/usr/bin:/vol/qiime-1.9/anaconda-2.1x64/bin:$PATH"}, memMb=1024, nCores=nCores)

# Generate taxonomical plots
		if self.stop_at > 6:
			for entry in qiime_output_dirs:
				cmd = cmd_qiime_base\
					+ "summarize_taxa_through_plots.py"\
					+ " -i %s" % (os.path.join(self.output_dir, "qiime", entry, "otu_table_mc2_w_tax_no_pynast_failures.biom"))\
					+ " -o %s" % (os.path.join(self.output_dir, "qiime", entry, "taxa_plots"))\
					+ " -m %s" % (os.path.join(self.output_dir, "qiime", "combined_mapping.txt"))
				
				if self.qiime_settings:
					cmd += " -p %s" % (self.qiime_settings)
				
				self.addTask(label="taxa_plots_"+entry, command=cmd, dependencies=qiime_tasks, env={"PATH": "/bin:/usr/local/bin:/usr/bin:/vol/qiime-1.9/anaconda-2.1x64/bin:$PATH"}, isForceLocal=True)



def check_args(args):
	if not os.path.exists(args.flash):
		sys.exit("Specified binary of flash under %s does not exist." % (args.flash))
	if not os.path.exists(os.path.abspath(args.input_dir)):
		sys.exit("Specified input directory %s does not exist." % (args.input_dir))
	if not os.path.exists(os.path.abspath(args.output_dir)):
		sys.exit("Specified output directory %s does not exist." % (args.output_dir))
	if args.qiime_settings:
		if not os.path.exists(os.path.abspath(args.qiime_settings)):
			sys.exit("Specified qiime settings file %s does not exist." % (args.qiime_settings))
	if args.qiime_metadata:
		if not os.path.exists(os.path.abspath(args.qiime_metadata)):
			sys.exit("Specified qiime metadata file %s does not exist." % (args.qiime_metadata))
		else:
			metadata = open(os.path.abspath(args.qiime_metadata))
			err_string = []
			header = metadata.readline()
			header = header.strip().split("\t")
			if not len(header) > 1:
				err_string.append("No tab '\t' as separator found in header")
			if not header[0][0] == "#":
				err_string.append("Header must start with a number sign '#'")
			if not header[0] == "#SampleID":
				err_string.append("Header must start with '#SampleID'")
			ind = 1
			for line in metadata:
				line = line.strip().split("\t")
				if not len(line) == len(header):
					err_string.append("Not enough values in line #%i. Must match with column count in header." % (ind))
				ind += 1
			metadata.close()
			if not len(err_string) == 0:
				sys.exit("Errors found in metadata file:\n  %s\n" % ("\n  ".join(err_string)))






def main():
	
	args = parse_arguments()
	check_args(args)

	samples = {}
	for entry in os.listdir(args.input_dir):
		if os.path.isdir(os.path.join(args.input_dir, entry)):
			files = os.listdir(os.path.join(args.input_dir, entry))
			files = filter(lambda s : s.endswith(".fastq.gz"), files)
			if len(files) == 2:
				samples[entry] = [os.path.join(args.input_dir, entry, file_) for file_ in files]

	wflow = RRNa16sWorkflow(args.stop_at, args.output_dir, samples, args.flash, args.flash_th, args.flash_min, args.flash_max, args.fastqc, args.qiime_steps, args.qiime_settings, args.qiime_metadata,  args.nCores)
	retval = wflow.run(mode="sge", nCores=args.nCores, memMb="unlimited", isDryRun=args.dry_run, schedulerArgList=['-l','linh=1'])
	sys.exit(retval)


if __name__ == "__main__":
	main()

import argparse
import cFilter
import sFilter
import bedFilter
import baffler
import svHelpers
import pandas as pd
import detexy
import numpy as np

REFERENCE_FASTA='/hpf/largeprojects/adam/local/reference/homosapiens/ucsc/hs37d5/fasta/hs37d5.fa'
CENTROMERES = '/hpf/largeprojects/adam/projects/ewing_sarcoma/analysis/sv/delly/0.5.9/Mar11_Andrej/code/pype2/centromeres.tab'
BAFFLER_LOG = '/hpf/largeprojects/adam/projects/ewing_sarcoma/analysis/sv/delly/0.5.9/Mar11_Andrej/logs/baffler.tmp.' + str(np.random.randint(0, 100000)) + '.log'
def main(tumor_tab, tumor_bam, normal_tab, normal_bam, config, gatk_path=None, run_cFilter=False, run_sFilter=False, run_bedFilter=False, run_baffler=False, run_detexy=False, output_tab=None):

	df_tumor = pd.DataFrame.from_csv(tumor_tab, sep="\t", index_col=False)
	
	if output_tab is not None:
		BAFFLER_LOG = output_tab.strip('.tab') + '.baffler.log'

	if run_cFilter is True:
		print("Running cFilter...")
		cf = cFilter.cFilter(centromere_list=CENTROMERES)
		cf.configure(config)
		cFilter_results = cf.run_cFilter(
			tumor_tab=df_tumor,
			passed_tumor_df=True,
			tumor_bam=tumor_bam,
			normal_bam=normal_bam,
			gatk_path=gatk_path
		)
	else:
		cFilter_results = df_tumor

	if run_bedFilter is True:
		print("Running bedFilter...")
		bdf = bedFilter.bedFilter()
		bdf.configure(config)
		bedFilter_results = bdf.run_BedFilter(
			tumor_tab=cFilter_results,
			passed_tumor_df=True,
			normal_tab=normal_tab
		)
	else:
		bedFilter_results = cFilter_results

	if run_sFilter is True:
		print("Running sFilter...")
		sf = sFilter.sFilter()
		sf.configure(config)
		sFilter_results = sf.run_sFilter(
			tumor_tab=bedFilter_results,
			passed_tumor_df=True,
			tumor_bam=tumor_bam,
			normal_bam=normal_bam
		)
	else:
		sFilter_results = bedFilter_results

	if run_baffler is True:
		print("Running baffler...")
		bf = baffler.baffler(reference_fasta=REFERENCE_FASTA, log=BAFFLER_LOG)
		svHelpers.configure(bf, config, 'baffler')
		baffler_results = bf.run_baffler(
			tumor_tab=sFilter_results,
			tumor_bam=tumor_bam,
			merge=True
		)
	else:
		baffler_results = sFilter_results

	if run_detexy is True:
		print("Running detexy...")
		dt = detexy.detexy()
		dt.configure(config)
		detexy_results = dt.run_detexy(
			tumor_tab = baffler_results,
			passed_tumor_df=True,
			plot_results=False
		)
	else:
		detexy_results = baffler_results
		
	if output_tab is not None:
		detexy_results.to_csv(output_tab, sep="\t", index=False)	

	return detexy_results

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--tumor_tab', required=True)
	parser.add_argument('--tumor_bam', required=True)
	parser.add_argument('--normal_tab', required=True)
	parser.add_argument('--normal_bam', required=True)
	parser.add_argument('--config', required=True)


	parser.add_argument('--output_tab', required=False, default=None)
	parser.add_argument('--sFilter', required=False, action='store_true', dest='sFilter', default=False)
	parser.add_argument('--cFilter', required=False, action='store_true', dest='cFilter', default=False)
	parser.add_argument('--bedFilter', required=False, action='store_true', dest='bedFilter', default=False)
	parser.add_argument('--baffler', required=False, action='store_true', dest='baffler', default=False)
	parser.add_argument('--detexy', required=False, action='store_true', dest='detexy', default=False)
	parser.add_argument('--gatk_path', required=False, default=None)

	args = parser.parse_args()
	main(
			tumor_tab=args.tumor_tab,
			normal_tab=args.normal_tab,
			normal_bam=args.normal_bam,
			tumor_bam=args.tumor_bam,
			output_tab=args.output_tab,
			run_cFilter=args.cFilter,
			run_sFilter=args.sFilter,
			run_bedFilter=args.bedFilter,
			run_baffler=args.baffler,
			run_detexy=args.detexy,
			config=args.config,
			gatk_path=args.gatk_path
	)

	

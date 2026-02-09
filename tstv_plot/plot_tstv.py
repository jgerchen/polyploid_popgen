import matplotlib.pyplot as plt
import argparse
import gzip
import numpy
import itertools
import sys

parser=argparse.ArgumentParser()
parser.add_argument('--n_sites', default=1000000)
parser.add_argument('--plot_n_windows', default=100)
parser.add_argument('--vcf')

args=parser.parse_args()

plot_n_windows=int(args.plot_n_windows)
n_sites=int(args.n_sites)
vcf_input=args.vcf
line_counter=0

nucs={"A":0,"T":1,"C":2,"G":3}
id_nucs={0:"A",1:"T",2:"C",3:"G"}
tv_ts_dict={("A","G"):0,("G","A"):0,("C","T"):0, ("T","C"):0, ("A","T"):1,("T","A"):1,("C","G"):1,("G","C"):1,("A","C"):1,("C","A"):1,("G","T"):1,("T","G"):1}

print("Making arrays...", file=sys.stderr)
GATK_index_array={"FS":0, "QD":1, "MQ":2, "MQRankSum":3, "ReadPosRankSum":4, "SOR":5, "QUAL":6}
GATK_thresholds={"FS":60, "QD":2, "MQ":40, "MQRankSum":-12.5, "ReadPosRankSum":-8, "SOR":3, "QUAL":"None"}
INFO_array=numpy.empty((n_sites,len(GATK_index_array)))
INFO_array.fill(numpy.nan)
tstv_array=numpy.empty(n_sites)
tstv_array.fill(numpy.nan)
print("Made arrays", file=sys.stderr)

curr_sites_end=n_sites


if vcf_input[-3:]==".gz":
	vcf_file=gzip.open(vcf_input, 'rt')
else:
	vcf_file=open(vcf_input)

for line in vcf_file:
	if line[0]=="#":
		if line[0:6]=="#CHROM":
			sample_list=line.strip().split("\t")[9:]
			n_samples=len(sample_list)		
	else:
		line_cats=line.strip().split("\t")
		assert(len(line_cats)>8),"Variant line sould be at least 8 categories, but line %s has only %s" % (line, len(line_cats))
		ref=line_cats[3]
		alt=line_cats[4]
		QUAL=line_cats[5]

		#only look at biallelic SNPs
		if ref in nucs and alt in nucs:
			#ref_array[line_counter]=nucs[ref]
			#alt_array[line_counter]=nucs[alt]
			tstv_array[line_counter]=tv_ts_dict[(ref, alt)]
			info=line_cats[7]
			info_cats=info.split(";")
			try:
				INFO_array[line_counter,GATK_index_array["QUAL"]]=float(QUAL)
			except:
				pass
			for info_cat in info_cats:
				info_cat_subs=info_cat.split("=")
				if info_cat_subs[0] in GATK_index_array:
					INFO_array[line_counter,GATK_index_array[info_cat_subs[0]]]=float(info_cat_subs[1])
			line_counter+=1
			if line_counter%100000==0:
				print("Read %s variants" % line_counter, file=sys.stderr)
			if line_counter==curr_sites_end:
				print("Read to the end of preallocated array size, which is %s sites. Increasing array size by 1000000 sites. To prevent this from happening and speed up runtime of this script increase the n_sites parameter." % curr_sites_end)
				INFO_array_add=numpy.empty((1000000,len(GATK_index_array)))
				INFO_array_add.fill(numpy.nan)
				INFO_array=numpy.vstack((INFO_array,INFO_array_add))
				tstv_array_add=numpy.empty(1000000)
				tstv_array_add.fill(numpy.nan)
				tstv_array=numpy.concatenate((tstv_array,tstv_array_add))
				curr_sites_end=line_counter+1000000
vcf_file.close()
tstv_Y=tstv_array[0:line_counter]
X_all=INFO_array[0:line_counter,]

#make ts/tv plotting function
def plot_tstv(pred_name,pred_out, tstv_values, plot_n_windows, logY, Ythreshold="None"):
	sort_indexes=numpy.argsort(pred_out)
	tstv_sorted=tstv_values[sort_indexes]
	pred_array_sorted=pred_out[sort_indexes]
	#estimate tstv in windows
	tstv_windows=numpy.array([(len(nwindow)-numpy.sum(nwindow))/numpy.sum(nwindow) for nwindow in numpy.array_split(tstv_sorted, plot_n_windows)])
	#estimate mean predictor per window
	mean_pred=numpy.array([numpy.mean(swindow) for swindow in numpy.array_split(pred_array_sorted, plot_n_windows)])
	x_prop=numpy.arange(plot_n_windows)
	fig, ax1 = plt.subplots()
	ax1.set_xlabel('Proportion data')
	ax1.set_ylabel('ts/tv', color='tab:red')
	ax1.set_ylim(0,max(tstv_windows))
	ax1.plot(x_prop, tstv_windows, color='tab:red')
	ax1.tick_params(axis='y', labelcolor='tab:red')
	
	ax2 = ax1.twinx()
	ax2.set_ylabel(pred_name, color='tab:blue')
	ax2.plot(x_prop, mean_pred, color='tab:blue')
	ax2.tick_params(axis='y', labelcolor='tab:blue')
	if logY:
		ax2.set_yscale('log')
	if Ythreshold!="None":
		ax2.axhline(y=Ythreshold, xmin=0, xmax=1, ls=":", c="g")
	n_total=len(pred_out)
	n_missing=numpy.sum(numpy.isnan(pred_out))
	prop_missing=round(float(n_missing)/n_total, 2)
	plt.title("%s, prop missing: %s (%s out of %s)" % (pred_name, prop_missing, n_missing, n_total))
	fig.tight_layout()
	plt.savefig(pred_name+".pdf")
	plt.close()
	print("Plotted %s.pdf" % pred_name, file=sys.stderr)
#plot ts/tv for GATK stats
for GATK_stat in GATK_index_array:
	plot_tstv(GATK_stat, X_all[:,GATK_index_array[GATK_stat]], tstv_Y, plot_n_windows, False, GATK_thresholds[GATK_stat])
plot_tstv(output_prefix+"QUAL_logY", X_all[:,GATK_index_array["QUAL"]], tstv_Y, plot_n_windows, True)

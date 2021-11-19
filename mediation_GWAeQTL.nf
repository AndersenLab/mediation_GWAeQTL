#! usr/bin/env nextflow


/*
~ ~ ~ > * USER INPUT PARAMETERS 
*/

date = new Date().format( 'yyyyMMdd' )


params.gwa = null
params.gwa_dir = null

params.transcripteQTL = "${workflow.projectDir}/bin/eQTL6545forMed.tsv"
params.transcript_exp = "${workflow.projectDir}/bin/tx5291exp_st207.tsv"
params.fix_names = "fix"


 




/*
~ ~ ~ > * OUTPUT DIRECTORY 
*/

params.out = "${params.gwa_dir}/mediation-${date}"






log.info ""
log.info "------------------------------------------"
log.info "  C. elegans mediation of GWAS pipeline   "
log.info "------------------------------------------"
log.info ""

log.info ""
log.info "Result Directory                        = ${params.out}"
log.info "GWA mapping pipeline                    = ${params.gwa}"
log.info "GWA mapping result directory            = ${params.gwa_dir}"
log.info "eQTL                                    = ${params.transcripteQTL}"
log.info "Input expression data of eQTL calling   = ${params.transcript_exp}"
log.info ""







/*
~ ~ ~ > * expression QTL  
*/

 


File transcripteqtl_all = new File("${params.transcripteQTL}")
transcript_eqtl = transcripteqtl_all.getAbsolutePath()



/*
~ ~ ~ > * transcript expression
*/

Channel
	.fromPath("${params.transcript_exp}")
	.into{ transcript_exp_file; 
			transcript_exp_file2}



 

 




/*
~ ~ ~ > * genome matrix used in mapping of target traits
*/

 

params.genoMatri = "${params.gwa_dir}/Genotype_Matrix/Genotype_Matrix.tsv"

Channel
	.fromPath("${params.genoMatri}")
	.into{ 	med_gm; 
			med_gm2}







if("${params.gwa}" == "cegwas2nf") {



/*
~ ~ ~ > * target trait QTL_peaks 
*/



params.qpeak = "${params.gwa_dir}/Mappings/Data/QTL_peaks.tsv"

Channel
	.fromPath("${params.qpeak}")
	.splitCsv(sep: '\t', skip: 1)
   	.into{peaks; 
   		  mediate_peaks}




} else if("${params.gwa}" == "nemascan"){ 


/*
~ ~ ~ > * target trait QTL_peaks 
*/



params.qpeak = "${params.gwa_dir}/Mapping/Processed/QTL_peaks.tsv"

Channel
	.fromPath("${params.qpeak}")
	.splitCsv(sep: '\t', skip: 1)
	.map { tch,marker,logPvalue,TRAIT,tstart,tpeak,tend,peak_id,h2 -> [TRAIT,tch,tstart,tpeak,tend,logPvalue,peak_id,h2] }
   	.into{peaks; 
   		  mediate_peaks}



} else{
	

	println("ERROR. Please choose --gwa (cegwas2nf or nemascan)")
				exit 1

				
}




Channel.fromPath("${params.gwa_dir}/Phenotypes/pr*")
    .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }
	.set{traits_to_mediate}






/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *      Mediation    * < ~ ~ ~
~ ~ > *          data         * < ~ ~
~ > *                           * < ~
=====================================
*/




mediate_peaks
.combine(traits_to_mediate, by: 0)
.set{medQTL_peaks}






process mediation_data {

 
	executor 'local'

	tag {TRAIT}


	input:

		set val(TRAIT),val(tch),val(tstart),val(tpeak),val(tend),val(logPvalue), val(peak_id), val(h2),file(t_file) from medQTL_peaks

	output:

		set val(TRAIT),val(tch),val(tpeak),val(tstart),val(tend), file("${TRAIT}_scaled_mapping.tsv"),file("${TRAIT}_${tch}_${tpeak}_eqtl.tsv") into QTL_phe, QTL_phe2

 

	"""


    Rscript --vanilla "${workflow.projectDir}/bin/mediaton_input.R" ${TRAIT} ${t_file} ${tch} ${tstart} ${tend} ${tpeak} ${transcript_eqtl}


	"""
}






/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *   Multi-Mediation * < ~ ~ ~
~ ~ > *                       * < ~ ~
~ > *                           * < ~
=====================================
*/





QTL_phe
.spread(med_gm)
.spread(transcript_exp_file)
.into{Gpeak_egene; Gpeak_egene_print}





process multi_mediation {


	cpus 1
	memory '2 GB'

    tag {"${TRAIT}_${tch}_${tpeak}"}

	input:
		set val(TRAIT),val(tch),val(tpeak), val(tstart),val(tend), file(pheno), file(tr_eqtl), file(geno), file(texpression) from Gpeak_egene


	output:
		file("${TRAIT}_${tch}_${tpeak}_medmulti.tsv") optional true into result_multi_mediate

		file("${TRAIT}_${tch}_${tpeak}_elist.tsv") optional true into eQTL_gene


	"""

    Rscript --vanilla "${workflow.projectDir}/bin/multi_mediation.R" ${geno} ${texpression} ${pheno} ${tch} ${tpeak} ${TRAIT} ${tr_eqtl}

	"""
}


 





// =====================================
// ~ > *                           * < ~
// ~ ~ > *                       * < ~ ~
// ~ ~ ~ > *    Mediation      * < ~ ~ ~
// ~ ~ > *                       * < ~ ~
// ~ > *                           * < ~
// =====================================






eQTL_gene
 .splitCsv(sep: '\t')
 .into{eQTL_gene_list ; eQTL_gene_list_print}


 eQTL_gene_list
.combine(QTL_phe2, by: [0,1,2])
.spread(med_gm2)
.spread(transcript_exp_file2)
.into{med_eachtx; med_eachtx_print}

 

process simple_mediation {

 
	cpus 1
	memory '2 GB'


    tag {"${TRAIT}_${gene}"}



	input:
		set val(TRAIT),val(tch),val(tpeak),val(gene), val(tstart),val(tend), file(pheno), file(tr_eqtl), file(geno), file(expression) from med_eachtx


	output:
		file("${TRAIT}_${tch}_${tpeak}_${gene}_med.tsv") into result_mediate

		

	"""

    Rscript --vanilla "${workflow.projectDir}/bin/simple_mediation.R" ${gene} ${geno} ${expression} ${pheno} ${tch} ${tpeak} ${TRAIT} ${tr_eqtl}

	"""
}





/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *    Mediation      * < ~ ~ ~
~ ~ > *      summary          * < ~ ~
~ > *                           * < ~
=====================================
*/







// peaks
//  .map {TRAIT,tch,tstart,tpeak,tend,logPvalue,var_exp,h2 -> TRAIT}.view()
// .unique()
// .set{traits_list}




process summary_mediation {

	cpus 2
	memory '32 GB'

	publishDir "${params.out}/mediation/file_summary", mode: 'copy', pattern: "*mediation.tsv"

	publishDir "${params.out}/mediation/plot_summary", mode: 'copy', pattern: "*plot.png"

	input:
	 set val(TRAIT),val(tch),val(tstart),val(tpeak),val(tend),val(logPvalue), val(peak_id),val(h2) from peaks
	 file(indimed) from result_mediate.collect()
	 file(multimed) from result_multi_mediate.collect()

	output:
	set val(TRAIT), file("${TRAIT}_mediation.tsv")  
	file("*plot.png") optional true


	"""

 
	cat ${TRAIT}_*medmulti.tsv > ${TRAIT}_multi_mediation_analysis.tsv
 
	cat ${TRAIT}_*med.tsv  > ${TRAIT}_indiv_mediation_analysis.tsv
 
	Rscript --vanilla "${workflow.projectDir}/bin/summary_mediation.R" ${TRAIT}_multi_mediation_analysis.tsv ${TRAIT}_indiv_mediation_analysis.tsv ${TRAIT}

 
	"""
}



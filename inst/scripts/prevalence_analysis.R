#!/usr/bin/env Rscript
library(optparse)

opt_list <- list(
	make_option("--asv_table", action = "store_true", type = "character",
				help = "Name of the file with the ASV table"),
	make_option("--idtaxa", action = "store_true", type = "character",
				help = "Name of the file with idtaxa labels"),
	make_option("--outprefix", action = "store_true", type = "character",
    			default = "prevalence",
            	help = "Name of the output file with the prevalence stats, the full file name will
        			be {outdir}/{outprefix}_prevalence.rds"),
	make_option("--outdir", action = "store_true", type = "character",
            	default = tempdir(),
            	help = "Location of the output directory"),
	make_option("--cores", action = "store_true", type = "numeric",
    	        default = 4,
        	    help = "Number of parallel cpus to use"))


opt <- parse_args(OptionParser(option_list = opt_list))


options(mc.cores = opt$cores)

stopifnot(file.exists(opt$asv_table), file.exists(opt$idtaxa))

out_file <- file.path(opt$outdir,"error_rates",paste0(opt$outprefix,"_prevalence.rds"))
summary_file <- str_replace(out_file,".rds","_summary.rds")
plot_file <- file.path(opt$outdir,"error_rates",paste0(opt$outprefix,"_diagnostics_prevelence.pdf"))

library(magrittr)
library(tidyverse)
library(microbiome.onglab)

message("Prevalence analysis:")
message("ASV table: ",opt$asv_table)
message("idtaxa results: ",opt$idtaxa)
message("Output file: ", out_file)

asv_table <- readRDS(opt$asv_table)
idtaxa <- readRDS(opt$idtaxa)

prevalence <- tibble(
	id = str_c("asv",seq_len(ncol(asv_table)),sep = "_"),
	seq = colnames(asv_table),
	prev = apply( asv_table, 2 , FUN = function(x)sum(x  > 0)),
	nreads = colSums(asv_table),
	ave_reads = colMeans(asv_table))
	
idtaxa %<>% as_tibble()
	
	
any_na <- function(x)any(is.na(x))
	
idtaxa %<>%
	mutate( id = str_c("asv",seq_len(ncol(asv_table)), sep = "_")) %>% 
	mutate_if(any_na, list( ~ fct_explicit_na(.,"unlabelled")))
	
nsamples <- nrow(asv_table)

prevalence %<>% inner_join(idtaxa,by = "id")
	


prevalence_summary <- prevalence %>% 
	group_by(phylum) %>% 
	summarize(
		n = n(),
		nreads = sum(nreads),
		max_reads = max(nreads),
		prev = sum(prev)) %>% 
	mutate(
		prev_ratio = prev / n)
		
prevalence_summary %>% arrange(prev_ratio) %>% as.data.frame %>% 
	mutate( cum_reads = cumsum(nreads),
			prop = round(cum_reads / max(cum_reads) * 100,4))
	
prev_check <- prevalence_summary %>% 
	filter(prev_ratio <= 1.2 ) %>% 
	pull(phylum) %>% as.character()
	
prevalence %<>%  
	mutate(
		to_check = if_else(phylum %in% prev_check,"yes","no"))
theme_set(theme_bw())
		
prev_plot <- prevalence %>% 
	ggplot()+geom_point(aes(nreads,prev / nsamples,colour = to_check), alpha = .7)+
	scale_x_log10( labels = scales::comma_format(accuracy = 1))+
	scale_y_continuous( labels = scales::percent_format(accuracy = 1))+
	labs(
		colour = "reads in \n single asv",
		x = "total abundance",
		y = "prevalence")+
	facet_wrap( . ~ phylum)+
	scale_color_manual(values = c(`yes` = "red",`no` = "navyblue"))+
	theme(
		legend.position = "top",
		strip.background = element_blank())
		
prevalence %>% saveRDS(out_file)
prevalence_summary %>% saveRDS(summary_file)

ggsave(
	filename = plot_file,
	plot = prev_plot,
	width = 10,
	height = 9.2,
	units = "in")


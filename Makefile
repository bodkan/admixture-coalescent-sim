SHELL := /bin/bash

output_dir := output
fig_dir := fig
script_dir := src

dirs := $(output_dir) $(fig_dir)

tree_output := $(output_dir)/trees.txt

simulator_bin := scrm

# number of individuals to simulate from each population
n_afr := 50
n_eur := 50
n_nea := 1

# number of loci to sample from each individual and their length
n_haplotypes := 5
hap_length := 1000000

# generation time
gen_time := 25

# effective population size (= the size of present-day European populations)
# used in different scaling equations bellow
N0 := 512000

# convert time in years to units of 4*N_0 generations as used by simulator
scale_time = $(shell echo '($(1) / $(gen_time)) / (4 * $(N0))' | bc -l)
# get a fraction of a given population size to N0
scale_Ne = $(shell echo '$(1) / $(N0)' | bc -l)
# scale the migration parameter m by $ * N_0
scale_m = $(shell echo '4 * $(1) * $(N0)' | bc -l)

# modern human-Neandertals split time (Vernot and Akey, Science 2014)
T_nea_modh_split := 700000
# African/non-African split and bottleneck
T_afr_nonafr_split := 60000

# starts of growth periods of modern human populations
T_gradual_growth := 23000 # Gravel et al., PNAS 2011
T_exp_growth := 5115 # Tennesse et al., Science 2012

# size of the founding Eurasian population (Gravel et al., PNAS 2011)
Ne_nonafr_bottleneck := 1861

# Neandertal population has a constant size in the sumulation
# (Vernot and Akey, Science 2014)
Ne_nea := 1500

# size of an Ancestral african population prior to 148kya
Ne_anc_afr := 7300
T_afr_step_growth := 148000

# population sizes after exponential growth at 5115 years ago until now
# (all from Tennessen et al., Science 2012)
Ne_afr_after_exp_growth := 424000 
Ne_eur_after_exp_growth := $(N0)
afr_exp_growth_rate := 0.0166
eur_exp_growth_rate := 0.0195

# population size of Africans and Europeans at the beginning of exp. growth
# 5115 years ago
Ne_afr_before_exp_growth := $(shell echo '$(Ne_afr_after_exp_growth) / e($(afr_exp_growth_rate) * $(T_exp_growth) / $(gen_time))' | bc -l)
Ne_eur_before_exp_growth := $(shell echo '$(Ne_eur_after_exp_growth) / e($(eur_exp_growth_rate) * $(T_exp_growth) / $(gen_time))' | bc -l)

# scaled exponential growth rates of Africans and Europeans (re-calculated for
# times given as units of 4*N0 generations)
scaled_afr_exp_growth_rate := $(shell echo 'l($(Ne_afr_after_exp_growth) / $(Ne_afr_before_exp_growth)) / $(call scale_time, $(T_exp_growth))' | bc -l)
scaled_eur_exp_growth_rate := $(shell echo 'l($(Ne_eur_after_exp_growth) / $(Ne_eur_before_exp_growth)) / $(call scale_time, $(T_exp_growth))' | bc -l)

# scaled gradual growth rate of Europeans between 23k and 5115 years ago
scaled_eur_gradual_growth_rate := $(shell echo 'l($(Ne_eur_before_exp_growth) / $(Ne_nonafr_bottleneck)) / $(call scale_time, ($(T_gradual_growth) - $(T_exp_growth)))' | bc -l)

# migration rates (Gravel et al., PNAS 2011)
m_afr_euras := 15*10^-5
m_afr_eur := 2.5*10^-5

# introgression parameters (Vernot and Akey, Science 2014)
T_1st_pulse_start := $(shell echo '(36000 + 55000) / 2' | bc)
T_1st_pulse_end := $(shell echo '$(T_1st_pulse_start) + 500' | bc)
m_first_pulse := 0.0015

# mutation rate per site per generation
per_gen_mut_rate := 2.5*10^-8
# cross-over probability between adjacent bases per generation 
per_gen_crossover_rate := 1*10^-8

# population scale parameters for an entire simulated locus
pop_mut_rate := $(shell echo '4 * $(N0) * $(per_gen_mut_rate) * $(hap_length)' | bc -l)
pop_crossover_rate := $(shell echo '4 * $(N0) * $(per_gen_crossover_rate) * ($(hap_length)-1)' | bc -l)

tree_files := $(addsuffix .newick, $(addprefix $(output_dir)/locus_, $(shell seq 1 $(n_haplotypes))))
fasta_files := $(addsuffix .fasta, $(addprefix $(output_dir)/locus_, $(shell seq 1 $(n_haplotypes))))
phylogeny_plots := $(addsuffix .pdf, $(addprefix $(fig_dir)/phylogeny_locus_, $(shell seq 1 $(n_haplotypes))))
coalescent_plots := $(addsuffix .pdf, $(addprefix $(fig_dir)/coalescent_locus_, $(shell seq 1 $(n_haplotypes))))

all: $(phylogeny_plots) $(coalescent_plots)

$(output_dir)/locus_%.fasta: $(output_dir)/locus_%.newick
	seq-gen -mHKY -l $(hap_length) -p `wc -l $< | cut -f1 -d' '` < $< > $<_tmp; \
	tail -n+2 $<_tmp | awk '{ print ">"$$1"\n"$$2 }' > $@; \
	rm $<_tmp

$(output_dir)/locus_%.newick: $(dirs) $(tree_output)
	python3 $(script_dir)/split_tree_output.py \
		--input=$(tree_output) \
		--output_prefix=locus_

$(fig_dir)/phylogeny_locus_%.pdf: $(output_dir)/locus_%.fasta
	Rscript $(script_dir)/plot_phylogeny.R \
		--fasta=$< \
		--n_afr=$(n_afr) \
		--n_eur=$(n_eur) \
		--n_nea=$(n_nea) \
		--output=$@

$(fig_dir)/coalescent_locus_%.pdf: $(output_dir)/locus_%.newick
	Rscript $(script_dir)/find_introgressed_segments.R \
		--tree_file=$< \
		--n_afr=$(n_afr) \
		--n_eur=$(n_eur) \
		--n_nea=$(n_nea) \
		--split_time=$(call scale_time, $(T_nea_modh_split)) \
		--output=$@
	
$(tree_output):
	$(simulator_bin) `expr $(n_afr) + $(n_eur) + $(n_nea)` $(n_haplotypes) -T \
	    `# generate samples from three populations: AFR, EUR, Neandertals` \
	    -I 3 $(n_afr) $(n_eur) $(n_nea) \
\
	    -r $(pop_crossover_rate) $(hap_length) \
\
	    `# constant Neandertal population throughout the whole simulation` \
	    -n 3 $(call scale_Ne, $(Ne_nea)) \
\
	    `# exponential growth in Africans starting 5115 years ago` \
	    -n 1 $(call scale_Ne, $(Ne_afr_after_exp_growth)) \
	    -g 1 $(scaled_afr_exp_growth_rate) \
\
	    `# exponential growth in Europeans starting 5115 years ago` \
	    -g 2 $(scaled_eur_exp_growth_rate) \
\
	    `# constant African population size until 5115 years ago` \
	    -en $(call scale_time, $(T_exp_growth)) 1 $(call scale_Ne, $(Ne_afr_before_exp_growth)) \
\
	    `# gradually growing European population between 23k and 5115 years ago` \
	    -eg $(call scale_time, $(T_exp_growth)) 2 $(scaled_eur_gradual_growth_rate) \
\
  	    `# Neandertal introgression into early modern humans` \
  	    -em $(call scale_time, $(T_1st_pulse_end)) 2 3 $(call scale_m, $(m_first_pulse)) \
  	    -em $(call scale_time, $(T_1st_pulse_start)) 2 3 0 \
\
		`# constant African population size until 148k years ago` \
		-en $(call scale_time, $(T_afr_step_growth)) 1 $(call scale_Ne, $(Ne_anc_afr)) \
\
	    `# split of Africans from non-Africans` \
	    -ej $(call scale_time, $(T_afr_nonafr_split)) 2 1 \
\
	    `# split with Neandertals` \
	    -ej $(call scale_time, $(T_nea_modh_split)) 3 1 \
	    > $@

# for debugging
show_params:
	@echo "Effective population size N0: $(N0)"
	@echo "Mutation rate per bp per generation: $(per_gen_mut_rate)"
	@echo "Cross-over rate per adjacent base pairs per generation: $(per_gen_crossover_rate)"
	@echo
	@echo "Whole locus population mutation rate (theta): $(pop_mut_rate)"
	@echo "Whole locus population cross-over rate (rho): $(pop_crossover_rate)"
	@echo
	@echo "Generation time: $(gen_time)"
	@echo
	@echo "Test of scale_time function: 4 * Ne * 25 years = $(call scale_time, 4*$(N0)*25) in units of 4*Ne"
	@echo 
	@echo "African Ne after exp growth: $(Ne_afr_after_exp_growth) (growth rate: $(afr_exp_growth_rate))"
	@echo "European Ne after exp growth: $(Ne_eur_after_exp_growth) (growth rate: $(eur_exp_growth_rate))"
	@echo
	@echo "African Ne before exp growth: $(Ne_afr_before_exp_growth)"
	@echo "European Ne before exp growth: $(Ne_eur_before_exp_growth)"
	@echo
	@echo "Split between Africans and Eurasians: $(call scale_time, $(T_afr_nonafr_split))"
	@echo "Split between AMH and Neandertals: $(call scale_time, $(T_nea_modh_split))"

$(dirs):
	mkdir $@

clean:
	rm -rf $(dirs)

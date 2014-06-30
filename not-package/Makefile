RHOS = 0 0.025 0.075 0.2 0.5   # values of Rho used in simulation 1
HISTOS = $(patsubst %, histogram_%.pdf, $(RHOS))  # basename of histograms for supplement
SUPPFIGS = $(patsubst %, ../supplements/images/%, $(HISTOS))  # final path for histograms
FIGURES = manuscript/images/rho_and_allele.pdf $(SUPPFIGS) 


all: simulation figures

simulation:  simulations/out_list_01.rda



simulations/out_list_01.rda: simulations/simulation_functions.R simulations/01_simulation_1.R
	cd simulations && Rscript --vanilla 01_simulation_1.R

.PHONY: figures
figures: simulations/02_makecharts_1.R simulations/out_list_01.rda
	cd simulations && \
	Rscript --vanilla 02_makecharts_1.R && \
	mv rho_and_allele.pdf ../manuscript/images/ && \
	mv $(HISTOS) ../supplements/images/



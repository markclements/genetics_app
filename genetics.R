library(tidyverse)
library(stringr)
library(googlesheets)
library(stringi)
library(shiny)
library(DT)
library(shinyjs)


#### helper functions #####

reproduction<-function(P1,P2){## chromosome, gene, alleles, rec_freq ##
	#gametes P1
	P1[[1]]%>%
		group_by(chromosome)%>%
		mutate(alleles=unlist(map2(.x=alleles,.y=rec_freq,~cross_over(.x,.y)))) %>% # recomb
		mutate(hap=sample(c(1,2),size = 1)) %>% 
		mutate(alleles=map2(.x=alleles,.y=hap,~(str_split(.x,",",simplify = T))[.y]))->P1 #ind assort
	
	#gametes P2
	P2[[1]]%>%
		group_by(chromosome)%>%
		mutate(alleles=unlist(map2(.x=alleles,.y=rec_freq,~cross_over(.x,.y)))) %>%
		mutate(hap=sample(c(1,2),size = 1)) %>%
		mutate(alleles=map2(.x=alleles,.y=hap,~(str_split(.x,",",simplify = T))[.y]))->P2
	
	#offspring
	P1$alleles<-map2_chr(.x=P1$alleles,
											 .y=P2$alleles,
											 ~str_c(.x,.y,sep=","))
	P1 %>%
		select(-hap)->P1
	return(P1)
}
### 
cross_over<-function(gen,rf){
	if(runif(1,0,1)>rf) gen<-stri_reverse(gen)
	else gen<-gen
	return(gen)
}


##### SHINY APP ##########	

ui<-fluidPage(tabsetPanel(
	tabPanel("Experiment",
					 fluidRow(
					 	column(
					 		wellPanel(
					 			fluidRow(p(),
					 							 p(),
					 							 h5("Selected Parents"),
					 							 tableOutput("P1")),
					 			fluidRow(
					 				h5("Controls"),
					 				column(
					 					p(),
					 					p(),
					 					actionButton("reset", "Clear Selection"),
					 					width = 3,
					 					offset = 0
					 				),
					 				column(
					 					p(),
					 					p(),
					 					uiOutput("mate"),
					 					width = 1,
					 					offset = 3
					 				)
					 			),
					 			fluidRow(
					 				p(),
					 				p(),
					 				numericInput(
					 					"num_matings",
					 					label = "Expected Number of Offspring",
					 					value = 10,
					 					min = 1,
					 					max = 100
					 				)
					 			)
					 		),
					 		width = 4,
					 		offset = 0.5
					 	),
					 	
					 	column(
					 		wellPanel(
					 			p(),
					 			h5("Available for Mating"),
					 			p(),
					 			DT::dataTableOutput("data")
					 		),
					 		width = 8,
					 		offset = 0.5
					 	)
					 )),
	####tabPanel Experiment####
	
	tabPanel("Results",
					 fluidRow(
					 	column(
					 		wellPanel(
					 			selectizeInput(
					 				"select_grp",
					 				'Select Traits to Summarise',
					 				choices = c("NA"),
					 				multiple = TRUE
					 			),
					 			
					 			numericInput(
					 				"select_gen",
					 				"Select Generation to Summarise",
					 				value = 1,
					 				min = 1,
					 				max = 1
					 			),
					 			#),
					 			#wellPanel(
					 			h5("Parents of this Generation"),
					 			tableOutput("parents")
					 		),
					 		width = 4,
					 		offset = 0.5
					 	),
					 	
					 	
					 	column(
					 		wellPanel(h5("Offspring Summary Table"),
					 							dataTableOutput("summary")),
					 		width = 8,
					 		offset = 0.5
					 	)
					 )) ####tabPanel Results####
))	
			 
	
	
server <- function(input, output, session){
	
####### Begin Initialize ########
	##### generate all possible genotypes ########
	read_csv("files/genes_alleles.csv")->ga ## save temporarily
	
	ga$alleles %>% ## select alleles column from user input sheet
		map(.,~str_split(.,",")) %>%
		flatten() %>%
		map(.,~cross(list(.,.))) %>%
		cross() %>%
		modify_depth(.,2,~str_c(.,collapse = ",")) %>%
		map(.,~unlist((.))) %>%
		map(.,~str_c(.,collapse = "-")) %>%
		tibble(genotype=unlist(.))%>%
		select(genotype)->table ## save as table
	
	
  #####assign allele (and eventually gene) interactions #####
	read_csv("files/genotype_phenotype.csv")->gb ## get other google sheet with user input
	
	gb %>%
		gather(gene,phenotype,2:5) %>% ## clean up, make tidy
		filter(phenotype != "NA")->pheno ## save as phenotype 

  ####generate look-up table for genotypes and phenotypes####
	for (i in 1:dim(pheno)[1]){
		table[str_detect(table$genotype,
										 str_c(gb$genotype[i],
										 			stri_reverse(gb$genotype[i]),sep="|")),
					as.character(pheno[i,2])]<-as.character(pheno[i,3])
	}
#### End Initialize #####

#### reactive values	####
	rv<-reactiveValues()

##### starting data table ####
	### save tibble as reactive value, so it can be modified/added to
	ga %>% ###
		select(-alleles) %>%
		gather(id,alleles,-chromosome,-gene,-rec_freq) %>%
		group_by(id)%>%
		nest() %>%
		mutate(gen=1) %>% 
		mutate(parent1=NA,parent2=NA)%>% ## below is code to assign phenotype to genotype
		mutate(genotype=map(.$data,~map(.$alleles,~(str_c(.))))) %>% #extract genotype from nested tibble
		mutate(genotype=map(.$genotype,~stri_join_list(.,collapse="-"))) %>% #format genotype chr string from nested tibble
		mutate(genotype=unlist(genotype))%>% #format genotype 
		left_join(.,table,by="genotype")->rv$data #join to lookup table = assign phenotype to genotype and save
	
#### Data table that stores all data and allows selection of parents####
	### Using DT for now, may change later 
	### access selected rows via input$data_rows_selected 
	### are integers if server = F; character list if server =T
	output$data<-DT::renderDataTable({rv$data %>%
																		select(-data,-genotype,-parent1,-parent2)
																		},server = F,
																		options = list(searchHighlight = TRUE),
																		filter = 'bottom',
																	 rownames= FALSE)
	
#### enable mate button when two selected ####
	output$mate<-renderUI({
		s=input$data_rows_selected
		if (length(s)<2) return()
		rv$data %>%
			slice(c(s[1],s[2])) %>%
			select(sex) %>%
			n_distinct()->ss
		
		if (ss==2) actionButton("mate","Mate")
		})
	
#### clear selections #### 
	proxy = dataTableProxy('data')
	
	observeEvent(input$reset,{
	proxy %>%
			selectRows(NULL)
	}) 
		
#### reactive to select parents, reproduce, assign phenotypes ####
	observeEvent(input$mate,{
		
		rv$gen<-max(rv$data$gen)
		P1 = input$data_rows_selected[1]
	  P2 = input$data_rows_selected[2]
	  #### IDEAS: 
	  #### save P1 and P2 in rv$data as new column variable (type?)
	  #### parents = P, offspring = O 
	rv$data %>%	
		slice(c(P1,P2)) %>% 
		select(id,data) %>%
		mutate(x=1) %>% #### not sure why this is needed, but it works ####
		spread(id,data) %>%
		summarise(data=map2(.x=.[,2],.y=.[,3],~rerun(input$num_matings,reproduction(.x,.y)))) %>%
		map_df(~.)%>%
		setNames(.,"data") %>%
		mutate(id=as.character(1:n()),gen=rv$gen+1) %>% 
		mutate(parent1=P1,parent2=P2)%>%## assign phenotypes (below)
		mutate(genotype=map(.$data,~map(.$alleles,~(str_c(.))))) %>% #extract genotype from nested tibble
		mutate(genotype=map(.$genotype,~stri_join_list(.,collapse="-"))) %>% #format genotype chr string from nested tibble
		mutate(genotype=unlist(genotype))%>% #format genotype
		left_join(.,table,by="genotype") %>%
		bind_rows(rv$data,.) %>%
		ungroup()->rv$data # join to original data table !!
})
	
##### selected parent UI table #####	
	output$P1<-renderTable({
				P1 = input$data_rows_selected[1]
				P2 = input$data_rows_selected[2]
				if (is.null(P1)) return()
		# Generation 1
		# Sex: Male
		# Body Color: Red
		# etc...

		rv$data %>%
			slice(P1) %>%
			select(-data,-genotype,-parent1,-parent2) %>%
			gather(key,P1)->x
		
		rv$data %>%
			slice(P2) %>%
			select(-data,-genotype,-parent1,-parent2) %>%
			gather(key,P2) %>%
			left_join(x,.,by="key")->y
		
		return(y)
			
	},spacing = 'xs',colnames = T,na = "_") 

	# ### parent 2 UI table ###	
	# output$P2<-renderTable({
	# 	P2 = input$data_rows_selected[2]
	# 	if (is.null(P2)) return(NULL)
	# 	# Generation 1
	# 	# Sex: Male
	# 	# Body Color: Red
	# 	# etc...
	# 	
	# 	rv$data %>%
	# 		slice(P2) %>%
	# 		select(-data,-genotype) %>%
	# 		gather()
	# 	
	# },spacing = 'xs',colnames = F) 
	
#### Reactive results summary table ####
	output$summary<-DT::renderDataTable({
		s<-input$select_gen
		if (is.null(s)) return()
		rv$data %>%
			filter(gen==s) %>%
			select(-data,-genotype,-id,-gen) %>%
			group_by_(.dots=input$select_grp) %>%
			summarise(N=n())
		},selection = 'none',
		rownames=FALSE)
	
##### Reactive results parent mated table #####	
	output$parents<-renderTable({
		s<-input$select_gen
		if (is.null(s)) return()
		# Generation 1
		# Sex: Male
		# Body Color: Red
		# etc...
		rv$data %>%
			filter(gen==s) %>%
			summarise(P1=first(parent1))%>%
			as.numeric()->P1

		rv$data %>%
			filter(gen==s) %>%
			summarise(P2=first(parent2))%>%
			as.numeric()->P2
		
		rv$data %>%
			slice(P1) %>%
			select(-data,-genotype,-parent1,-parent2) %>%
			gather(key,P1)->x
		
		rv$data %>%
			slice(P2) %>%
			select(-data,-genotype,-parent1,-parent2) %>%
			gather(key,P2) %>%
			left_join(x,.,by="key")->y
		
		return(y)
		
	},spacing = 'xs',colnames = T,na = "_")

#### update the selectize, trait summary choices ####
observe({
	ga$gene %>%
		as.list() %>%
		updateSelectizeInput(session,"select_grp", choices = .,selected = .)
})

	
observe({
	rv$data %>%
		select(gen) %>%
		summarise(max=max(gen))%>%
		as.numeric()%>%
		updateNumericInput(session,inputId = "select_gen",max = .,value=.)
})	
			
}
#### RUN ####
shinyApp(ui = ui, server = server)


tibble(sex=c("M","M"))%>%
	n_distinct()

dd %>%
	gather() %>%
	filter() %>%
	knitr::kable()


### generate all possible genotypes
gs_key("1CxFLwMxoW06J7vxWgMFiYTQ6hc6GheLETI8G0aYtPFk") %>%
	gs_read(ws = "genes_alleles") %>%
	write_csv(.,"files/genes_alleles.csv")

ga<-read_csv("files/genes_alleles.csv")
ga$alleles %>%
	map(.,~str_split(.,",")) %>%
	flatten() %>%
	map(.,~cross(list(.,.))) %>%
	cross() %>%
	modify_depth(.,2,~str_c(.,collapse = ",")) %>%
	map(.,~unlist((.))) %>%
	map(.,~str_c(.,collapse = "-")) %>%
	tibble(genotype=unlist(.))%>%
	select(genotype)->table 


#assign allele and gene interactions
gs_key("1CxFLwMxoW06J7vxWgMFiYTQ6hc6GheLETI8G0aYtPFk") %>%
	gs_read(ws = "genotype_phenotype") %>%
	write_csv("files/genotype_phenotype.csv")
gb<-read_csv("files/genotype_phenotype.csv")

gb %>%
	gather(gene,phenotype,2:5) %>%
	filter(phenotype != "NA")->pheno


##generate look-up table for genotypes and phenotypes
for (i in 1:dim(pheno)[1]){
	table[str_detect(table$genotype,
									 str_c(gb$genotype[i],
									 			stri_reverse(gb$genotype[i]),sep="|")),
				as.character(pheno[i,2])]<-as.character(pheno[i,3])
}

## starting data table
ga %>%
	select(-alleles) %>%
	gather(id,alleles,-chromosome,-gene,-rec_freq) %>%
	group_by(id)%>%
	nest() %>%
	mutate(gen=1) %>% ## below is code to assign phenotype to genotype
	mutate(genotype=map(.$data,~map(.$alleles,~(str_c(.))))) %>% #extract genotype from nested tibble
	mutate(genotype=map(.$genotype,~stri_join_list(.,collapse="-"))) %>% #format genotype chr string from nested tibble
	mutate(genotype=unlist(genotype))%>% #format genotype 
	left_join(.,table,by="genotype")->data #join to lookup table = assign phenotype to genotype and save

# select parents, reproduce, assign phenotypes
gen<- max(data$gen)+1
data %>%	
	sample_n(size=2,replace = F) %>% ## select at random for now
	select(id,gen,data) %>%
	spread(id,data) %>%
	summarise(data=map2(.x=.[,2],.y=.[,3],~rerun(100,reproduction(.x,.y)))) %>%
	map_df(~.)%>%
	setNames(.,"data") %>%
	mutate(id=as.character(1:n()),gen=gen) %>% ## assign phenotypes (below)
	mutate(genotype=map(.$data,~map(.$alleles,~(str_c(.))))) %>% #extract genotype from nested tibble
	mutate(genotype=map(.$genotype,~stri_join_list(.,collapse="-"))) %>% #format genotype chr string from nested tibble
	mutate(genotype=unlist(genotype))%>% #format genotype
	left_join(.,table,by="genotype") %>%
	bind_rows(data)->data1 # join to original data table !!

data1 %>%
	select(-data,-id,-gen,-genotype)%>%
	group_by(sex,blood_type,eye_color,body_color) %>%
	summarise(N=n())



summarise(N=n()) %>%
	spread(trait,phenotype)


ggplot(aes(x=phenotype,y=N))+
	geom_bar(stat="identity") +
	facet_grid(sex~trait)
	

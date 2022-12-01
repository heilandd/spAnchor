pb <- progress_estimated(length(spots))


memory <- list()

data.new <-map_dfr(.x=1:20,
                   .f=function(i){
                     memory[[i]] <- system2("memory_pressure", stdout = TRUE)
                     bc_run <- spots[i]
                     data <-
                       scDF %>%
                       dplyr::filter(barcodes==bc_run) %>%
                       dplyr::mutate(celltypes=as.character(celltypes)) %>%
                       dplyr::arrange(celltypes)
                     nr_cells <- nrow(data)
                     n_select <- data %>% count(celltypes)
                     #Initiate population
                     pop <- initiate_Population(nr_of_random_spots, n_select,nested_ref_meta,cell_type_var,mat.ref)
                     qc <- lapply(1:iter_GA, function(zz) run(zz))
                     gc()
                     validate_randoms_select <-
                       map(.x=1:nr_of_random_spots, function(j){fitness(pop[j,], nr_cells, bc_run)}) %>%
                       unlist() %>%
                       as.data.frame() %>%
                       rownames_to_column("order") %>%
                       rename("cor":=.) %>%
                       arrange(desc(cor))
                     pop_select <- pop[as.numeric(validate_randoms_select$order[1]), ]
                     select_cells <- names(which(pop_select==1))
                     data$best_match <- select_cells
                     return(data)
                   })

report.memory <- function(file="report.mem.csv"){
  if(!file.exists(file)){
    mem <- data.frame(x=mem_used() %>% as.numeric()*10e-10 )
    write.csv(mem, file, row.names = F)
  }else{
    mem <- read.csv(file)
    mem <- rbind(mem, data.frame(x=mem_used() %>% as.numeric()*10e-10))
    write.csv(mem, file, row.names = F)
  }
  }
run <- function(zz){
  validate_randoms_select <- lapply(1:nr_of_random_spots, function(j) fitness(pop[j,], nr_cells,bc_run)) %>% unlist()
  names(validate_randoms_select) <- 1:nr_of_random_spots
  validate_randoms_select <- validate_randoms_select[order(-validate_randoms_select)]

  #select parents
  parents <- pop[as.numeric(names(validate_randoms_select)[1:2]), ]

  #Create Children
  offspring_2 <- cross_over(parents,cross_over_point)
  offspring <- mutation_GA(offspring_2,
                           nr_mut,
                           nr_offsprings,
                           ref_meta=ref_meta,
                           cell_type_var = cell_type_var)

  # remove old parents
  remove <- as.numeric(tail(validate_randoms_select, dim(offspring)[1]) %>% names())
  pop.new <- rbind(pop[-remove, ], offspring)

  #update pop
  pop <<- pop.new

  return(mean(validate_randoms_select[1:2]))

}


gc()
base::options(future.fork.enable = TRUE)
future::plan("multisession", workers = 32)
future::supportsMulticore()
base::options(future.globals.maxSize = 20* 10* 1024^2)
message("... Run multicore ... ")

system.time(data.new <-furrr::future_map(.x=1:20,
                   .f=function(i){
                     bc_run <- spots[i]
                     data <-
                       scDF %>%
                       dplyr::filter(barcodes==bc_run) %>%
                       dplyr::mutate(celltypes=as.character(celltypes)) %>%
                       dplyr::arrange(celltypes)
                     nr_cells <- nrow(data)
                     n_select <- data %>% count(celltypes)
                     #Initiate population
                     pop <- initiate_Population(nr_of_random_spots, n_select,nested_ref_meta,cell_type_var,mat.ref)
                     qc <- lapply(1:iter_GA, function(zz) run(zz))
                     gc()
                     validate_randoms_select <-
                       map(.x=1:nr_of_random_spots, function(j){fitness(pop[j,], nr_cells, bc_run)}) %>%
                       unlist() %>%
                       as.data.frame() %>%
                       rownames_to_column("order") %>%
                       rename("cor":=.) %>%
                       arrange(desc(cor))
                     pop_select <- pop[as.numeric(validate_randoms_select$order[1]), ]
                     select_cells <- names(which(pop_select==1))
                     data$best_match <- select_cells
                     return(data)
                   },
                   .options = furrr::furrr_options(seed = TRUE), .progress=T))





gc()
base::options(future.fork.enable = TRUE)
future::plan("multicore", workers = 16)
future::supportsMulticore()
base::options(future.globals.maxSize = 8* 10* 1024^2)
message("... Run multicore ... ")
memory <- seq(5,20, by=5)
system.time(data.new <-furrr::future_map_dfr(.x=1:20,
                             .f=function(i){
                               if(any(i==memory)){gc()}
                               bc_run <- spots[i]
                               data <-
                                 scDF %>%
                                 dplyr::filter(barcodes==bc_run) %>%
                                 dplyr::mutate(celltypes=as.character(celltypes)) %>%
                                 dplyr::arrange(celltypes)
                               nr_cells <- nrow(data)
                               n_select <- data %>% count(celltypes)
                               #Initiate population
                               pop <- initiate_Population(nr_of_random_spots, n_select,nested_ref_meta,cell_type_var,mat.ref)
                               qc <- lapply(1:iter_GA, function(zz) run(zz))
                               gc()
                               validate_randoms_select <-
                                 map(.x=1:nr_of_random_spots, function(j){fitness(pop[j,], nr_cells, bc_run)}) %>%
                                 unlist() %>%
                                 as.data.frame() %>%
                                 rownames_to_column("order") %>%
                                 rename("cor":=.) %>%
                                 arrange(desc(cor))
                               pop_select <- pop[as.numeric(validate_randoms_select$order[1]), ]
                               select_cells <- names(which(pop_select==1))
                               data$best_match <- select_cells
                               return(data)
                             },
                             .options = furrr::furrr_options(seed = TRUE),.progress=T))




gc()
base::options(future.fork.enable = TRUE)
future::plan("sequential", workers = 16)
future::supportsMulticore()
base::options(future.globals.maxSize = 30* 10* 1024^2)
message("... Run multicore ... ")
memory2 <- list()
system.time(data.new <-furrr::future_map(.x=1:20,
                             .f=function(i){
                               memory2[[i]] <<- system2("memory_pressure", stdout = TRUE)
                               bc_run <- spots[i]
                               data <-
                                 scDF %>%
                                 dplyr::filter(barcodes==bc_run) %>%
                                 dplyr::mutate(celltypes=as.character(celltypes)) %>%
                                 dplyr::arrange(celltypes)
                               nr_cells <- nrow(data)
                               n_select <- data %>% count(celltypes)
                               #Initiate population
                               pop <- initiate_Population(nr_of_random_spots, n_select,nested_ref_meta,cell_type_var,mat.ref)
                               qc <- lapply(1:iter_GA, function(zz) run(zz))
                               gc()
                               validate_randoms_select <-
                                 map(.x=1:nr_of_random_spots, function(j){fitness(pop[j,], nr_cells, bc_run)}) %>%
                                 unlist() %>%
                                 as.data.frame() %>%
                                 rownames_to_column("order") %>%
                                 rename("cor":=.) %>%
                                 arrange(desc(cor))
                               pop_select <- pop[as.numeric(validate_randoms_select$order[1]), ]
                               select_cells <- names(which(pop_select==1))
                               data$best_match <- select_cells
                               return(data)
                             },
                             .options = furrr::furrr_options(seed = TRUE),.progress=T))




mem.rep <- read.csv("report.mem2.csv")
plot(mem.rep$x)




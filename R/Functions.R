###############################################################
###########   Deconvolution methods
###############################################################

#' @title  runRCTD cell-type deconvolution
#' @author Dieter Henrik Heiland
#' @description Run the a cell type deconvolution using the SpaceXR/RCTD algorithm from SPATA objects
#' @inherit
#' @param ref The reference seurat object
#' @param cell_type_var The variale in the Seurat meta.data file containing the cell type group
#' @param object The variale in the Seurat meta.data file containing the cell type group
#' @return
#' @examples
#' @export
#'
runRCTD <- function(object,ref, cell_type_var, overwrite=T, return.RCTD=F){

  #Some check up
  if(!any(cell_type_var %in% names(ref@meta.data))) stop(paste0("The variable: ",cell_type_var, " was not found in the seurat object"))
  SPATA2::check_object(object)
  cluster <- cell_type_var
  #Load some packages
  library(spacexr)
  library(Matrix)

  #Get reference data
  counts <- ref@assays$RNA@counts %>% as.matrix()
  meta_data <-
    ref@meta.data %>%
    as.data.frame() %>%
    dplyr::select({{cluster}}, nFeature_RNA)


  #Creat RCTD object
  cell_types <- meta_data[, cluster]
  names(cell_types) <- rownames(meta_data)
  cell_types <- as.factor(cell_types)
  nUMI <- colSums(counts)
  reference <- Reference(counts, cell_types, nUMI)
  counts <- SPATA2::getCountMatrix(object) %>% as.matrix()
  coords <- SPATA2::getCoordsDf(object) %>% dplyr::select(barcodes,x, y) %>% as.data.frame()
  rownames(coords) <- coords$barcodes
  coords <- coords[, c(2:3)]
  names(coords) <- c("xcoord", "ycoord")
  nUMI <- colSums(counts)

  #Run Deconv
  puck <- spacexr::SpatialRNA(coords, counts, nUMI)
  myRCTD <- spacexr::create.RCTD(puck, reference, max_cores = 5)
  myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = "doublet")

  if(return.RCTD==T){
    object <- myRCTD
  }else{


    # Get Results
    results <- myRCTD@results
    norm_weights = spacexr::normalize_weights(results$weights)
    cell_type_names <- myRCTD@cell_type_info$info[[2]]
    spatialRNA <- myRCTD@spatialRNA

    #Create output

    norm_weights <-
      norm_weights %>%
      as.data.frame() %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(object %>%
                         SPATA2::getCoordsDf() %>%
                         dplyr::select(barcodes),., by = "barcodes")

    norm_weights[is.na(norm_weights)] <- 0

    object <- object %>% SPATA2::addFeatures(., norm_weights, overwrite = overwrite)


  }



  return(object)
}


###############################################################
###########   Cell Type Deconv from scores
###############################################################


#' @title  runRCTD cell-type deconvolution
#' @author Dieter Henrik Heiland
#' @description Run the a cell type deconvolution using the SpaceXR/RCTD algorithm from SPATA objects
#' @inherit
#' @param object SPATA2 object
#' @param deconv_cell_types Character value of all celltypes (from the feature data.frame) should be used
#' @param spot_extension Numeric value between 0 and 0.7. Specifies extension of the spot radius to include more cells. For example, a spot extension of 0.2 will increase the spot radius of 20%)
#' @param multicore Logical. If TRUE will be run on all cores (recommended)
#' @param flip.y Logical. If TRUE y axis will be flipped
#' @param scale_factor Numeric value, if adjustment of the SPATAwrappers::getNucleusPosition() is required
#' @return
#' @examples
#' @export
#'


getSingleCellDeconv <- function(object,
                                deconv_cell_types,
                                spot_extension=NULL,
                                multicore=T,
                                flip.y=T,
                                workers=16,
                                scale_factor=NULL){

  SPATA2::check_method(object)

  if(spot_extension>0.7) stop("The spot extension will cause overlap in segmentation ")



  # First quantify the numbers of cell types in each spot and quantify the celltypes
  sc_dat <- SPATAwrappers::getNucleusPosition(object)

  #If required flipping and scaling of single cell coords
  if(flip.y==T){
    message("Y-axis will be fliped")
    yrange <- range(sc_dat$y)
    sc_dat$y <- yrange[2] - sc_dat$y + yrange[1]

  }
  if(!is.null(scale_factor)){
    sc_dat[,c("x","y")]=sc_dat[,c("x", "y")]*scale_factor


  }

  #Get cells per spot

  #-define the spotsize radius
  d <- SPATAwrappers::getSurroundedSpots(object) %>% filter(distance!=0) %>%  pull(distance) %>% min()
  r=(d*c(55/100))/2

  if(!is.null(spot_extension)){r=r+(r*spot_extension)}

  grid.plot <- SPATA2::getCoordsDf(object)


  if(multicore==T){

    # Run multicore
    base::options(future.fork.enable = TRUE)
    future::plan("multiprocess", workers = workers)
    future::supportsMulticore()
    base::options(future.globals.maxSize = 600 * 1024^2)
    message("... Run multicore ... ")

    #plot(sc_dat$x, sc_dat$y, pch=".")
    segments <- furrr::future_map_dfr(.x=1:nrow(grid.plot),
                                      .f=function(x){



                                        #print(x)
                                        #Create Segment
                                        segment <- swfscMisc::circle.polygon(x=grid.plot$x[x],
                                                                             y=grid.plot$y[x],
                                                                             radius=r,
                                                                             poly.type= "cartesian") %>% as.data.frame()
                                        #polygon(segment, border="red")

                                        #Get Objects in Segment
                                        nuc <- sp::point.in.polygon(pol.x=segment$x,
                                                                    pol.y=segment$y,
                                                                    point.x=sc_dat$x,
                                                                    point.y=sc_dat$y)

                                        cells <- sc_dat[nuc==1, ]$Cell

                                        if(is_empty(cells)){
                                          return <- data.frame(barcodes = grid.plot$barcodes[x],
                                                               Nr_of_cells = 0, cells = NA)
                                        }else{
                                          cells_in_spot <- sum(nuc)
                                          return <- data.frame(barcodes = grid.plot$barcodes[x],
                                                               Nr_of_cells = cells_in_spot, cells = cells)
                                        }
                                        return(return)


                                      },
                                      .progress = T)



  }else{

    #plot(sc_dat$x, sc_dat$y, pch=".")
    segments <- map_dfr(.x=nrow(grid.plot), .f=function(x){



      #print(x)
      #Create Segment
      segment <- swfscMisc::circle.polygon(x=grid.plot$x[x],
                                           y=grid.plot$y[x],
                                           radius=r,
                                           poly.type= "cartesian") %>% as.data.frame()
      polygon(segment, border="red")

      #Get Objects in Segment
      nuc <- sp::point.in.polygon(pol.x=segment$x,
                                  pol.y=segment$y,
                                  point.x=sc_dat$x,
                                  point.y=sc_dat$y)

      cells <- sc_dat[nuc==1, ]$Cell
      cells_in_spot<-sum(nuc)
      return <-
        data.frame(barcodes=grid.plot$barcodes[x],
                   Nr_of_cells=cells_in_spot,
                   cells=cells)



      return(return)


    })

  }


  deconv_df <- SPATA2::getFeatureDf(object) %>% dplyr::select(barcodes,{{deconv_cell_types}})

  # Run a cell type quantification per spot
  Cell_types <- furrr::future_map_dfr(.x=1:nrow(grid.plot),
                                      .f=function(x){

                                        #print(x)

                                        out <-
                                          segments %>%
                                          dplyr::filter(barcodes==grid.plot$barcodes[x])

                                        #Get the best cells
                                        nr.cells <- nrow(out)

                                        scores <-
                                          deconv_df %>%
                                          dplyr::filter(barcodes==grid.plot$barcodes[x]) %>%
                                          dplyr::select(-barcodes) %>%
                                          t() %>% as.data.frame() %>%
                                          dplyr::arrange(desc(V1)) %>%
                                          dplyr::mutate(score=scales::rescale(V1, c(0,1)) %>% round(.,digits = 3)) %>%
                                          dplyr::mutate(cells=round(score*c(nr.cells/sum(score))+0.5 ,digits = 0)) %>%
                                          dplyr::select(-V1)

                                        i <- 1
                                        cells_select=0
                                        while (cells_select < nr.cells) {
                                          cells_select <- cells_select+scores$cells[i]
                                          #print(cells_select)
                                          i=i+1
                                        }
                                        scores <- scores[1:i-1, ]

                                        if(sum(scores$cells)>nr.cells){
                                          neg <- sum(scores$cells)-nr.cells
                                          scores$cells[nrow(scores)]= scores$cells[nrow(scores)]-neg
                                        }

                                        cells_add <- map(.x=1:nrow(scores),.f=function(i){
                                          rep(rownames(scores)[i], scores$cells[i])}) %>% unlist()

                                        out$celltypes <- cells_add

                                        return(out)
                                      },
                                      .progress = T)


  #Align Coords

  Cell_types <- Cell_types %>% dplyr::left_join(., sc_dat %>% dplyr::rename("cells":=Cell), by="cells")

  return(Cell_types)

}



###############################################################
###########   Seurat add on Functions
###############################################################

#' @title  Run single-cell mapping
#' @author Dieter Henrik Heiland
#' @description Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset
#' @return
#' @examples
#' @export
DownScaleSeurat <- function(seurat,
                            maintain_var="cell_states",
                            max=10000,
                            min=50,
                            only_var=T,
                            factor=5,
                            n_feature=3000){


  tab_quant <-
    seurat@meta.data %>%
    as.data.frame() %>%
    pull(!!sym(maintain_var)) %>%
    table() %>%
    as.data.frame() %>%
    mutate(nr=Freq/min(Freq)) %>%
    mutate(cells = round(nr*factor)) %>%
    mutate(cells = ifelse(cells>max, max, cells)) %>%
    mutate(cells = ifelse(cells<min, Freq, cells)) %>%
    dplyr::rename("cluster":=.) %>%
    dplyr::select(cluster,cells) %>%
    mutate(cluster=as.character(cluster))

  cells <- map(.x=1:nrow(tab_quant), .f=function(i){
    sample(rownames(seurat@meta.data[seurat@meta.data[,maintain_var]==tab_quant$cluster[i], ]),tab_quant$cells[i])
  }) %>% unlist()

  seurat.out <- subset(seurat, cells=cells)
  seurat.out <- Seurat::SCTransform(seurat.out, return.only.var.genes=only_var, variable.features.n = n_feature)

}






###############################################################
###########   Full Mapping
###############################################################

#' @title  Run single-cell mapping
#' @author Dieter Henrik Heiland
#' @description Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset
#' @inherit
#' @param object SPATA2 object
#' @param scDF Data.frame; Output of the getSingleCellDeconv()
#' @param reference Seurat Object; Reference dataset used for runRCTD()
#' @param cell_type_var Character value; The col of the Seurat meta data indicating the cell type annotations
#' @param metaSpace The metaSpace definition from the function "defineMetaSpace"
#' @param workers Integer, Number of cors
#' @param iter Integer: For model BF: Number of random spot compositions; For model oGA size of initial population
#' @param iter_GA Integer: Number of iterations of the oGA model
#' @param nr_mut Integer: Number of mutations
#' @param nr_offsprings Integer: Number of offsprings
#' @param cross_over_point Numeric value: Percentage of cross-over cutt-off
#' @param ram Integer: GB of ram can be used for multicore session
#' @return
#' @examples
#' @export
#'


runMappingGA <- function (object,
                          scDF,
                          reference,
                          cell_type_var,
                          metaSpace,
                          workers = 16,
                          iter = 200,
                          iter_GA = 20,
                          nr_mut = 2,
                          nr_offsprings = 7,
                          cross_over_point = 0.5,
                          ram = 8)
{
  message(paste0("Start: ", print(Sys.time())))
  spots <- scDF %>% filter(Nr_of_cells>0) %>% pull(barcodes) %>% unique()
  message("--- Load data ----")
  mat.ref <- reference %>% Seurat::GetAssayData()
  genes.ref <- rownames(mat.ref)
  object <- SPATA2::setActiveExpressionMatrix(object, "scaled")
  mat.spata <- SPATA2::getExpressionMatrix(object)
  genes.spata <- rownames(mat.spata)
  metaSpace.genes <- as.data.frame(do.call(rbind, metaSpace))
  genes <- unique(metaSpace.genes$gene)
  message("--- Merge Data ----")
  mat.ref <- mat.ref[genes, ]
  mat.spata <- mat.spata[genes, ]
  message("--- Run Mapping ----")
  ref_meta <- reference@meta.data[colnames(mat.ref), c("nCount_RNA", cell_type_var)] %>% as.data.frame()
  nested_ref_meta <- ref_meta %>% rownames_to_column("cells") %>% group_by(!!sym(cell_type_var)) %>% nest()

  ### Loop
  gc()
  base::options(future.fork.enable = TRUE)
  future::plan("multisession", workers = 3)
  future::supportsMulticore()
  base::options(future.globals.maxSize = 600* 10* 1024^2)
  message("... Run multicore ... ")
  data.new <- furrr::future_map(.x=1:length(spots),.f=function(aa){
    #print(aa)
    bc_run <- spots[aa]
    # Cells in spot
    data <-
      scDF %>%
      dplyr::filter(barcodes==bc_run) %>%
      dplyr::mutate(celltypes=as.character(celltypes)) %>%
      dplyr::arrange(celltypes) %>%
      count(celltypes) %>%
      ungroup() %>%
      as.data.frame()
    nr_cells <- sum(data$n)


    # 1. Get genes of cell population
    genes.use <- as.data.frame(do.call(rbind, metaSpace[data$celltypes]))$gene

    # 2. Remove genes that are not helpful
    exp.spot <- mat.spata[genes.use, bc_run]
    mat.ref.use <- mat.ref[genes.use, ] %>% as.matrix()

    dim(mat.ref.use)[1]==length(exp.spot)
    type <- reference@meta.data[colnames(mat.ref.use), cell_type_var]
    names(type) <- colnames(mat.ref.use)

    # 3. Initiate a Population of cells with similar distrubution

    names(nested_ref_meta)[1] <- "celltypes"
    select <- nested_ref_meta[nested_ref_meta$celltypes %in% data$celltypes, ]




    if(length(unique(data$celltypes)) %in% c(1,2)){
      initial.pop <- map(.x=1:1000, .f=function(a){
        map(.x=1:nrow(select), .f=function(.x){
          s <- NA
          while(any(is.na(s)==T)){
            s <- sample(select$data[[.x]]$cells, size = data$n[.x])
          }
          return(s)

        }) %>% unlist()
      })
      #print(map(.x=initial.pop, .f=~length(.x)) %>% unlist())
      #print(map(.x=initial.pop, .f=~type[names(type)%in% .x]))

      fitness <- map(.x=initial.pop, .f=function(.x){
        .x <- .x[.x %in% colnames(mat.ref.use)]
        if(length(.x)==1){out <- cor(mat.ref.use[,.x], exp.spot)}else{
          out <- cor(rowMeans(mat.ref.use[,.x]), exp.spot)
        }
      } ) %>% unlist()
      names(fitness) <- 1:1000
      fitness <- fitness[order(-fitness)]
      top <- initial.pop[names(fitness[1:data$n %>% sum()]) %>% as.numeric()]
    }else{

      #init pop
      initial.pop <- map(.x=1:iter, .f=function(a){
        map(.x=1:nrow(select), .f=function(.x){
          s <- NA
          while(any(is.na(s)==T)){
            s <- sample(select$data[[.x]]$cells, size = data$n[.x])
          }
          return(s)

        }) %>% unlist()
      })

      if(length(unique(data$celltypes))>5){mutation.vec <- seq(3,1,length.out=iter_GA) %>% round()}
      if(length(unique(data$celltypes))<=5){mutation.vec <- rep(1,iter_GA)}
      for(z in 1:iter_GA){
        #print(z)
        fitness <- map(.x=initial.pop, .f=function(.x){
          .x <- .x[.x %in% colnames(mat.ref.use)]
          if(length(.x)==1){out <- cor(mat.ref.use[,.x], exp.spot)}else{
            out <- cor(rowMeans(mat.ref.use[,.x]), exp.spot)
          }
        } ) %>% unlist()
        names(fitness) <- 1:iter
        fitness <- fitness[order(-fitness)]
        parents <- initial.pop[names(fitness[1:2]) %>% as.numeric()]

        cross <- parents

        # 5. Mutation
        nr_mut=mutation.vec[z]
        offspring <- map(.x=1:7, .f=function(.x){
          
          # random parent
          gender <- runif(1, 1, 2) %>% round()
          
          mut <- cross[gender][[1]]
          non_mutual <- type[!c(names(type) %in% mut)]
          
          mut.pos <- sample(mut, nr_mut)
          subtype <- type[names(type) %in% mut.pos]
          
          # get the random new cells
          new.cells <- map(1:nr_mut, .f=function(.x){
            if(length(non_mutual[non_mutual == subtype[.x]])==0){
              sample(type[type == subtype[.x]], 1)
            }else{
              sample(non_mutual[non_mutual == subtype[.x]], 1)
            }
            
          }) %>% unlist()
          
          final <- c(mut[!c(mut %in% mut.pos)], names(new.cells))
          
          return(final)
        })

        #map(.x=offspring, .f= ~ unique(type[names(type) %in% .x]) %in% select$celltypes)

        # 6. Remove the tail and add offsprings
        remove <- names(fitness[c(length(fitness)-length(offspring)+1):length(fitness)]) %>% as.numeric()
        initial.pop[remove] <- offspring

      }
      # After GA loop eval. top comb
      fitness <- map(.x=initial.pop, .f=function(.x){
        .x <- .x[.x %in% colnames(mat.ref.use)]
        if(nr_cells==1){out <- cor(mat.ref.use[,.x], exp.spot)}else{
          out <- cor(rowMeans(mat.ref.use[,.x]), exp.spot)
        }
      } ) %>% unlist()
      names(fitness) <- 1:iter
      fitness <- fitness[order(-fitness)]
      top <- initial.pop[names(fitness[1]) %>% as.numeric()] %>% unlist()

    }



    data.out <-
      scDF %>%
      dplyr::filter(barcodes==bc_run) %>%
      dplyr::mutate(celltypes=as.character(celltypes)) %>%
      dplyr::arrange(celltypes)
    diff <- nrow(data.out)-length(top)
    if(diff!=0){top <- top %>% unlist()}
    diff <- nrow(data.out)-length(top)
    if(diff!=0){top <- c(top,rep("Empty", diff))}

    data.out <-
      data.out %>%
      dplyr::mutate(match=top)#, ref.type=type[names(type)%in%top])





  },.options = furrr::furrr_options(seed = TRUE), .progress=T)

  message(paste0("End: ", print(Sys.time())))
  return(data.new)
}

#' @title  Run single-cell mapping
#' @author Dieter Henrik Heiland
#' @description Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset
#' @inherit
#' @param object SPATA2 object
#' @param scDF Data.frame; Output of the getSingleCellDeconv()
#' @param cell_type_var Character value; The col of the Seurat meta data indicating the cell type annotations
#' @param MetaVol Number of genes in meta space
#' @param ref.new Seurat Object; Reference dataset used for runRCTD()
#' @return
#' @examples
#' @export
#'
defineMetaSpace <- function(object,
                            reference,
                            cell_type_var="annotation_level_4",
                            scDF,
                            MetaVol=20){
  cell_types <- scDF
  # getSPATA genes
  st.genes <- object %>% SPATA2::getGenes()

  # Subset shared genes
  genes.shared <- intersect(reference@assays$RNA@counts %>% rownames(), st.genes)
  reference <- subset(reference, features=genes.shared)

  #Subset annotations
  annotations <- unique(cell_types$celltypes)
  cells <- reference@meta.data %>% as.data.frame() %>% filter(!!sym(cell_type_var) %in% annotations) %>% rownames()
  reference <- subset(reference, cells=cells)

  message("Get DE per cluster")
  Idents(reference) <- cell_type_var
  diff_gene_exp <- Seurat::FindAllMarkers(reference, max.cells.per.ident = 500)

  message("Get PC1 for variance")
  pca <- map(.x=1:length(annotations), .f=function(i){

    cells <- reference@meta.data %>% as.data.frame() %>% filter(!!sym(cell_type_var)==annotations[i] ) %>% rownames()
    message(paste0("Size of the explored subcluster: ", length(cells)))
    sub <- subset(reference, cells=cells)

    sub <- sub %>% Seurat::RunPCA(verbose=F, npcs=2)
    pca1 <- sub@reductions$pca@feature.loadings[,1] %>% as.data.frame()
    names(pca1) <- "PCA1"
    pca1$cluster <- annotations[i]




    return(pca1)

  })
  names(pca) <- annotations

  # Get the top genes of spe/var

  optimal <- map(.x=1:length(annotations), .f=function(i){

    anno <- as.character(annotations[i])
    DE <- diff_gene_exp %>% filter(cluster==anno)
    pca.genes <- pca[[anno]] %>% rownames_to_column("gene")

    merge <- DE %>% dplyr::select(gene, avg_log2FC) %>% left_join(., pca.genes, "gene") %>% na.omit()
    merge <- merge %>% arrange(desc(avg_log2FC),desc(PCA1)) %>% head(MetaVol)

    return(merge)

  })
  names(optimal) <- annotations

  return(optimal)

}


###############################################################
###########   GA Functions
###############################################################

#' @title  fitness
#' @author Dieter Henrik Heiland
#' @description fitness
#' @inherit
#' @return
#' @examples
#' @export
#'
fitness <- function(x,nr_cells){

  if(nr_cells==1){
    y <- cor(as.numeric(mat.ref[,x==1]),
             as.numeric(mat.spata[,bc_run]))
  }else{
    y <- cor(as.numeric(as.matrix(mat.ref[,x==1]) %>% rowMeans()),
             as.numeric(mat.spata[,bc_run]))
  }
  return(y)
}

#' @title  initiate_Population
#' @author Dieter Henrik Heiland
#' @description initiate_Population
#' @inherit
#' @return
#' @examples
#' @export
#'
#'
initiate_Population <- function(nr_of_random_spots,
                                n_select,
                                nested_ref_meta,
                                cell_type_var){

  select <- nested_ref_meta %>% filter(!!sym(cell_type_var) %in%  n_select$celltypes)

  cells.select <- map(1:nrow(select), .f= ~select$data[[.x]]$cells) %>% unlist()
  mat.select <- mat.ref[,cells.select]
  dim(mat.select)

  mat <- matrix(0, nrow=nr_of_random_spots, ncol=ncol(mat.select))
  colnames(mat)=colnames(mat.select)


  for(x in 1:nr_of_random_spots){
    selected_cells <- map(.x=1:nrow(select), .f=function(i){
      sample(select$data[i] %>% as.data.frame() %>% pull(cells), n_select$n[i])
    }) %>% unlist()
    mat[x,selected_cells] <- 1
  }

  return(mat)

}

#' @title  cross_over
#' @author Dieter Henrik Heiland
#' @description cross_over
#' @inherit
#' @return
#' @examples
#' @export
#'
cross_over <- function(parents,cross_over_point=0.5){

  cross_over_select <- c(ncol(parents)*cross_over_point) %>% round()

  cross <- parents
  cross[1, cross_over_select:ncol(parents)] <- parents[2, cross_over_select:ncol(parents)]
  cross[2, cross_over_select:ncol(parents)] <- parents[1, cross_over_select:ncol(parents)]

  return(cross)
}

#' @title  mutation_GA
#' @author Dieter Henrik Heiland
#' @description mutation_GA
#' @inherit
#' @return
#' @examples
#' @export
#'
mutation_GA <- function(offspring_2,nr_mut, nr_offsprings){

  #Create random selection
  selector <- runif(nr_offsprings, 1, 2) %>% round(digits = 0)

  # Create a mutated and non-mutated output
  offsprings_mut <- matrix(0, ncol = ncol(offspring_2), nrow=nr_offsprings)
  colnames(offsprings_mut) <- colnames(offspring_2)


  #Run loop for offsprings

  for(u in 1:nr_offsprings){


    #Select random cell and mutate from same cell type
    cells_off <- which(offspring_2[selector[u], ]==1) %>% names()
    cells_mut <- sample(cells_off, nr_mut)
    # Get cell type

    ref_sub <- ref_meta[!rownames(ref_meta) %in% cells_off, ]
    cell_type_select <- ref_meta[cells_mut, cell_type_var]

    new <- list()
    for(z in 1:nr_mut){
      new[[z]] <- sample(ref_sub %>%
                           filter(!!sym(cell_type_var)==cell_type_select[z]) %>%
                           rownames(),1)

    }

    offsprings_inter <- offspring_2

    #message(length(cells_mut)==length(unlist(new)))

    offsprings_inter[selector[u],cells_mut]=0
    offsprings_inter[selector[u],unlist(new)]=1

    #message(length(which(offsprings_inter[selector[u], ]==1)))
    #message(any((cells_off %in% unlist(new)) == T))

    offsprings_mut[u, ] <- offsprings_inter[selector[u], ]


  }

  return(offsprings_mut)
}



#' @title  fitness_V1
#' @author Dieter Henrik Heiland
#' @description fitness_V1
#' @inherit
#' @return
#' @examples
#' @export
#'
fitness_V1 <- function(x,nr_cells){

  if(nr_cells==1){
    y <- cor(as.numeric(mat.ref[,x==1]),
             as.numeric(mat.spata[,bc_run]))
  }else{
    y <- cor(as.numeric(mat.ref[,x==1] %>% rowMeans()),
             as.numeric(mat.spata[,bc_run]))
  }
  return(y)
}


#' @title  initiate_Population_V1
#' @author Dieter Henrik Heiland
#' @description initiate_Population_V1
#' @inherit
#' @return
#' @examples
#' @export
#'
initiate_Population_V1 <- function(nr_of_random_spots, n_select,nested_ref_meta){

  mat <- matrix(0, nrow=nr_of_random_spots, ncol=ncol(mat.ref))
  colnames(mat)=colnames(mat.ref)
  select <- nested_ref_meta %>% filter(!!sym(cell_type_var) %in%  n_select$celltypes)

  for(x in 1:nr_of_random_spots){
    selected_cells <- map(.x=1:nrow(select), .f=function(i){
      sample(select$data[i] %>% as.data.frame() %>% pull(cells), n_select$n[i])
    }) %>% unlist()
    mat[x,selected_cells] <- 1
  }

  return(mat)

}

#' @title  cross_over_V1
#' @author Dieter Henrik Heiland
#' @description cross_over_V1
#' @inherit
#' @return
#' @examples
#' @export
#'
cross_over_V1 <- function(parents,cross_over_point=0.5){

  parents_out <-
    parents %>%
    apply(1, function(x) which(x == 1)) %>% t()
  #as.data.frame() %>%
  #filter(!V1 %in% intersect(V1,V2)) %>%
  #filter(!V2 %in% intersect(V1,V2)) %>%
  #t()

  #message(length(intersect(parents_out[1,], parents_out[2,])))

  cross_over_select <- c(ncol(parents_out)*cross_over_point) %>% round()
  parents_out <- parents_out[, 1:cross_over_select]


  # cross values
  parents[1, parents_out[1,]]=0
  parents[1, parents_out[2,]]=1

  parents[2, parents_out[2,]]=0
  parents[2, parents_out[1,]]=1

  #message(length(which(parents[1, ]==1)))
  #message(length(which(parents[2, ]==1)))



  return(parents)
}

#' @title  mutation_GA_V1
#' @author Dieter Henrik Heiland
#' @description mutation_GA_V1
#' @inherit
#' @return
#' @examples
#' @export
#'
mutation_GA_V1 <- function(offspring_2,nr_mut, nr_offsprings){

  #Create random selection
  selector <- runif(nr_offsprings, 1, 2) %>% round(digits = 0)

  # Create a mutated and non-mutated output
  offsprings_mut <- matrix(0, ncol = ncol(offspring_2), nrow=nr_offsprings)
  colnames(offsprings_mut) <- colnames(offspring_2)


  #Run loop for offsprings

  for(u in 1:nr_offsprings){


    #Select random cell and mutate from same cell type
    cells_off <- which(offspring_2[selector[u], ]==1) %>% names()
    cells_mut <- sample(cells_off, nr_mut)
    # Get cell type

    ref_sub <- ref_meta[!rownames(ref_meta) %in% cells_off, ]
    cell_type_select <- ref_meta[cells_mut, cell_type_var]

    new <- list()
    for(z in 1:nr_mut){
      new[[z]] <- sample(ref_sub %>%
                           filter(!!sym(cell_type_var)==cell_type_select[z]) %>%
                           rownames(),1)

    }

    offsprings_inter <- offspring_2

    #message(length(cells_mut)==length(unlist(new)))

    offsprings_inter[selector[u],cells_mut]=0
    offsprings_inter[selector[u],unlist(new)]=1

    #message(length(which(offsprings_inter[selector[u], ]==1)))
    #message(any((cells_off %in% unlist(new)) == T))

    offsprings_mut[u, ] <- offsprings_inter[selector[u], ]


  }

  return(offsprings_mut)
}



###############################################################
###########   Other Functions
###############################################################

#' @title  Run single-cell mapping
#' @author Dieter Henrik Heiland
#' @description Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset
#' @inherit
#' @param tab Input data of the class data.frame
#' @param class Character value; The col containing the main class
#' @param subclass Character value; The col containing the subclass
#' @param pal Character value; Color pal from brewer.pal()
#' @param random Logical. If TRUE random mixing colors from pal
#' @param seed Integer value. set.seed() for constant color alignment
#' @param into Character value; The color for non-classified samples
#' @param add_n Integer value. Adopting the contrast of the subclass colors
#' @return
#' @examples
#' @export
#'
getSubColors <- function(tab,
                         class="annotation_level_2",
                         subclass="annotation_level_4",
                         pal="Set3",
                         max_pal=12,
                         random=T,
                         seed=200,
                         into="#EDEDED",
                         add_n=1){


  out <-
    tab %>%
    as.data.frame() %>%
    dplyr::group_by(!!sym(class), !!sym(subclass)) %>%
    dplyr::summarise(n=length(!!sym(class)))

  class_unique <- unique(tab %>% pull(!!sym(class))) %>% as.character()

  if(length(class_unique)>max_pal){
    color <- RColorBrewer::brewer.pal(max_pal,pal)
    color_L1 <- colorRampPalette(color)(length(class_unique))
  }else{
    color_L1 <- RColorBrewer::brewer.pal(length(class_unique),pal)
  }

  if(random==T){
    set.seed(seed)
    color_L1 <- sample(color_L1)}


  out2 <- map_dfr(.x=1:length(class_unique), .f=function(i){
    x <-
      out %>%
      dplyr::filter(!!sym(class)==class_unique[i])
    a <- nrow(x)
    x <-
      x %>%
      dplyr::mutate(colors=c(colorRampPalette(color=c(into,color_L1[i]))(a+add_n)[c(1+add_n):c(a+add_n)] %>% rev()) )
  })

  all <- out %>% dplyr::pull(!!sym(subclass)) %>% unique()
  withcolor <- out2 %>% dplyr::pull(!!sym(subclass)) %>% unique()
  inter <- intersect(withcolor, all)

  if(length(all[!all %in% inter])>0){
    out3 <-data.frame(a="NaN", b=all[!all %in% inter], n=1, colors=into)
    names(out3) <- names(out2)
    out_4 <- rbind(out2,out3) %>% dplyr::ungroup()
  }else{
    out_4 <- out2
  }


  return(out_4)

}


###############################################################
###########   Infer Cell Position
###############################################################

#' @title  Run single-cell mapping
#' @author Dieter Henrik Heiland
#' @description Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset
#' @inherit
#' @param object SPATA object
#' @param sample.folder path to Visium sample folder
#' @return SPATA object
#' @examples
#' @export
#'

inferCellPositionsfromHE <- function(object,sample.folder){

  message("Create ST Object")
  sample_name <- SPATA2::getSampleNames(object)


  input.df <- data.frame(samples=paste0(sample.folder, "/outs/filtered_feature_bc_matrix.h5"),
                         spotfiles=paste0(sample.folder, "/outs/spatial/tissue_positions_list.csv"),
                         imgs=paste0(sample.folder, "/outs/spatial/tissue_hires_image.png"),
                         json=paste0(sample.folder, "/outs/spatial/scalefactors_json.json"))
  library(imager)
  library(Seurat)
  se <- STutility::InputFromTable(infotable = input.df,
                                  min.gene.count = 1,
                                  min.gene.spots = 1,
                                  min.spot.count = 5,
                                  platform =  "Visium")


  se <- STutility::LoadImages(se, time.resolve = FALSE, xdim=1000, verbose = F)
  se <- STutility::MaskImages(se, verbose = F)
  se <- STutility::AlignImages(se, verbose = F)

  message("Infer Cell Position")
  se <- STutility::Create3DStack(se, verbose = F)
  stack_3d <- setNames(STutility::GetStaffli(se)@scatter.data, c("x", "y", "z", "grid.cell"))
  img <- EBImage::Image(se@tools$Staffli@rasterlists$processed$`1` %>% as.matrix() %>% t(), colormode = "Color")
  #plot(img)
  #points(stack_3d$x, stack_3d$y)

  dim.org <- object@images[[sample_name]]@image %>% dim()
  dim.proc <- img %>% dim()
  scale.f <- dim.org[1]/dim.proc[1]
  message("Transfer to SPATA")
  nuc.data <-
    stack_3d %>%
    dplyr::mutate(Cell=paste0("Cell_", 1:nrow(.))) %>%
    dplyr::select(Cell,x,y) %>%
    dplyr::mutate(x=x*scale.f,
                  y=y*scale.f)
  object@spatial[[sample_name]]$Cell_coords <- nuc.data

  return(object)
}









###############################################################
###########   Image or DF sum per SPOT
###############################################################

#' @title  runSegmentfromImage
#' @author Dieter Henrik Heiland
#' @description Get summary of image intensity per spot
#' @inherit
#' @param object SPATA object
#' @param Image EBImage object
#' @param spot_extension Extent the area of the spot (from SPATA2 Image processing)
#' @param multicore Use multicore !!recommended!!
#' @param workers Number of workers for multisession
#' @return data.frame with intensity per spot
#' @examples
#' @export
#'
runSegmentfromImage<- function(object, Image, spot_extension=0, multicore=T, workers = 16){

  #message fine correct spot size

  SPATA2::check_method(object)
  if (spot_extension > 0.7) stop("The spot extension will cause overlap in segmentation ")

  # Get spot radius
  getSpotRadius <- function(object){
    of_sample <- SPATA2::getSampleNames(object)
    coords <- SPATA2::getCoordsDf(object)
    bc_origin <- coords$barcodes
    bc_destination <- coords$barcodes
    d <-
      tidyr::expand_grid(bc_origin, bc_destination) %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
      dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
      dplyr::filter(distance!=0) %>%
      dplyr::pull(distance) %>%
      min()

    r = (d * c(55/100))/2

    return(r)
  }

  r <- getSpotRadius(object)
  if (!is.null(spot_extension)) {r = r + (r * spot_extension)}
  grid.plot <- SPATA2::getCoordsDf(object)

  dim(Image)
  if(dim(Image)[3]!=1){ Image <- Reduce(`+`, map(.x=1:dim(Image)[3], ~Image[,,.x])) }

  sc_dat <- reshape2::melt(Image)
  head(sc_dat)

  names(sc_dat)=c("x","y","var")

  yrange <- range(sc_dat$y)
  sc_dat$y <- (yrange[2] - sc_dat$y) + yrange[1]


  if (multicore == T) {
    base::options(future.fork.enable = TRUE)
    future::plan("multiprocess", workers = workers)
    future::supportsMulticore()
    base::options(future.globals.maxSize = 600 * 1024^2)
    message("... Run multicore ... ")

    #plot(sc_dat$x, sc_dat$y, pch = ".")
    segments <- furrr::future_map_dfr(.x = 1:nrow(grid.plot),
                                      .f = function(x) {
                                        segment <- swfscMisc::circle.polygon(x = grid.plot$x[x],
                                                                             y = grid.plot$y[x], radius = r, poly.type = "cartesian") %>% as.data.frame()


                                        nuc <- sp::point.in.polygon(pol.x = segment$x,
                                                                    pol.y = segment$y,
                                                                    point.x = sc_dat$x,
                                                                    point.y = sc_dat$y)


                                        intensity <- sc_dat[nuc == 1, ]$var

                                        out <- data.frame(barcodes = grid.plot$barcodes[x],
                                                          mean=mean(intensity),
                                                          median=median(intensity),
                                                          max=max(intensity),
                                                          min=min(intensity),
                                                          sd=sd(intensity))

                                        return(out)
                                      }, .progress = T)
  }
  else {
    segments <- map_dfr(.x = 1:nrow(grid.plot), .f = function(x) {
      segment <- swfscMisc::circle.polygon(x = grid.plot$x[x],
                                           y = grid.plot$y[x], radius = r, poly.type = "cartesian") %>%
        as.data.frame()
      polygon(segment, border = "red")
      nuc <- sp::point.in.polygon(pol.x = segment$x, pol.y = segment$y,
                                  point.x = sc_dat$x, point.y = sc_dat$y)
      intensity <- sc_dat[nuc == 1, ]$var

      out <- data.frame(barcodes = grid.plot$barcodes[x],
                        mean=mean(intensity),
                        median=median(intensity),
                        max=max(intensity),
                        min=min(intensity),
                        sd=sd(intensity))

      return(out)
    })
  }

  # Summarize on barcode level

  return(segments)

}

#' @title  runSegmentfromCoords
#' @author Dieter Henrik Heiland
#' @description Get summary of image intensity per spot
#' @inherit
#' @param object SPATA object
#' @param Coord_file data.frame with cell or other coordinates if NULL the SPATAwrappers::getNucleusPosition will look into saves scCoords
#' @param spot_extension Extent the area of the spot (from SPATA2 Image processing)
#' @param multicore Use multicore !!recommended!!
#' @param workers Number of workers for multisession
#' @return data.frame with intensity per spot
#' @examples
#' @export
#'
runSegmentfromCoords <- function(object, Coord_file=NULL, spot_extension=0, multicore=T,workers = 16){

  #message fine correct spot size

  SPATA2::check_method(object)
  if (spot_extension > 0.7) stop("The spot extension will cause overlap in segmentation ")
  if(is.null(Coord_file)){sc_dat <- SPATAwrappers::getNucleusPosition(object)}else{
    if (!any(names(Coord_file) %in% c("Cell","x","y"))) stop("Coord_file does not contain the variables: Cell,x,y ... ")
    sc_dat <- Coord_file[, c("Cell","x","y")]
  }

  # Get spot radius
  getSpotRadius <- function(object){
    of_sample <- SPATA2::getSampleNames(object)
    coords <- SPATA2::getCoordsDf(object)
    bc_origin <- coords$barcodes
    bc_destination <- coords$barcodes
    d <-
      tidyr::expand_grid(bc_origin, bc_destination) %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
      dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
      dplyr::filter(distance!=0) %>%
      dplyr::pull(distance) %>%
      min()

    r = (d * c(55/100))/2

    return(r)
  }

  r <- getSpotRadius(object)
  if (!is.null(spot_extension)) {r = r + (r * spot_extension)}
  grid.plot <- SPATA2::getCoordsDf(object)


  if (multicore == T) {
    base::options(future.fork.enable = TRUE)
    future::plan("multiprocess", workers = workers)
    future::supportsMulticore()
    base::options(future.globals.maxSize = 600 * 1024^2)
    message("... Run multicore ... ")

    #plot(sc_dat$x, sc_dat$y, pch = ".")
    segments <- furrr::future_map_dfr(.x = 1:nrow(grid.plot),
                                      .f = function(x) {
                                        segment <- swfscMisc::circle.polygon(x = grid.plot$x[x],
                                                                             y = grid.plot$y[x], radius = r, poly.type = "cartesian") %>%
                                          as.data.frame()
                                        nuc <- sp::point.in.polygon(pol.x = segment$x,
                                                                    pol.y = segment$y, point.x = sc_dat$x, point.y = sc_dat$y)
                                        cells <- sc_dat[nuc == 1, ]$Cell
                                        if(is_empty(cells)){
                                          return <- data.frame(barcodes = grid.plot$barcodes[x],
                                                               Nr_of_cells = 0, cells = NA)
                                        }else{
                                          cells_in_spot <- sum(nuc)
                                          return <- data.frame(barcodes = grid.plot$barcodes[x],
                                                               Nr_of_cells = cells_in_spot, cells = cells)
                                        }
                                        return(return)
                                      }, .progress = T)
  }
  else {
    plot(sc_dat$x, sc_dat$y, pch = ".")
    segments <- map_dfr(.x = 1:nrow(grid.plot), .f = function(x) {
      segment <- swfscMisc::circle.polygon(x = grid.plot$x[x],
                                           y = grid.plot$y[x], radius = r, poly.type = "cartesian") %>%
        as.data.frame()
      polygon(segment, border = "red")
      nuc <- sp::point.in.polygon(pol.x = segment$x, pol.y = segment$y,
                                  point.x = sc_dat$x, point.y = sc_dat$y)
      cells <- sc_dat[nuc == 1, ]$Cell
      if(is_empty(cells)){
        return <- data.frame(barcodes = grid.plot$barcodes[x],
                             Nr_of_cells = 0, cells = NA)
      }else{
        cells_in_spot <- sum(nuc)
        return <- data.frame(barcodes = grid.plot$barcodes[x],
                             Nr_of_cells = cells_in_spot, cells = cells)
      }

      return(return)
    })
  }

  # Summarize on barcode level

  out <-
    segments %>%
    filter(Nr_of_cells!=0) %>%
    group_by(barcodes) %>%
    summarise(Cells=length(barcodes)) %>%
    left_join(grid.plot %>% dplyr::select(barcodes), ., "barcodes") %>%
    replace_na(list(Cells = 0))

  out <- list(out, segments)
  names(out) <- c("Feature_cells", "DF_Segments")
  return(out)

}










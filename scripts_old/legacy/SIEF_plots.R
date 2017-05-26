# plots for Steve SIEF

#source('scripts/transformations.R')

# Calvin/Light Rxns vs % canopy openness (to replace vs mean annual irradiance)

#deprecated# agg_plot_save('calv_per_photo', 'gap_mean', logx=FALSE, proportion = TRUE, indepvarType = 'gap', 
             # labs =  c('Calvin cycle / photosystems abundance', 'Canopy openness (%)'), outDir = 'output/figures/20170217/tiff', fileType = 'tiff', goldenRatio = FALSE)



# RCA/PGK vs MAT or another version of T that is a stronger predictor (MAT is R2 0.28)

        source('scripts/transformations.R')
        
        source('scripts/prep_data_mg_per_mm2.R')
        
        # pull out particular proteins from protein_samples_D14 and arrange
        
        rbact <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.j02030.1.p', 'eucgr.j01234.1.p','eucgr.l03031.1.p','eucgr.b02310.1.p','eucgr.b02532.1.p'),] # rubisco activase, chaperonin
        rb <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.c03525.1.p', # rubisco large subunit
                                                                     'tr|t1qkk4|t1qkk4_eucgl'),]
        rb <- rbind(rb, rbact)
        rb <- rbind(rb, protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.e01261.1.p', 'eucgr.b01439.1.p', 'eucgr.f01476.1.p'),]) # PRK, PGLP, PGK
        
        rb <- rbind(rb, protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.j01502.1.p', 'eucgr.j01502.2.p', 'eucgr.c00150.1.p', 'eucgr.b03013.1.p', 'eucgr.k02223.1.p'),]) # small subunit of rubisco
        
        
        rb <- arrange(rb, Protein)
        names <- rb$Protein
        rb <- as.data.frame(t(rb[,2:314]))
        
        names(rb) <- names
        rb$sample <- rownames(rb)
        rownames(rb) <- NULL
        
        rb <- merge(rb, total_protein_D14)
        
        rb[,2:(ncol(rb)-1)] <- lapply(rb[,2:(ncol(rb)-1)],function(x) as.numeric(as.character(x)))
        
        #str(rb)
        
        # sum components of multi-protein complexes
        
        rb$rb_largesub <- rowSums(rb[,c('eucgr.c03525.1.p',
                                        'tr|t1qkk4|t1qkk4_eucgl'),])
        
        rb$rb_act <- rowSums(rb[,c('eucgr.j02030.1.p', 
                                   'eucgr.j01234.1.p',
                                   'eucgr.l03031.1.p',
                                   'eucgr.b02310.1.p')])
        
        
        # merge with full protein and env data sets
        
        bla <- merge(rb, protein_D14)
        bla <- merge(bla, climate_locs)
        bla <- merge(bla, replicates)
        
        # aggregate to ID-wise means
        
        data_means <- bla %>% group_by(ID) %>% summarise(Rubisco_mean = mean(Rubisco, na.rm=TRUE),
                                                         Calvin_cycle_mean = mean(Calvin_cycle, na.rm=TRUE),
                                                         Calvin_per_lightrxns_mean = mean((Calvin_cycle/Light_reactions), na.rm=TRUE),
                                                         Calvin_per_lightrxns_SE = SE(Calvin_cycle/Light_reactions),
                                                         PRK_mean = mean(eucgr.e01261.1.p, na.rm=TRUE),
                                                         PRK_SE = SE(eucgr.e01261.1.p),
                                                         PGLP_mean = mean(eucgr.b01439.1.p, na.rm=TRUE),
                                                         PGLP_SE = SE(eucgr.b01439.1.p),
                                                         PGK_mean = mean(eucgr.f01476.1.p, na.rm=TRUE),
                                                         PGK_SE = SE(eucgr.f01476.1.p),
                                                         Photorespiration_mean = mean(Photorespiration, na.rm=TRUE),
                                                         rbact_mean = mean(rb_act, na.rm=TRUE),
                                                         rbact_SE = SE(rb_act),
                                                         rbL_mean = mean(rb_largesub, na.rm=TRUE),
                                                         chaperonin_mean = mean(eucgr.b02532.1.p, na.rm=TRUE),
                                                         total_protein_mean = mean(total_protein),
                                                         Photosystems_mean = mean(Photosystems, na.rm=TRUE),
                                                         Light_rxns_mean = mean(Light_reactions, na.rm=TRUE),
                                                         rb_small1 = mean(eucgr.j01502.1.p, n.rm=TRUE),
                                                         rb_small2 = mean(eucgr.j01502.2.p, na.rm=TRUE),
                                                         rb_small3 = mean(eucgr.c00150.1.p, na.rm=TRUE),
                                                         rb_small4 = mean(eucgr.b03013.1.p, na.rm=TRUE),
                                                         rb_small5 = mean(eucgr.k02223.1.p, na.rm=TRUE),
                                                         rbact_per_PGK_mean = mean((rb_act/eucgr.f01476.1.p), na.rm=TRUE),
                                                         rbact_per_PGK_SE = SE(rb_act/eucgr.f01476.1.p),
                                                         Calv_ex_rbact_mean = mean((Calvin_cycle - Rubisco), na.rm=TRUE),
                                                         Calv_ex_rbact_SE = SE(Calvin_cycle - Rubisco)) %>%
          full_join(bla, by = 'ID') %>%
          mutate(Calv_ex_rbact = Calvin_cycle_mean - rbact_mean) %>%
          mutate(Calv_ex_rbact_ex_rbc = Calvin_cycle_mean - rbact_mean - Rubisco_mean) %>%
          mutate(Calv_ex_rbact_ex_rbc_ex_chap = Calv_ex_rbact_ex_rbc - chaperonin_mean)
        
  model <- summary(lm(rbact_per_PGK_mean ~ tavg, data_means))
    
  y = as.numeric(coef(model)[1]) + as.numeric(coef(model)[2]) * min(data_means$tavg)
  x = as.numeric(coef(model)[1]) + as.numeric(coef(model)[2]) * max(data_means$tavg)
  
  rise = x - y 
  
  run = max(data_means$tavg) - min(data_means$tavg)
  
  slope = rise/run
  
  outDir <- 'output/figures/20170217/tiff'
  depvar <-  'rbact_per_PGK_mean'
  indepvar <- 'tavg'

  tiff(paste(outDir, '/', depvar, '_vs_', indepvar, '_R2-', round(model$r.squared,3), '_pval-', round(model$coefficients[,4][2],2),'.tiff', sep = ""), height = 6, width = 6, units = 'in', res = 600)
    p <- ggplot(data_means, aes(y = rbact_per_PGK_mean, x = tavg)) 
    p <- p + geom_point(size = 4) 
    p <- p + geom_smooth(method = 'lm', size=0.5, colour = 'black', se=FALSE)
    p <- p + geom_errorbar(aes(ymin = rbact_per_PGK_mean - rbact_per_PGK_SE, ymax = rbact_per_PGK_mean + rbact_per_PGK_SE), width = 0.3)
    p <- p + expand_limits(y=0)
    p <- p + ylab('Rbc_act per PGK abundance') + xlab(expression(paste("Mean annual temperature (",degree,"C)", sep = "")))
    p <- p + theme_bw()
    p <- p + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   legend.title=element_blank(),
                   legend.position="bottom",
                   axis.line = element_line(colour = "black"),
                   text = element_text(size = 20))
    p
  dev.off()


# PGLP vs PRK colored by T. Could you please do a version without the CI? Also, units are mg/m^2, is that right?
  
  model <- summary(lm(PGLP_mean ~ PRK_mean, data_means))
  
  outDir <- 'output/figures/20170217/tiff'
  depvar <-  'PGLP_mean'
  indepvar <- 'PRK_mean'
  
  
  tiff(paste(outDir, '/', depvar, '_vs_', indepvar, '_R2-', round(model$r.squared,3), '_pval-', round(model$coefficients[,4][2],2),'.tiff', sep = ""), height = 6, width = 6, units = 'in', res = 600)
  p <- ggplot(data_means, aes(y = PGLP_mean, x = PGK_mean)) 
  p <- p + geom_point(size = 4, aes(colour = tavg)) 
  p <- p + scale_colour_gradientn(colours=c('cyan', 'red'), name = expression(paste("Mean annual temperature (",degree,"C)", sep = "")))
  p <- p + geom_smooth(method = 'lm', size=1, colour = 'blue', se=FALSE)
  p <- p + geom_errorbar(aes(ymin = PGLP_mean - PGLP_SE, ymax = PGLP_mean + PGLP_SE, colour = tavg), width = 0)
  p <- p + geom_errorbarh(aes(xmin = PRK_mean - PRK_SE, xmax = PRK_mean + PRK_SE, colour = tavg), width = 0)
  p <- p + expand_limits(y=0)
  p <- p + ylab('Mean PGLP abundance (mg/m2)') + xlab('Mean PRK abundance (mg/m2)')
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.position="top",
                 axis.line = element_line(colour = "black"),
                 text = element_text(size = 20))
  p
  dev.off()


  y = as.numeric(coef(model)[1]) + as.numeric(coef(model)[2]) * min(data_means$PRK_mean)
  x = as.numeric(coef(model)[1]) + as.numeric(coef(model)[2]) * max(data_means$PRK_mean)
  
  rise = x - y 
  
  run = max(data_means$PRK_mean) - min(data_means$PRK_mean)
  
  rise/run
  
  
  # rbact per PGK colored by T. 
  
  model <- summary(lm(rbact_mean ~ PGK_mean, data_means))
  
  outDir <- 'output/figures/20170217/tiff'
  depvar <-  'rbact_mean'
  indepvar <- 'PGK_mean'
  
  
  tiff(paste(outDir, '/', depvar, '_vs_', indepvar, '_R2-', round(model$r.squared,3), '_pval-', round(model$coefficients[,4][2],2),'.tiff', sep = ""), height = 6, width = 6, units = 'in', res = 600)
  p <- ggplot(data_means, aes(y = rbact_mean, x = PGK_mean)) 
  p <- p + geom_point(size = 4, aes(colour = tavg)) 
  p <- p + scale_colour_gradientn(colours=c('cyan', 'red'), name = expression(paste("Mean annual temperature (",degree,"C)", sep = "")))
  p <- p + geom_smooth(method = 'lm', size=1, colour = 'blue', se=FALSE)
  p <- p + geom_errorbar(aes(ymin = rbact_mean - rbact_SE, ymax = rbact_mean + rbact_SE, colour = tavg), width = 0)
  p <- p + geom_errorbarh(aes(xmin = PGK_mean - PGK_SE, xmax = PGK_mean + PGK_SE, colour = tavg), width = 0)
  p <- p + expand_limits(y=0)
  p <- p + ylab('Mean rbact abundance (mg/m2)') + xlab('Mean PGK abundance (mg/m2)')
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.position="top",
                 axis.line = element_line(colour = "black"),
                 text = element_text(size = 20))
  p
  dev.off()
  

  
  # rbact per Calvin ex rbact colored by T. 
  
  model <- summary(lm(rbact_mean ~ Calv_ex_rbact_mean, data_means))
  
  outDir <- 'output/figures/20170217/tiff'
  depvar <-  'rbact_mean'
  indepvar <- 'Calv_ex_rbact_mean'
  
  
  tiff(paste(outDir, '/', depvar, '_vs_', indepvar, '_R2-', round(model$r.squared,3), '_pval-', round(model$coefficients[,4][2],2),'.tiff', sep = ""), height = 6, width = 6, units = 'in', res = 600)
  p <- ggplot(data_means, aes(y = rbact_mean, x = Calv_ex_rbact_mean)) 
  p <- p + geom_smooth(method = 'lm', size=1, colour = 'blue', se=FALSE)
  p <- p + geom_point(size = 4, aes(colour = tavg)) 
  p <- p + scale_colour_gradientn(colours=c('cyan', 'red'), name = expression(paste("Mean annual temperature (",degree,"C)", sep = "")))
  p <- p + geom_errorbar(aes(ymin = rbact_mean - rbact_SE, ymax = rbact_mean + rbact_SE, colour = tavg), width = 0)
  p <- p + geom_errorbarh(aes(xmin = Calv_ex_rbact_mean - Calv_ex_rbact_SE, xmax = Calv_ex_rbact_mean + Calv_ex_rbact_SE, colour = tavg), width = 0)
  p <- p + expand_limits(y=0)
  p <- p + ylab('Mean rbact abundance (mg/m2)') + xlab('Mean Calvin ex rbact abundance (mg/m2)')
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.position="top",
                 axis.line = element_line(colour = "black"),
                 text = element_text(size = 20))
  p
  dev.off()
  
  
  
# proportion of total protein for Calvin cycle and Photosystems on same plot vs % canopy openness
 
  source('scripts/prep_data.R') 
  
  
  data_means <- data %>% group_by(ID) %>% summarise(Photosystems_mean = mean(Photosystems, na.rm=TRUE),
                                                   Photosystems_SE = SE(Photosystems),
                                                   Calvin_cycle_mean = mean(Calvin_cycle, na.rm=TRUE),
                                                   Calvin_cycle_SE = SE(Calvin_cycle)) %>%
    full_join(data, by = 'ID')

  tiff(paste(outDir,'/calv_photosystems_vs_gap.tiff', sep = ""), height = 6, width = 6, units = 'in', res = 600)

  p <- ggplot(data_means, aes(x=gap_mean)) 
  p <- p + geom_point(size = 4, colour = 'forestgreen', shape = 17, aes(x = gap_mean, y = Photosystems_mean))
  p <- p + geom_errorbar(colour = 'forestgreen', aes(x = gap_mean, ymin = Photosystems_mean - Photosystems_SE, ymax = Photosystems_mean + Photosystems_SE))
  p <- p + geom_errorbarh(colour = 'forestgreen', aes(x = gap_mean, y = Photosystems_mean, xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE))
  
  p <- p + geom_point(size = 4, colour = 'blue', shape = 16, aes(x = gap_mean, y = Calvin_cycle_mean))
  p <- p + geom_errorbar(colour = 'blue', aes(x = gap_mean, ymin = Calvin_cycle_mean - Calvin_cycle_SE, ymax = Calvin_cycle_mean + Calvin_cycle_SE))
  p <- p + geom_errorbarh(colour = 'blue', aes(x = gap_mean, y = Calvin_cycle_mean, xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE))
  
  p <- p + geom_smooth(aes(x = gap_mean, y = Photosystems_mean), method = 'lm', size=1, colour = 'forestgreen', se=FALSE)
  #p <- p + geom_smooth(aes(x = gap_mean, y = Calvin_cycle_mean), method = 'lm', size=1, colour = 'blue', se=FALSE)
  p <- p + expand_limits(y=0)
  p <- p + ylab('Mean proportion of total protein') + xlab('Canopy openness (%)')
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.position="top",
                 axis.line = element_line(colour = "black"),
                 text = element_text(size = 20))
  p
  dev.off()
  
#-DP14 Eucs and SIEF Eucs across E Australia--same as UTS application fig, just with SIEF points added.
#-Calvin/Light Rxns vs % canopy openness (to replace vs mean annual irradiance)
#-RCA/PGK vs MAT or another version of T that is a stronger predictor (MAT is R2 0.28)
#-PGLP vs PRK colored by T. Could you please do a version without the CI? Also, units are mg/m^2, is that right?

  
  
  
  

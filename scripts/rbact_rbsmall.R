# rubisco activase

rb <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.j02030.1.p', 'eucgr.j01234.1.p',
                                                             'eucgr.j01502.1.p',
                                                             'eucgr.j01502.2.p',
                                                             'eucgr.c00150.1.p',
                                                             'eucgr.b03013.1.p',
                                                             'eucgr.k02223.1.p'),]
rb <- arrange(rb, Protein)

names <- rb$Protein

rb <- as.data.frame(t(rb[,2:314]))

names(rb) <- names
rb$sample <- rownames(rb)
rownames(rb) <- NULL

rb <- merge(rb, total_protein_D14)

rb[,2:8] <- lapply(rb[,2:8],function(x) as.numeric(as.character(x)))

#rb[,2:8] <- rb[,2:8]/rb$total_protein
str(rb)

rb$rb_smallsub_mean <- rowMeans(rb[,c('eucgr.j01502.1.p', 'eucgr.j01502.2.p', 'eucgr.c00150.1.p', 'eucgr.b03013.1.p', 'eucgr.k02223.1.p')])

rb$rb_act_mean <- rowMeans(rb[,c('eucgr.j02030.1.p', 'eucgr.j01234.1.p')])

bla <- merge(rb, climate_locs)

plot(bla$eucgr.j01234.1.p ~ bla$tavg, ylab = 'eucgr.j01234.1.p', xlab = 'tavg (oC)')
abline(lm(bla$eucgr.j01234.1.p ~ bla$tavg))
summary(lm(bla$eucgr.j01234.1.p ~ bla$tavg))

blax <- merge(bla, replicates)

blax <- blax %>% group_by(ID) %>% summarise(rbact_mean = mean(eucgr.j01234.1.p, na.rm=TRUE),
                                            rbact_SE = SE(eucgr.j01234.1.p))

blax <- merge(blax, merge(climate_locs, replicates))

plot(blax$rbact_mean ~ blax$tavg)

rbact.lm <- summary(lm(blax$rbact_mean ~ blax$tavg))

round(rbact.lm$coefficients[,4][2],2)


h <- ggplot(blax, aes(x = tavg, y = rbact_mean)) + geom_point() 
h <- h + geom_errorbar(aes(ymin = rbact_mean - rbact_SE, ymax = rbact_mean + rbact_SE), width = 0.2, alpha = 0.8)
h <- h + geom_smooth(method = 'lm')
h <- h + expand_limits(y=0)
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean), label = paste('R2 =', round(rbact.lm$r.squared,3), sep = " "))
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean)*0.9, label = paste('pval =', round(rbact.lm$coefficients[,4][2],2), sep = " "))
h <- h + ylab('mean eucgr.j01234.1.p')
h <- h + theme_bw()
h <- h + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
h



plot(bla$eucgr.j02030.1.p ~ bla$tavg)
abline(lm(bla$eucgr.j02030.1.p ~ bla$tavg))
summary(lm(bla$eucgr.j02030.1.p ~ bla$tavg))



blax <- merge(bla, replicates)

blax <- blax %>% group_by(ID) %>% summarise(rbact_mean = mean(eucgr.j02030.1.p, na.rm=TRUE),
                                            rbact_SE = SE(eucgr.j02030.1.p))

blax <- merge(blax, merge(climate_locs, replicates))

plot(blax$rbact_mean ~ blax$tavg)

rbact.lm <- summary(lm(blax$rbact_mean ~ blax$tavg))

round(rbact.lm$coefficients[,4][2],2)


h <- ggplot(blax, aes(x = tavg, y = rbact_mean)) + geom_point() 
h <- h + geom_errorbar(aes(ymin = rbact_mean - rbact_SE, ymax = rbact_mean + rbact_SE), width = 0.2, alpha = 0.8)
h <- h + geom_smooth(method = 'lm')
h <- h + expand_limits(y=0)
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean), label = paste('R2 =', round(rbact.lm$r.squared,3), sep = " "))
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean)*0.9, label = paste('pval =', round(rbact.lm$coefficients[,4][2],2), sep = " "))
h <- h + ylab('mean eucgr.j02030.1.p')
h <- h + theme_bw()
h <- h + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
h






blax <- merge(bla, replicates)

blax <- blax %>% group_by(ID) %>% summarise(rb_small_mean = mean(eucgr.k02223.1.p, na.rm=TRUE),
                                            rb_small_SE = SE(eucgr.k02223.1.p))

blax <- merge(blax, merge(climate_locs, replicates))

plot(blax$rb_small_mean ~ blax$tavg)

rb_small.lm <- summary(lm(blax$rb_small_mean ~ blax$tavg))

round(rb_small.lm$coefficients[,4][2],2)


h <- ggplot(blax, aes(x = tavg, y = rb_small_mean)) + geom_point() 
h <- h + geom_errorbar(aes(ymin = rb_small_mean - rb_small_SE, ymax = rb_small_mean + rb_small_SE), width = 0.2, alpha = 0.8)
h <- h + geom_smooth(method = 'lm')
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rb_small_mean), label = paste('R2 =', round(rb_small.lm$r.squared,3), sep = " "))
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rb_small_mean)*0.9, label = paste('pval =', round(rb_small.lm$coefficients[,4][2],2), sep = " "))
h <- h + ylab('mean eucgr.k02223.1.p')
h <- h + theme_bw()
h <- h + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
h

plot(blax$rb_small_SE ~ blax$tavg)

plot(bla$rb_smallsub_mean ~ bla$tavg)

plot(bla$rb_act_mean/bla$rb_smallsub_mean ~ bla$tavg)

plot(eucgr.j01502.1.p ~ tavg, bla[bla$eucgr.j01502.1.p > max(bla$eucgr.j01502.1.p) * 0.01,])
plot(eucgr.j01502.2.p ~ tavg, bla[bla$eucgr.j01502.2.p > max(bla$eucgr.j01502.2.p) * 0.01,])
plot(eucgr.c00150.1.p ~ tavg, bla[bla$eucgr.c00150.1.p > max(bla$eucgr.c00150.1.p) * 0.01,])
plot(eucgr.b03013.1.p ~ tavg, bla[bla$eucgr.b03013.1.p > max(bla$eucgr.b03013.1.p) * 0.01,])
plot(eucgr.k02223.1.p ~ tavg, bla[bla$eucgr.k02223.1.p > max(bla$eucgr.k02223.1.p) * 0.01,])

 


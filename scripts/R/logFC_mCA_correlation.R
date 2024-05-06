

logFC_mCA_corr = function(df,plotType){
  logFC_mCA.nuclear.bound = df %>% 
    filter(Binding == 'Bound') %>%
    filter(DEG!='Unchanged') %>% 
    filter(Meth_CA=='CA') %>%
    #mutate(MECP2 = ifelse((distanceToTSS <=3000 & distanceToTSS >=-3000),'Promoter','Distal')) %>% 
    dplyr::select(logFC,CA_meth) %>% 
    arrange(desc(CA_meth)) %>%
    distinct()
  logFC_mCA.nuclear.all = df %>% 
    filter(DEG!='Unchanged') %>% 
    filter(Meth_CA=='CA') %>%
    #mutate(MECP2 = ifelse((distanceToTSS <=3000 & distanceToTSS >=-3000),'Promoter','Distal')) %>% 
    dplyr::select(logFC,CA_meth) %>% 
    arrange(desc(CA_meth)) %>%
    distinct()
  
  
  df.bin =list()
  summarize_logFC_mCA = function(x,binSize = 19,stepSize =4){
    for (i in seq(from=1, to=nrow(x), by=stepSize)){
      print(i)
      start = i
      stop  = i+binSize
      if(stop <= nrow(x)){
        df = x[start:stop,]
        df.bin[[i]] = df  %>% summarize(logFC.mean=mean(logFC),logFC.sd = sd(logFC),CA_meth.mean=mean(CA_meth))
      }
      else{
        return(df.bin)
      }
    }
  }
  x = summarize_logFC_mCA(logFC_mCA.nuclear.bound, binSize = 39,stepSize =4)
  y = summarize_logFC_mCA(logFC_mCA.nuclear.all, binSize = 79,stepSize =8)
  MECP2_bound = do.call(rbind,x) %>% mutate(MECP2 ='MECP2 Bound')
  MECP2_all   = do.call(rbind,y) %>% mutate(MECP2 ='All')
  
  plt = list()
  plt.boxplot = list()
  plt[[1]] = rbind(MECP2_bound,MECP2_all)  %>% 
    ggplot(aes(x=CA_meth.mean,y=logFC.mean,color=MECP2))+
    geom_point(size=0.5)+
    geom_line()+
    geom_ribbon(aes(x=CA_meth.mean,y=logFC.mean,
                    ymin=logFC.mean-logFC.sd,
                    ymax=logFC.mean+logFC.sd,fill=MECP2),
                alpha=0.1,color='NA') +
    scale_fill_manual(name = "DEGs (MECP2 KO vs WT)",values = c('black','red'),aesthetics = c("colour", "fill"))+
    stat_cor()+
    theme_test(20)+
    xlab('mCA/CA')+
    ylab('Mean log2 Fold Change')+
    guides(color = guide_legend(override.aes = aes(label = "")))
  
  
  plt[[1]]
  
  plt.boxplot[[1]] = rbind(MECP2_bound,MECP2_all) %>% 
    ggplot(aes(x=MECP2,y=logFC.mean,color=MECP2))+
    geom_boxplot()+stat_compare_means(comparisons = list(c('All','MECP2 Bound'))) + 
    scale_color_manual(values = c('black','red'))+
    theme_test(20) + 
    xlab('DEGs (MECP2 KO vs WT)')
  
  
  #########################################################################################################################################
  MECP2_bound.pp = df %>% 
    filter(Binding == 'Bound') %>%
    filter(DEG!='Unchanged') %>% 
    filter(Meth_CA=='CA') %>% 
    mutate(MECP2 = ifelse(distanceToTSS <3000 & distanceToTSS >-3000,'Promoter proximal','Distal')) %>% 
    filter(MECP2 == 'Promoter proximal') %>% 
    dplyr::select(logFC,CA_meth) %>% 
    arrange(desc(CA_meth)) %>% 
    distinct()
  
  MECP2_bound.distal = df %>% 
    filter(Binding == 'Bound') %>%
    filter(DEG!='Unchanged') %>% 
    filter(Meth_CA=='CA') %>% 
    mutate(MECP2 = ifelse(distanceToTSS <3000 & distanceToTSS >-3000,'Promoter proximal','Distal')) %>% 
    filter(MECP2 == 'Distal') %>% 
    dplyr::select(logFC,CA_meth) %>% 
    arrange(desc(CA_meth)) %>% 
    distinct()
  
  x = summarize_logFC_mCA(MECP2_bound.pp, binSize = 19,stepSize =2)
  z = summarize_logFC_mCA(MECP2_bound.distal, binSize = 39,stepSize =4)
  MECP2_bound.pp = do.call(rbind,x) %>% mutate(MECP2 ='MECP2 Bound (Promoter proximal)')
  MECP2_bound.distal  = do.call(rbind,z) %>% mutate(MECP2 ='MECP2 Bound (Distal)')
  
  plt[[2]] = rbind(MECP2_bound.pp,MECP2_bound.distal,MECP2_all)  %>% 
    ggplot(aes(x=CA_meth.mean,y=logFC.mean,color=MECP2))+
    geom_point(size=0.5)+
    geom_line()+
    geom_ribbon(aes(x=CA_meth.mean,y=logFC.mean,
                    ymin=logFC.mean-logFC.sd,
                    ymax=logFC.mean+logFC.sd,fill=MECP2),
                alpha=0.1,color='NA') +
    scale_fill_manual(name = "DEGs (MECP2 KO vs WT)",values = c('black','red','blue'),aesthetics = c("colour", "fill"))+
    stat_cor()+
    theme_test(20)+
    xlab('mCA/CA')+
    ylab('Mean log2 Fold Change')+
    guides(color = guide_legend(override.aes = aes(label = "")))
  plt[[2]]
  
  plt.boxplot[[2]] = rbind(MECP2_bound.pp,MECP2_bound.distal,MECP2_all)  %>% 
    ggplot(aes(x=MECP2,y=logFC.mean,color=MECP2))+
    geom_boxplot()+stat_compare_means(comparisons = list(c('All','MECP2 Bound (Distal)'),
                                                         c('All','MECP2 Bound (Promoter proximal)'),
                                                         c('MECP2 Bound (Distal)','MECP2 Bound (Promoter proximal)'))) + 
    scale_color_manual(values = c('black','red','blue'))+
    theme_test(20) +
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))+
    xlab('DEGs (MECP2 KO vs WT)')
  plt.boxplot[[2]] 
  
  #########################################################################################################################################
  logFC_mCA.nuclear.bound_with_count = df %>% 
    filter(Binding == 'Bound') %>%
    filter(DEG!='Unchanged') %>% 
    filter(Meth_CA=='CA') %>%
    dplyr::select(logFC,CA_meth) %>% 
    group_by(logFC,CA_meth) %>% 
    summarize(Count=dplyr::n()) %>% 
    mutate(MECP2 = case_when(Count <7~'1-6 MBHs',
                             Count >=7~'>=7 MBHs')) 
  
  logFC_mCA.nuclear.bound_with_count.list = split(logFC_mCA.nuclear.bound_with_count, f =logFC_mCA.nuclear.bound_with_count$MECP2)
  
  logFC_mCA.nuclear.bound_with_count.list = lapply(logFC_mCA.nuclear.bound_with_count.list, function(x) x %>% arrange(desc(CA_meth)))
  
  tmp = lapply(logFC_mCA.nuclear.bound_with_count.list, function(x) summarize_logFC_mCA(x %>% as.data.frame(),binSize = 19,stepSize = 4))
  tmp1 = lapply(tmp, function(x) do.call(rbind,x))
  tmp2 = plyr::ldply(tmp1,data.frame)
  colnames(tmp2)[1] = 'MECP2'
  logFC_mCA.nuclear.bound_with_count.avg = tmp2[,c(2,3,4,1)]
  
  
  plt[[3]] = rbind(logFC_mCA.nuclear.bound_with_count.avg,MECP2_all)  %>% 
    filter(MECP2 !='1-6 MBHs') %>% 
    ggplot(aes(x=CA_meth.mean,y=logFC.mean,color=MECP2))+
    
    geom_point(size=0.5)+
    geom_line()+
    geom_ribbon(aes(x=CA_meth.mean,y=logFC.mean,
                    ymin=logFC.mean-logFC.sd,
                    ymax=logFC.mean+logFC.sd,fill=MECP2),
                alpha=0.1,color='NA') +
    scale_fill_manual(name = "DEGs (MECP2 KO vs WT)",values = c('red','black'),aesthetics = c("colour", "fill"))+
    stat_cor()+
    theme_test(20)+
    xlab('mCA/CA')+
    ylab('Mean log2 Fold Change')+
    #ylim(-0.9,0.9)+
    #scale_y_continuous(limits = c(-0.6, 0.7), breaks = seq(-0.6, 0.7, by = 0.3))+
    guides(color = guide_legend(override.aes = aes(label = "")))
  plt[[3]]
  
  plt.boxplot[[3]] = rbind(logFC_mCA.nuclear.bound_with_count.avg,MECP2_all)  %>% 
    mutate(MECP2 = factor(MECP2, levels=c('All','1-6 MBHs','>=7 MBHs'))) %>% 
    filter(MECP2 !='1-6 MBHs') %>% 
    ggplot(aes(x=MECP2,y=logFC.mean,color=MECP2))+
    geom_boxplot()+stat_compare_means(comparisons = list(c('All','>=7 MBHs')))+
    #c('All','1-6 MBHs'),
    
    #c('1-6 MBHs','>=7 MBHs'))) + 
    scale_color_manual(values = c('black','red'))+
    theme_test(20) + 
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))+
    xlab('DEGs (MECP2 KO vs WT)')
  plt.boxplot[[3]]
  if(plotType == 'linePlot'){
    return(plt)
  }
  if(plotType == 'boxplot'){
    return(plt.boxplot)
  }
  
}
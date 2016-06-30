
forest.plot <- function(dat_lis,mode='single',x_lim=c(0,10),use_log=NA,phen=NA,module=NA,p_thresh=NA,point_width_scale=4,line_width=0.1,ci_bar_height=0.35){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 needed for this function to work. Please install it.",
    call. = FALSE)
  }
  #library(ggplot2)
  cat('\tWARNING : currently does not work for a single module\n\n')

#  colvec=c('red','darkblue','darkorange','magenta','orange')
   colvec=colmix
   cat('\tphenotypes detected :',paste(names(dat_lis),collapse=', '),'\n')
dat_lis=dat_lis[rev(names(dat_lis))]
options(warn=-1)
  if(!is.na(phen)){
     cat('\tdataset selected:',paste(phen,collapse=', '),'\n')
    dat_lis=dat_lis[names(dat_lis)%in%phen]
  }
    if(!is.na(module)){
       cat('\tmodules selected:',paste(module,collapse=', '),'\n')
  }
  plot_dat=list()
  for(ilis in 1:length(dat_lis)){
    if(!is.na(module)){
      dat_lis[[names(dat_lis)[ilis]]]=dat_lis[[names(dat_lis)[ilis]]][rownames(dat_lis[[names(dat_lis)[ilis]]])%in%module,]
   }

    plot_dat[[names(dat_lis)[ilis]]]=as.data.frame(make.numeric(dat_lis[[names(dat_lis)[ilis]]][,c('module size','FET p.value','OR','[95% OR CI] inf','OR [95% OR CI] sup')]))
      colnames(plot_dat[[names(dat_lis)[ilis]]])=c('n.genes','FETp','fetOR','lowerCI','upperCI')

    plot_dat[[names(dat_lis)[ilis]]]$module=as.factor(paste(rownames(plot_dat[[names(dat_lis)[ilis]]]),names(dat_lis)[ilis],sep='_'))
    plot_dat[[names(dat_lis)[ilis]]]$color=colvec[ilis]#c(rep(colvec[ilis],nrow(plot_dat[[names(dat_lis)[ilis]]])-1),'black')
#    plot_dat[[names(dat_lis)[ilis]]]$color=c(rep(colvec[ilis],nrow(plot_dat[[names(dat_lis)[ilis]]])-1),'black')
  plot_dat[[names(dat_lis)[ilis]]]$point_width=(plot_dat[[names(dat_lis)[ilis]]]$fetOR*3 )+point_width_scale
  }


   if(!is.na(p_thresh)){
    dummy=list()
    for(ilis in 1:length(plot_dat)){
      tempr=plot_dat[[names(dat_lis)[ilis]]][plot_dat[[names(dat_lis)[ilis]]]$FETp<p_thresh,]
      if(nrow(tempr)>0){
        dummy[[names(plot_dat)[ilis]]]=plot_dat[[names(plot_dat)[ilis]]][plot_dat[[names(plot_dat)[ilis]]]$FETp<p_thresh,]
      }
    }
    plot_dat=dummy
   }
options(warn=0)

#  if(!is.na(use_log)){
#    plot_dat$fetOR=log(as.numeric(plot_dat$fetOR),base=use_log)
#    plot_dat$lowerCI=log(as.numeric(plot_dat$lowerCI),base=use_log)
#    plot_dat$upperCI=log(as.numeric(plot_dat$upperCI),base=use_log)

#  }

####### since gplots insists on reversing all the phenotypes, pre-empt it to keep original rownames
  for(ilis in 1:length(plot_dat)){
      plot_dat[[names(plot_dat)[ilis]]]=plot_dat[[names(plot_dat)[ilis]]][rev(rownames(plot_dat[[names(plot_dat)[ilis]]])),]
   }


gplots_dat=list()
  if(mode=='single'){
    for(idat in 1:length(plot_dat)){
    myplot=ggplot() +                                                     # initiate the plot space, saved to myplot variable for later display
      #  - subsequent additions / modifications of the plot can be done by manipulating this space
      #  geom_point(data = df, aes(y = Module, x = OR),colour = 'red', size = 3) +
      #  geom_errorbar(data=pdat,aes(x=Module,y=OR,ymin=lower.95..CI,ymax=upper.95..CI),colour="grey60",width=0.2)
      geom_vline(xintercept=1, linetype="dashed")+                         # add vertical line
      #    geom_vline(xintercept=(1), linetype="dashed")+                         # add vertical line
      geom_errorbarh(data=plot_dat[[idat]], aes(y=module,x=fetOR, xmin=lowerCI, xmax=upperCI), colour="black", height=ci_bar_height,size=2)+  # add HORISONTAL error bars, geom_errorbar used for vertical..
      geom_point(data = plot_dat[[idat]], aes(y = module, x = fetOR), shape=15, cex=plot_dat[[idat]]$point_width, colour=plot_dat[[idat]]$color)+    # add plotting points
      #    xlim(0, 5)                                                         # removes datapoints outside the range == scale_x_continuous(limits = c(-5000, 5000))
      #
      #    coord_cartesian(xlim = c(0, 60)) +                                  # moves the window only
      coord_cartesian(xlim = x_lim) +    #        coord_cartesian(xlim = c(-0.5, 1))                        # moves the window only
      scale_x_continuous(breaks=c(0,1:9,seq(10,2000,by=10))) +
      #    scale_x_continuous(breaks=c(0,1,5,seq(10,2000,by=10))) +
      #  scale_size_area() +
      xlab("FET odds ratio and CI") +
      ylab("module") 
    ##  plot line styles      ||  http://www.cookbook-r.com/Graphs/Shapes_and_line_types/  
    #  ggtitle(names(table(dat$serie.name)[idat]))

    ## essential if using the ugly mess that is ggplot2, for more options : http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/ 
    #  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #
      gplots_dat[[names(plot_dat)[idat]]]=(myplot + theme_bw())
    }
  }

  if(mode=='combined'){
    plot_merge=plot_dat[[1]]
    if(length(plot_dat)>1){
        for(ilis in 2:length(plot_dat)){
          plot_merge=rbind(plot_merge,plot_dat[[ilis]])
        }}
    myplot=ggplot() +                                                     # initiate the plot space, saved to myplot variable for later display
      #  - subsequent additions / modifications of the plot can be done by manipulating this space
      #  geom_point(data = df, aes(y = Module, x = OR),colour = 'red', size = 3) +
      #  geom_errorbar(data=pdat,aes(x=Module,y=OR,ymin=lower.95..CI,ymax=upper.95..CI),colour="grey60",width=0.2)
      geom_vline(xintercept=1, linetype="dashed")+                         # add vertical line
      #    geom_vline(xintercept=(1), linetype="dashed")+                         # add vertical line
      geom_errorbarh(data=plot_merge, aes(y=module,x=fetOR, xmin=lowerCI, xmax=upperCI), colour="black", height=ci_bar_height,size=2)+  # add HORISONTAL error bars, geom_errorbar used for vertical..
      geom_point(data=plot_merge, aes(y = module, x = fetOR), shape=15, cex=plot_merge$point_width, colour=plot_merge$color)+    # add plotting points
      #    xlim(0, 5)                                                         # removes datapoints outside the range == scale_x_continuous(limits = c(-5000, 5000))
      #
      #    coord_cartesian(xlim = c(0, 60)) +                                  # moves the window only
      coord_cartesian(xlim = x_lim) +    #        coord_cartesian(xlim = c(-0.5, 1))                        # moves the window only
      scale_x_continuous(breaks=c(0,1:9,seq(10,2000,by=10))) +
      #    scale_x_continuous(breaks=c(0,1,5,seq(10,2000,by=10))) +
      #  scale_size_area() +
      xlab("FET odds ratio and CI") +
      ylab("module") 
    ##  plot line styles      ||  http://www.cookbook-r.com/Graphs/Shapes_and_line_types/  
    #  ggtitle(names(table(dat$serie.name)[idat]))



    ## essential if using the ugly mess that is ggplot2, for more options : http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/ 
    #  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #
      gplots_dat=(myplot + theme_bw())
  }
  return(gplots_dat)
}



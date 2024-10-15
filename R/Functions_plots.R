#Function that calculates the yearly cumulative incidence (per 100,000 patients)
#output : output of the model
#n_sim: number of simulations
yearly_cum_inc = function(output,n_sim) {
  
  yearly_inc_baseline = c()
  
  for(i in 1:n_sim) {
    
    yearly_inc_baseline[i] = 100000*(sum(output[[i]]$incidence)/(output[[i]]$count_patients - sum(output[[i]]$cont_new_patients)))
    
  }
  
  mean_inc = mean(yearly_inc_baseline)
  inf_inc = quantile(yearly_inc_baseline,probs = c(0.025))
  sup_inc = quantile(yearly_inc_baseline,probs = c(0.975))
  
  return(cum_inc = c(round(mean_inc,2),round(inf_inc,2),round(sup_inc,2)))
  
}

#Function that plots the daily incidence over a given period 
#output : output of the model
#n_sim: number of simulations
#time:tmax for plotting the daily incidence
#ylim = min.and max.of the y-axis
plot_incidence = function(output, n_sim, time, time_step, ylim = c(0,1)) {
  
  time=time*(1/time_step)
  oord_lim =ylim/100
  
  seq_days = seq(0,time,288)
  seq_days[1]=1
  
  mat_inc_sim = matrix(nrow = n_sim,ncol = length(seq_days[-1]))
  
  for(s in 1:n_sim) {
    
    inc_rate = c()
    
    for(i in 1:(length(seq_days)-1)) {
      
      numerator = sum(output[[s]]$incidence[1,][seq_days[i]:seq_days[i+1]])
      denominator  = output[[s]]$s_patients[seq_days[i]]+(sum(output[[s]]$new_patients[1,][seq_days[i]:seq_days[i+1]])-sum(output[[s]]$cont_new_patients[1,][seq_days[i]:seq_days[i+1]]))
      
      inc_rate[i] = numerator/denominator  
      
    }
    
    mat_inc_sim[s,] = inc_rate
    
  }
  
  #Prediction interval
  df_inc_rate = matrix(ncol = 3,nrow = length(seq_days[-1]))
  
  for(i in 1:length(seq_days[-1])) {
    
    df_inc_rate [i,1] = mean(mat_inc_sim[,i])
    df_inc_rate[i,2] = quantile(mat_inc_sim[,i],probs = c(0.025))
    df_inc_rate[i,3] = quantile(mat_inc_sim[,i],probs = c(0.975))
    
  }
  
  df_inc_rate = data.frame(df_inc_rate)
  colnames(df_inc_rate)= c("mean","inf_pred","sup_pred")
  
  df_melt = melt(mat_inc_sim)
  colnames(df_melt) = c("sim","day","inc_rate")
  
  plot=ggplot(data =df_melt,aes(x=day,y=inc_rate,group=sim))+
    geom_point(aes(x=day,y=inc_rate,group=day),stat = "summary", fun = "mean",fill="#AA4499",size=2,shape=21)+
    geom_ribbon(data=df_inc_rate,aes(ymin = inf_pred, ymax = sup_pred,group=1,x=1:nrow(df_inc_rate),y=mean), alpha = 0.2,fill="#AA4499")+
    #scale_fill_manual(values="green")+
    ylab("Daily incidence rate (%)")+
    scale_y_continuous(labels = function(a) paste0(a*100, "%"),limits =oord_lim)+
    xlab("Time")+
    scale_x_continuous(labels = function(a) paste0(a, " days"))+
    #coord_trans( x = "day")+
    theme_classic()+
    theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"))
  
  return(list(plot,df_inc_rate,mat_inc_sim))
  
}

#Function that plots the yearly cumulative incidence in each ward
#output : output of the model
#n_sim: number of simulations
#n_wards: number of wards
#rm_0 : logical, should the wards with an incidence of 0 be removed 
#type : "mean" or "median", should the central plotted value of each ward be the mean or the median
plot_inc_wards = function(output,n_sim,n_wards,rm_0=T,type="mean") {
  
  mat_inc=matrix(ncol = n_sim,nrow = n_wards)
  
  mat_cases=matrix(ncol = n_sim,nrow = n_wards)
  
  for(w in 1:n_wards) {
    
    yearly_inc_ward = c()
    
    for(i in 1:n_sim) {
      
      mat_cases[,i] = output[[i]]$ward_event;
      
      # yearly_inc_ward[i] = (output[[i]]$ward_event[w]/(output[[i]]$count_patient_ward[w]))
      yearly_inc_ward[i] = output[[i]]$ward_event[w]/output[[i]]$count_patient_ward[w]
      
      if(is.nan(yearly_inc_ward[i]) | is.infinite(yearly_inc_ward[i])) {yearly_inc_ward[i]=0}
      
    }
    
    mat_inc[w,] = yearly_inc_ward
    
  }
  
  df_inc = matrix(ncol = 6, nrow = n_wards)
  
  for (w in 1:n_wards) {
    
    if(type=="median"){df_inc[w,1] = median(mat_inc[w,])}
    if(type=="mean"){df_inc[w,1] = mean(mat_inc[w,])}
    
    df_inc[w,2] = quantile(mat_inc[w,],prob=0.025)
    df_inc[w,3] = quantile(mat_inc[w,],prob=0.975)
    
    if(df_inc[w,2] == 0 & df_inc[w,3] ==0 & df_inc[w,1] != 0) {df_inc[w,3] =df_inc[w,1] }
    
    if(type=="median"){df_inc[w,4] = median(mat_cases[w,])}
    if(type=="mean"){df_inc[w,4] = mean(mat_cases[w,])}
    
    df_inc[w,5] = quantile(mat_cases[w,],prob=0.025)
    df_inc[w,6] = quantile(mat_cases[w,],prob=0.975)
  }
  
  df_inc = data.frame(ward=departments,df_inc)
  df_inc = data.frame(ward=as.character(1:28),df_inc)
  colnames(df_inc) = c("ward_nb","ward","mean","pred_inf","pred_sup","mean_cases","pred_inf_cases","pred_sup_cases")
  
  df_test = melt(mat_inc) ;
  df_test = data.frame(df_test,mean_cases = rep(df_inc$mean_cases,n_sim),as.character(rep(df_hosp_ward$department,n_sim)),
                       ward_name=rep(departments,n_sim))
  colnames(df_test) = c("ward","sim","y_inc","mean_cases","dpt","ward_name")
  
  print(df_test)
  print(df_inc)
  
  plot=ggplot(data =df_test,aes(y=reorder(ward_name, y_inc, FUN=mean),x=y_inc,group=ward_name))+
    geom_jitter(shape=23,fill="grey",alpha=0.3)+
    geom_errorbar(data=df_inc,aes(xmin=pred_inf,xmax=pred_sup,color=mean,y=reorder(ward,mean),x=mean),width=.3,size=0.75,color="black",inherit.aes=FALSE)+
    guides(fill=guide_legend(title="Average number \n of cases"))+
    ggnewscale::new_scale_fill()+
    geom_point(stat = "summary", fun = "mean",aes(y=reorder(ward_name, y_inc, FUN=mean),x=y_inc,group=ward,fill=mean_cases),size=5,shape=23)+
    scale_x_continuous(labels = function(x) paste0(x*100, "%"))+
    xlab("Yearly cumulative incidence(95% PI)")+
    ylab("Ward")+
    scale_fill_gradient2(midpoint = mean(df_test$mean_cases), high = "#0B4151" , mid = "#2A9882",
                         low= "#D4F3A3", space = "Lab" )+
    theme_bw()+
    guides(fill=guide_legend(title="Average number \n of cases"))+
    theme(axis.text.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 22, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 24, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 24, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          legend.title=element_text(size=24),
          legend.text = element_text(size=22) )
  
  return(list(plot,df_inc, mat_inc))
  
}

#Function that plots the attributable portion to new cases for each type of device
#output : output of the model
#n_sim: number of simulations
#n_eq: number of devices
#plot_log: logical, should the calculated portion be transformed into log 
plot_att_eq = function(output,n_sim,n_eq,plot_log=F) {
  
  mat_eq_inf=matrix(ncol = n_sim,nrow = n_eq)
  
  mat_mean_inf=matrix(ncol = n_sim,nrow = n_eq)
  
  index_no_inf = c()
  
  for(s in 1:n_sim) {
    
    #mat_eq_prop[,s] = output_baseline[[s]]$cont_mat/output_baseline[[s]]$used_mat
    if(sum(output[[s]]$inf_eq)!= 0) {mat_eq_inf[,s] = output[[s]]$inf_eq/sum(output[[s]]$inf_eq)}
    
    if(sum(output[[s]]$inf_eq) == 0) {index_no_inf = c(index_no_inf,s)}
    
    mat_mean_inf[,s] = output[[s]]$inf_eq;

  }
  
  mat_eq_inf[which(is.nan(mat_eq_inf))]=0
  
  if(is.null(index_no_inf)==F) {mat_eq_inf = mat_eq_inf[,-index_no_inf]}
  
  print(mat_eq_inf)
  
  df_inf_cont = matrix(ncol = 6, nrow = n_eq)
  
  for (w in 1:n_eq) {
    
    df_inf_cont [w,1] = mean(mat_eq_inf[w,])
    df_inf_cont [w,2] = quantile(mat_eq_inf[w,],prob=0.025)
    df_inf_cont [w,3] = quantile(mat_eq_inf[w,],prob=0.975)
    df_inf_cont [w,4] = mean(mat_mean_inf[w,])
    df_inf_cont [w,5] = quantile(mat_mean_inf[w,],prob=0.025)
    df_inf_cont [w,6] = quantile(mat_mean_inf[w,],prob=0.975)
  }
  
  if(plot_log == T) {
    df_inf_cont[,1:3] = log(df_inf_cont[,1:3])
    df_inf_cont[is.infinite(df_inf_cont)] = log(0.0001)
  }
  
  df_inf_cont = data.frame(eq=supplies,df_inf_cont)
  colnames(df_inf_cont) = c("eq","mean","pred_inf","pred_sup","mean_inf","pred_inf_inf","pred_sup_inf")
  
  if(plot_log==F) {
    
    plot=ggplot(df_inf_cont,aes(y=reorder(eq,mean),x=mean))+
      geom_errorbar(aes(xmin=pred_inf,xmax=pred_sup,color=mean),width=.3,size=1,color="black")+
      scale_x_continuous(labels = function(x) paste0(x*100, "%"))+
      geom_point(aes(fill=mean_inf),size=6,shape=21)+
      scale_fill_gradient2(midpoint = mean(df_inf_cont$mean_inf), low = "white", mid = "orange",
                           high = "red", space = "Lab" )+
      scale_x_continuous(labels = function(x) paste0(x*100, "%"))+
      xlab("Yearly attributable portion to new cases (95% PI)")+
      ylab("Medical equipment")+
      theme_bw()+
      guides(fill=guide_legend(title="Average number \n of cases"))+
      theme(axis.text.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 22, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
            axis.title.x = element_text(color = "grey20", size = 24, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 24, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            legend.title=element_text(size=22),
            legend.text = element_text(size=20) ) 
    
    
  }
  
  if(plot_log==T) {
    
    plot=ggplot(df_inf_cont,aes(y=reorder(eq,mean),x=mean))+
      geom_errorbar(aes(xmin=pred_inf,xmax=pred_sup,color=mean),width=.3,size=1,color="black")+
      geom_point(aes(fill=mean_inf),size=6,shape=21)+
      scale_fill_gradient2(midpoint = mean(df_inf_cont$mean_inf), low = "white", mid = "orange",
                           high = "red", space = "Lab" )+
      xlab("log(Yearly attributable portion to new cases) (95% PI)")+
      ylab("Medical equipment")+
      theme_classic()+
      guides(fill=guide_legend(title="Average number \n of cases"))+
      theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
            axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            legend.title=element_text(size=14),
            legend.text = element_text(size=12) ) 
    
    
  } 
  
  return(plot)
  
}

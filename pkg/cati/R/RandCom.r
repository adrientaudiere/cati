
############################################## FUNCTION #########################################"

RandCom = function(com, sp, Nind.com, sdlog, min_value_traits, max_value_traits,
                   cv_intra_sp, cv_intra_com, Int_Filter_Strength, Ext_Filter_Strength, Filter='None'){

## Filter=c('None', 'External', 'Internal', 'Both')

        Ncom <- length(com)
        Nsp <- length(sp)
        Nind <- Ncom * Nind.com

        Data = data.frame(matrix(0,ncol=4, nrow=Nind)); colnames(Data) = c("com","sp","trait1","trait2")

        # 1 - Assign a community to each individual
        ### HYP 1: all communities have the same number of individuals
        Data$com = as.factor(rep(com, rep(Nind.com,length(com))))

        # 2 - Assign a species to each individual
        ### HYP 1: species abudances follow a lognormal distribution within each comm
        ### HYP 2: all species can be as abundant as each other within the regional pool (no rare/common species)
        ### HYP 3: species abundance is not linked to any type of gradient -> could we influence the distribution by both lognormal and gradient?
        for(c in 1:Ncom){
              ex.sp = sample(sp, size = Nind.com, prob = rlnorm(Nsp, 0, sdlog), replace = T)
              Data$sp[((c-1)*Nind.com+1):(c*Nind.com)] = ex.sp}

        # 3 - Assign a trait value to each individual
        ## OPTION 1: random ##
        if(Filter=='None'){
                  #trait1: normal distribution
                  ### Pourquoi pas simplement  rnorm(Nind, 0,1)? Ici le trait n'a pas une distribution normale?
                  #Data$trait1 <- rnorm(Nind, rlnorm(Nind, 0, 1), rlnorm(Nind, 0, 1))
                  Data$trait1 = rnorm(Nind, (max_value_traits-min_value_traits), (max_value_traits-min_value_traits)*cv_intra_sp)
                  #trait2: uniform distribution
                  Data$trait2 <- runif(Nind,min_value_traits,max_value_traits)
                  }

        ## OPTION 2: internal filtering ##
        if(Filter=='Internal'){

                  # Parameter for the distance between species mean trait values:
                  ### HYP 1: the most extreme case (if Int_Filter_Strength==100) species have equally distributed mean values along the trait gradient
                  Init_sp_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Nsp)), 2)

                  # Defining traits mean by species
                  mean.sp <- sample(rnorm(Nsp, mean = Init_sp_mean, sd = 100-Int_Filter_Strength), replace=F)
                  mean.sp[mean.sp<0] = runif(sum(mean.sp<0),min_value_traits, max_value_traits) ## to avoid negative values!!!

                  # Draw the individual traits values depending on species attributes
                  for(i in 1:Nind){
                        #trait 1 : normal distribution
                        Data$trait1[i] = rnorm(1, mean.sp[which(Data$sp[i] == sp)], mean.sp[which(Data$sp[i] == sp)]*cv_intra_sp)
                        #trait 2 : uniform distribution
                        Data$trait2[i] = runif(1, mean.sp[which(Data$sp[i] == sp)]*(1-cv_intra_sp),
                                           mean.sp[which(Data$sp[i] == sp)]*(1+cv_intra_sp))
                 }
        }

        ## OPTION 3: external filtering ##
        if(Filter=='External'){

                  # Parameter for the distance between communities mean trait values:
                  ### HYP 1: the most extreme case (if Ext_Filter_Strength==100) communities have equally distributed mean values along the trait gradient
                  Init_com_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Ncom)), 2)

                  # Defining traits mean by community
                  mean.com <- sample(rnorm(Ncom, mean = Init_com_mean, sd = 100-Ext_Filter_Strength), replace=F)
                  mean.com[mean.com<0] = runif(sum(mean.com<0),min_value_traits, max_value_traits) ## to avoid negative values!!!

                  # Draw the individual traits depending on communities attributes
                  for(i in 1:Nind){
                        #trait 1 : normal distribution
                        Data$trait1[i] = rnorm(1, mean.com[which(Data$com[i] == com)], mean.com[which(Data$com[i] == com)]*cv_intra_com)
                        #trait 2 : uniform distribution
                        Data$trait2[i] = runif(1, mean.com[which(Data$com[i] == com)]*(1-cv_intra_com),
                             mean.com[which(Data$com[i] == com)]*(1+cv_intra_com))
                  }
            }

        ## OPTION 4: external and internal filtering ##
        if(Filter=='Both'){

                  # Parameter for the distance between species mean trait values:
                  ### HYP 1: the most extreme case (if Int_Filter_Strength==100) species have equally distributed mean values along the trait gradient
                  Init_sp_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Nsp)), 2)

                  # Defining traits mean by species
                  mean.sp <- sample(rnorm(Nsp, mean = Init_sp_mean, sd = 100-Int_Filter_Strength), replace=F)
                  mean.sp[mean.sp<0] = runif(sum(mean.sp<0),min_value_traits, max_value_traits) ## to avoid negative values!!!

                  # Parameter for the distance between communities mean trait values:
                  ### HYP 1: the most extreme case (if Ext_Filter_Strength==100) communities have equally distributed mean values along the trait gradient
                  Init_com_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Ncom)), 2)

                  # Defining traits mean by community
                  mean.com <- sample(rnorm(Ncom, mean = Init_com_mean, sd = 100-Ext_Filter_Strength), replace=F)
                  mean.com[mean.com<0] = runif(sum(mean.com<0),min_value_traits, max_value_traits) ## to avoid negative values!!!

                  # Draw the individual traits values depending on species attributes
                  for(i in 1:Nind){
                        #trait 1 : normal distribution
                        Data$trait1[i] = 1/2*(rnorm(1, mean.sp[which(Data$sp[i] == sp)], mean.sp[which(Data$sp[i] == sp)]*cv_intra_sp) +
                                         rnorm(1, mean.com[which(Data$com[i] == com)], mean.com[which(Data$com[i] == com)]*cv_intra_com))# - mean(mean.com)
                        #trait 2 : uniform distribution
                        Data$trait2[i] = 1/2*(runif(1, mean.sp[which(Data$sp[i] == sp)]*(1-cv_intra_sp),
                                                  mean.sp[which(Data$sp[i] == sp)]*(1+cv_intra_sp)) +
                                         runif(1, mean.com[which(Data$com[i] == com)]*(1-cv_intra_com),
                                                  mean.com[which(Data$com[i] == com)]*(1+cv_intra_com)))# - mean(mean.com)
                 }
        }

             return(Data)
}

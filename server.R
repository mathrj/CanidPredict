# Define server logic to read selected file ----
####NOTE: In the server version, as.charac before as.numeric for mlr impute
shinyServer(function(input, output) {
  
    options(shiny.maxRequestSize = 30*1024^2)
    sample<- reactive({
        ext <- tools::file_ext(input$sample_up$name)
        switch(ext,
               csv = as.data.frame(read.csv(input$sample_up$datapath)),
               CGmap = make_methyl("dog_probes_canfam4.bed", input$sample_up$datapath),
               validate("Invalid file; Please upload a .csv or .cgmap file")
        )
    })
        output$sample_disp <- renderPrint({
        dimensions <- dim(sample()) 
        print(paste0("Your sample has ", dimensions[2], " methylation sites"))

        
    })
        output$sample_summary <- renderPrint({
          ts <- t(sample())
          ts <- ts[-1,]
          ts <-sapply(ts, as.numeric)
          print("Your Sample Average Methylation:")
          summary(ts)
        })
    
    # output$sample_hist <- renderPlot({
    #     validate(
    #         need(input$sample_up != "", "No data has been uploaded")
    #     )
    #     #show a histogram without the NA
    #     X2 <- as.matrix(sample())
    #     X3 <- X2[ , colSums(is.na(X2)) == 0]
    #     sample_plot <- ggplot(as.data.frame(X3), aes(X3))
    #     sample_plot +geom_bar(stat = "identity")
    #     
    # })
    
    output$age_sample  <- renderPlot({
        Sys.sleep(2)
        load("pred_age.RData")
        load("y_age.RData")
        load("x_methyl.RData")
        load("age_model.RData")
        
        d <-sample()

        #filter sample
        d <- d[,-1]
        remove_cols <- names(as.data.frame(x_methyl))

        new_sample = subset(d, select = (names(d) %in% remove_cols))
        
        no_na_df <- mlr::impute(new_sample)
        no_na_df <- no_na_df$data
        no_na_df[] <- lapply(no_na_df, as.numeric)

        #predict sample
        smpl_y <- round(predict(model1, newx = as.matrix(no_na_df)), digits = 2)
        plot.data <- data.frame(actual = y_age, predicted = pred_age)
        #abline(lm(predicted~actual,data=plot.data),col='red')
        gp <- ggplot(plot.data, aes(actual, predicted)) + geom_point() + xlab("Actual Age") + ylab("Predicted Age") + ggtitle("Actual Age vs Predicted Age")
        gp + geom_hline(aes(yintercept=smpl_y, linetype=paste("Predicted Age:", smpl_y)), color = "red") +
          geom_smooth(method='lm') +
          scale_linetype_manual(name = "", values = c(2), guide = guide_legend(override.aes = list(color = c("red"))))
    })
    
    
    
    output$sex <- renderPlot({
        load("pred_sex.RData")
        load("y_sex.RData")
        load("X_sex.RData")
        load("sex_model.RData")
        f_or_m <- ifelse(y_sex == 0, "MALE", "FEMALE")
        
        plot.sex <- data.frame(predicted_sex = pred_sex, actual_sex = f_or_m)
        validate(
            need(input$sample_up != "", "No data has been uploaded")
        )
        d <- sample()

        #filter sample
        remove_cols <- names(as.data.frame(X_sex))
        sample_na <- d[,-1]
        
        new_df = subset(sample_na, select = (names(sample_na) %in% remove_cols))
        X_sex <- as.matrix(X_sex)
        
        no_na_df <- mlr::impute(new_df)
        no_na_df <- no_na_df$data
        no_na_df[] <- lapply(no_na_df, as.numeric)
        
        sample_predict <- round(as.numeric(predict(fit, newx=as.matrix(no_na_df),s="lambda.min", type = "response")), digits = 3)
        sample_sex <- ifelse(sample_predict < 0.5, "MALE", "FEMALE")
        
        gs <- ggplot(plot.sex, aes(x=predicted_sex, y=actual_sex)) + geom_violin() + coord_flip() + xlab("Actual Sex") + ylab("Predicted Sex") + ggtitle("Actual Sex vs Predicted Sex")
        gs + geom_vline(aes(xintercept=sample_predict, linetype=paste("Your Sample Predicted Sex:", sample_predict, "(", sample_sex, ")")), color = "red") +
          scale_linetype_manual(name = "", values = c(2), guide = guide_legend(override.aes = list(color = c("red")))) 
    })
    
    output$wght <- renderPlot({
        load("pred_weight.RData")
        load("X_weight.RData")
        load("y_weight.RData")
        load("weight_model.RData")
        plot.weight <- data.frame(pred_weight,y_weight)
        validate(
            need(input$sample_up != "", "No data has been uploaded")
        )
        
        d <- sample()

        #filter sample
        remove_cols <- names(as.data.frame(X_weight))
        sample_na <- d[,-1]
        
        new_df = subset(sample_na, select = (names(sample_na) %in% remove_cols))
        
        no_na_df <- mlr::impute(new_df)
        no_na_df <- no_na_df$data
        no_na_df[] <- lapply(no_na_df, as.numeric)
        
        smpl_y <- round(predict(fit, newx = as.matrix(no_na_df)), digits = 2)
        plot.weight <- data.frame(pred_weight,y_weight)
        gw <- ggplot(plot.weight, aes(x= y_weight, y=pred_weight)) + geom_point() + xlab("Actual Value") + ylab("Predicted Value")
        gw + geom_hline(aes(yintercept=smpl_y, linetype=paste("Predicted weight: ", smpl_y)), color = "red") +
          scale_linetype_manual(name = "", values = c(2), guide = guide_legend(override.aes = list(color = c("red"))))
        
    })
    
    output$strst <- renderPlot({
      load("x_svm.RData") #filter_svm
      load("svmworks.RData")
      load("predvsobs_svm.RData")
      load("plotspay.RData")
      # load("predvsobs_svm.RData") #cm
      # Sys.setenv(Reticulate_Python="usr/bin/python3")
      
      p <- sample()
      sample_na <- p[,-1]
      # 
      # align colnames that were shifted
      for(i in seq(ncol(sample_na)-1)){
        sample_na[,i] <- sample_na[,i+1]
      }
      keep_cols <- names(as.data.frame(filter_svm))
      new_df = subset(sample_na, select = (names(sample_na) %in% keep_cols))
      
      no_na_df <- Hmisc::impute(new_df,mean)
      no_na_df<-sapply(no_na_df, as.character)
      
      no_na_df<-sapply(no_na_df, as.numeric)
      no_na_df <- t(as.data.frame(no_na_df))
      # # source_python('/Users/ponmathi/Desktop/DogAgePred/svm_onesample.py')
      # # 
      # # pred_sample <- svm_onesample(no_na_df)
      
      test_pred<-predict(svmfit,no_na_df)
      test_pred <- ifelse(test_pred == 1, "Spayed",
                          "Not Spayed")
      plot.spay$ypred <- ifelse(plot.spay$ypred == 1, "Spayed",
                                "Not Spayed")
      plot.spay$targets_female <- ifelse(plot.spay$targets_female == 1, "Spayed",
                                         "Not Spayed")
      map_data <- as.data.frame(table(plot.spay))
      
      colnames(map_data) <- c("Predicted.Spay", "Observed.Spay", "Count")
      
      gsp <- ggplot(data = map_data,mapping = aes(x = Predicted.Spay,y = Observed.Spay)) + geom_tile(aes(fill = Count)) +geom_text(aes(label = sprintf("%1.0f", Count)), vjust = 1) +
        scale_fill_gradient( low="#D6EAF8",high = "#2E86C1",trans = "log") +ggtitle(paste("Sample Prediction is", test_pred), subtitle = "Please note: Model uses only female samples (77 samples used)")
      
      gsp
      
    })
    
    output$phylot <- renderPlot({
        load("geno_file.RData")
        load("geno_names.RData")
        validate(
            need(input$geno_up != "", "No data has been uploaded")
        )
        genotypes_csv <- read.csv(input$geno_up$datapath)[2] #sample
        sample_geno <- as.matrix((genotypes_csv))
        
        #genotypesvar <- matrix(as.numeric(genotypesvar), ncol = ncol(genotypesvar), nrow = nrow(genotypesvar))
        #t_geno <- t(genotypesvar)
        #t_sample_geno <- t(sample_geno)
        full_geno_df <- cbind(sample_geno, genotypesvar)
        full_geno_df <- matrix(as.numeric(full_geno_df), ncol = ncol(full_geno_df), nrow = nrow(full_geno_df))
        
        full_geno_df[full_geno_df == 0] <- 3
        
        genotype_mat <- knncatimpute(full_geno_df)
        
        genotype_mat[genotype_mat == 3] <- 0
        
        genotype_df <- data.frame(genotype_mat)
        #colnames(genotype_df) <- colnames((genotypes_csv))[-c(1:2)]
        
        #find genetic distance
        full_names <- c("sample", geno_names)
        colnames(genotype_mat) <- full_names
        
        #transpose matrix 
        genotype_mat <- t(genotype_mat)
        
        gdist <- NAM::Gdist(genotype_mat, method = 1)

        hier_cluster_ml <- hclust(gdist, method = 'ward.D')
        
        d <- as.phylo(hier_cluster_ml)
        bi_subset <- tree_subset(d, "sample", levels_back = input$levels)
        
        p <- bi_subset %>% 
            ggtree(aes(color = group)) + 
            geom_tiplab() + 
            theme_tree2() + 
            scale_color_manual(values = c(`1` = "red", `0` = "black")) +
            xlim(0, 0.75)
        p + lims(x = c(0, max(p$data$x) * input$tree_width))
    })
    
    
    output$behav <- renderPlot({
      #four plots in one?
      load("atten_pred.RData")
      load("atten_model.RData")
      load("social_model.RData")
      load("nonsocial_model.RData")
      load("energy_model.RData")
      load("n_fear_plot.RData")
      load("energy_plot.RData")
      load("s_fear_plot.RData")
      
      load("atten_y.RData")
      load("all_behav.RData")
      
      
      #repeat for each behaviour
      
      #for one sample
      s <- sample()
      #
      sample_na <- s[,-1]
      sample_na <- sample_na[,-1]
      remove_cols <- names(as.data.frame(all_behav))
      new_df = subset(sample_na, select = (names(sample_na) %in% remove_cols))
      no_na_df <- mlr::impute(new_df)
      no_na_df <- no_na_df$data
      no_na_df[] <- lapply(no_na_df, as.character)
      no_na_df[] <- lapply(no_na_df, as.numeric)
      #issue is that cgmap turns it into character
      #no_na_df <- as.data.frame(lapply(no_na_df,as.numeric))
      
      a_pred <- round(as.numeric(predict(atten_model, no_na_df, ncomp = 2)),digits=3)
      
      e_pred <- round(as.numeric(predict(energy_model, no_na_df, ncomp = 2)),digits=3)
      
      s_pred <- round(as.numeric(predict(social_model, no_na_df, ncomp = 2)),digits=3)
      
      n_pred <- round(as.numeric(predict(nonsocial_model, no_na_df, ncomp = 2)),digits=3)
      
      # energy_plot <- data.frame(energy_plot$energy_pred,energy_plot$energy_y) 
      attention_plot <- data.frame(atten_pred,atten_y)
      # n_fear_plot <- data.frame(n_fear_plot$n_fear_pred,n_fear_plot$n_fear_y)
      # s_fear_plot <- data.frame(s_fear_plot$s_fear_pred,s_fear_plot$s_fear_y)
      
      n <- 1
      p1 <- ggplot(energy_plot, aes(x = energy_y,y = energy_pred)) + geom_point() +
        ggtitle(paste("Actual vs Predicted Energy for all the Dogs in our Database (Corr coeff = ", round(cor(energy_plot$energy_pred,energy_plot$energy_y), digits = 3), " )"), subtitle = "Energy level: Energetic, “always on the go”, and/or playful.")+
        xlab("Actual Energy ") + ylab("Predicted Energy")+
        geom_hline(aes(yintercept=e_pred, linetype=paste("Predicted Sample Energy: ", e_pred, "(",if_else(e_pred > 2, "high", "low"), ")")), color = "red") +
        scale_linetype_manual(name = "", values = c(rep(2,44)),  guide = guide_legend(override.aes = list(color = c("red")))) +
        coord_fixed(ratio=n)+ theme(aspect.ratio=1) +
        geom_smooth(method='lm') + scale_x_continuous(breaks=c(0,2, 4), labels=c("Low", "Mid", "High"))
      
      p2 <- ggplot(attention_plot, aes(atten_y, atten_pred)) + geom_point()+
        ggtitle(paste("Actual vs Predicted Attention for all the Dogs in our Database (Corr coeff = ", round(cor(attention_plot$atten_pred,attention_plot$atten_y), digits =3), " )"), subtitle = "Attachment and attention-seeking: Maintaining close proximity to the owner or other members of the household, soliciting affection or attention, and displaying agitation when the owner gives attention to third parties.") + 
        xlab("Actual Attention ") + ylab("Predicted Attention")+
        geom_hline(aes(yintercept=a_pred, linetype=paste("Predicted Sample Attention: ", a_pred, "(",if_else(a_pred > 2, "high", "low"), ")")), color = "red") + 
        scale_linetype_manual(name = "",values = c(rep(2,44)), guide = guide_legend(override.aes = list(color = c("red")))) +
        coord_fixed(ratio=n)+ theme(aspect.ratio=1)+
        geom_smooth(method='lm')  + scale_x_continuous(breaks=c(0,2, 4), labels=c("Low", "Mid", "High"))
      # p3 <- ggplot(n_fear_plot, aes(n_fear_y,n_fear_pred)) + geom_point()+
      #     ggtitle(paste("Actual vs Predicted Nonsocial fear for all the Dogs in our Database (Corr coeff = ", round(cor(n_fear_plot$n_fear_pred,n_fear_plot$n_fear_y),digits = 3), " )")) +
      #     xlab("Actual Nonsocial fear") + ylab("Predicted Nonsocial fear")+
      #     geom_hline(aes(yintercept=n_pred, linetype=paste("Predicted Sample Nonsocial fear: ", n_pred)), color = "red") +
      #     scale_linetype_manual(name = "", values = c(rep(2,44)),  guide = guide_legend(override.aes = list(color = c("red")))) +
      #     coord_fixed(ratio=n)+ theme(aspect.ratio=1)+
      #     geom_smooth(method='lm') 
      
      p4 <- ggplot(s_fear_plot, aes(s_fear_y, s_fear_pred)) + geom_point()+ 
        ggtitle(paste("Actual vs Predicted Stranger-directed fear for all the Dogs in our Database (Coff coeff = ", round(cor(s_fear_plot$s_fear_pred,s_fear_plot$s_fear_y), digits=3), " )"), subtitle = "Stranger-directed fear: Fearful or wary responses when approached by strangers.")+
        xlab("Actual Stranger-directed fear ") + ylab("Predicted Stranger-directed fear")+
        geom_hline(aes(yintercept=s_pred, linetype=paste("Predicted Sample Stranger-directed fear: ", s_pred, "(",if_else(s_pred > 2, "high", "low"), ")")), color = "red") +
        scale_linetype_manual(name = "", values = c(rep(2,44)), guide = guide_legend(override.aes = list(color = c("red")))) +
        coord_fixed(ratio=n)+ theme(aspect.ratio=1)+
        geom_smooth(method='lm')   + scale_x_continuous(breaks=c(0,2, 4), labels=c("Low", "Mid", "High"))
      
      
      
      grid.arrange(p1,p2,p4, ncol=1)
      
    }, height = 1500, width = 1000)

})



